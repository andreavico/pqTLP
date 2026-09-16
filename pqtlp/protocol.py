"""Setup, base puzzles, and the paper's flat sequential compiler."""

import multiprocessing
import os
import pickle
from typing import Optional

from sage.all import ZZ, gcd, set_random_seed

from .backend import params
from .crypto import mask_next_puzzle, mask_solution_with_curve, unmask_next_puzzle
from .curves import (
    check_kernel_order_2e,
    deterministic_basis_2e,
    evaluate_2d_on_torsion,
    ideal_to_2d_repr,
    isogeny_from_kernel_2e_x_only,
    kernel_from_basis_2e,
    montgomery_isogeny_from_kernel,
    point_coords_in_basis_2e,
    push_2d_through_cgl_chunk,
    seeded_cgl_scalar,
)
from .delay import build_seeded_delay
from .ideals import (
    ideal_side_pushforward_codomain,
    sample_psi_ideal,
    sample_random_base_curve,
    short_equivalent_ideal,
)
from .parameters import (
    CHUNK_SIZE,
    SHORTCUT_SEARCH_ATTEMPTS,
    SOLUTION_BYTES,
    VERTICAL_DEGREE,
    WALK_SEED_BYTES,
)
from .types import BasePuzzle, CompiledPuzzle, Isogeny2D, SetupState

_PGEN_WORKER_STATE = None


def _initialise_pgen_worker(setup: SetupState) -> None:
    """Install shared setup state and independently reseed a forked worker."""

    global _PGEN_WORKER_STATE
    _PGEN_WORKER_STATE = setup
    set_random_seed(int.from_bytes(os.urandom(32), "big"))


def _pgen_base_worker(solution: bytes) -> bytes:
    if _PGEN_WORKER_STATE is None:
        raise RuntimeError("PGen worker has no setup state")
    # Sage/PARI finite-field objects are not safe to unpickle in the pool's
    # background result-handler thread.  Return an opaque byte string and
    # deserialize it on the calling thread below.
    return pickle.dumps(
        pgen_base(_PGEN_WORKER_STATE, solution),
        protocol=pickle.HIGHEST_PROTOCOL,
    )


def generate_base_puzzles(
    setup: SetupState,
    solutions: list[bytes],
    *,
    workers: int = 1,
) -> list[BasePuzzle]:
    """Generate independent base puzzles serially or in forked workers."""

    workers = int(workers)
    if workers <= 0:
        raise ValueError("workers must be positive")
    if workers == 1 or len(solutions) <= 1:
        return [pgen_base(setup, solution) for solution in solutions]
    worker_count = min(workers, len(solutions))
    # Base-puzzle generation is independent after setup.  Forking keeps the
    # large immutable Sage parameter objects copy-on-write and avoids paying a
    # fresh Sage import in every worker.
    context = multiprocessing.get_context("fork")
    with context.Pool(
        worker_count,
        initializer=_initialise_pgen_worker,
        initargs=(setup,),
    ) as pool:
        encoded = pool.map(_pgen_base_worker, solutions, chunksize=1)
    return [pickle.loads(item) for item in encoded]


def setup_delay(
    delay_length: int,
    *,
    walk_seed: Optional[bytes] = None,
) -> SetupState:
    """Set up a reusable base puzzle for ``delay_length`` degree-two steps."""

    delay_length = int(delay_length)
    if delay_length <= 0 or delay_length % CHUNK_SIZE != 0:
        raise ValueError(
            "delay_length must be a positive multiple of the pushforward "
            f"chunk size {CHUNK_SIZE}"
        )
    if walk_seed is None:
        walk_seed = os.urandom(WALK_SEED_BYTES)
    walk_seed = bytes(walk_seed)
    if len(walk_seed) != WALK_SEED_BYTES:
        raise ValueError(f"walk_seed must be exactly {WALK_SEED_BYTES} bytes")

    base_ideal, base_curve, base_transport_basis = sample_random_base_curve()
    base_full_basis = deterministic_basis_2e(base_curve, int(params.f))
    base_change_matrix = tuple(
        point_coords_in_basis_2e(
            X,
            base_transport_basis[0],
            base_transport_basis[1],
            int(params.f),
        )
        for X in base_full_basis
    )
    if gcd(ZZ(base_ideal.norm()), ZZ(2) * VERTICAL_DEGREE) != 1:
        raise RuntimeError("randomized base ideal has an incompatible norm")

    delay_ideal = build_seeded_delay(
        walk_seed,
        delay_length,
        base_ideal=base_ideal,
        base_curve=base_curve,
        base_transport_basis=base_transport_basis,
        chunk_size=CHUNK_SIZE,
        vertical_degree=VERTICAL_DEGREE,
    )
    shortcut_ideal = short_equivalent_ideal(
        delay_ideal,
        congruence_modulus=VERTICAL_DEGREE,
        degree_coprime_to=ZZ(2) * VERTICAL_DEGREE,
        attempts=SHORTCUT_SEARCH_ATTEMPTS,
    )
    shortcut_push_lattice = shortcut_ideal.right_order() * VERTICAL_DEGREE
    return SetupState(
        delay_seed=walk_seed,
        delay_length=delay_length,
        base_curve=base_curve,
        base_full_basis=base_full_basis,
        base_change_matrix=base_change_matrix,
        base_ideal=base_ideal,
        shortcut_ideal=shortcut_ideal,
        shortcut_push_lattice=shortcut_push_lattice,
    )


def public_pushforward_codomain(
    delay_seed: bytes,
    delay_length: int,
    psi: Isogeny2D,
):
    """Push the public 2D representation through the seeded CGL walk."""

    delay_seed = bytes(delay_seed)
    delay_length = int(delay_length)
    if len(delay_seed) != WALK_SEED_BYTES:
        raise ValueError(f"walk seed must be exactly {WALK_SEED_BYTES} bytes")
    if delay_length <= 0 or delay_length % CHUNK_SIZE != 0:
        raise ValueError("invalid delay length")

    chunk_count = delay_length // CHUNK_SIZE
    sigma = psi
    forbidden_backtrack = None
    for i in range(chunk_count):
        scalar, domain_basis = seeded_cgl_scalar(
            delay_seed,
            sigma.domain,
            i,
            CHUNK_SIZE,
            forbidden_backtrack,
        )
        if i == chunk_count - 1:
            domain_kernel = kernel_from_basis_2e(domain_basis, scalar)
            check_kernel_order_2e(
                domain_kernel,
                CHUNK_SIZE,
                label="final CGL chunk kernel",
            )
            kernel = evaluate_2d_on_torsion(sigma, domain_kernel)
            sigma.codomain.set_order((params.p + 1) ** 2)
            try:
                # The Kummer route is faster, but Sage may reject some models;
                # the fallback below keeps the public algorithm robust.
                _, output_curve = isogeny_from_kernel_2e_x_only(
                    sigma.codomain,
                    kernel,
                    CHUNK_SIZE,
                )
                return output_curve
            except (ValueError, NotImplementedError, AssertionError):
                return montgomery_isogeny_from_kernel(
                    sigma.codomain, kernel, CHUNK_SIZE
                ).codomain()

        sigma, forbidden_backtrack = push_2d_through_cgl_chunk(
            sigma,
            scalar,
            CHUNK_SIZE,
            basis=domain_basis,
        )

    raise RuntimeError("unreachable")


def _solution_bytes(message: bytes) -> bytes:
    """Normalize and validate a fixed-length solution before doing setup."""
    if isinstance(message, str):
        message = message.encode()
    message = bytes(message)
    if len(message) != SOLUTION_BYTES:
        raise ValueError(f"message must be exactly {SOLUTION_BYTES} bytes")
    return message


def pgen_base(setup: SetupState, message: bytes) -> BasePuzzle:
    """Generate one base puzzle from reusable setup state."""

    message = _solution_bytes(message)

    psi_ideal = sample_psi_ideal(base_ideal=setup.base_ideal)
    if gcd(VERTICAL_DEGREE, ZZ(setup.shortcut_ideal.norm())) != 1:
        raise ValueError("vertical degree is not coprime to the shortcut norm")

    key_curve = ideal_side_pushforward_codomain(
        setup.shortcut_ideal,
        psi_ideal,
        base_ideal=setup.base_ideal,
        shortcut_push_lattice=setup.shortcut_push_lattice,
    )
    psi_representation = ideal_to_2d_repr(
        psi_ideal,
        base_curve=setup.base_curve,
        base_full_basis=setup.base_full_basis,
        base_change_matrix=setup.base_change_matrix,
        base_ideal=setup.base_ideal,
    )
    masked_solution = mask_solution_with_curve(key_curve, message)
    return BasePuzzle(
        delay_seed=setup.delay_seed,
        delay_length=setup.delay_length,
        psi=psi_representation,
        masked_solution=masked_solution,
    )


def psolve_base(puzzle: BasePuzzle) -> bytes:
    """Solve one base puzzle using only its public data."""

    E_key = public_pushforward_codomain(
        puzzle.delay_seed,
        puzzle.delay_length,
        puzzle.psi,
    )
    return mask_solution_with_curve(E_key, puzzle.masked_solution)


def pgen(
    root_delay: int,
    message: bytes,
    *,
    tlp_count: Optional[int] = None,
    workers: int = 1,
) -> CompiledPuzzle:
    """Generate the flat compiler sequence ``(Z_1, ct_2, ..., ct_T1)``."""

    message = _solution_bytes(message)
    workers = int(workers)
    if workers <= 0:
        raise ValueError("workers must be positive")
    root_delay = int(root_delay)
    if root_delay <= 0:
        raise ValueError("root_delay must be positive")
    if tlp_count is None:
        tlp_count = root_delay
    tlp_count = int(tlp_count)
    if tlp_count <= 0:
        raise ValueError("tlp_count must be positive")

    if root_delay % CHUNK_SIZE != 0:
        raise ValueError(
            f"root_delay must be a multiple of the pushforward chunk size {CHUNK_SIZE}"
        )

    setup = setup_delay(root_delay)

    intermediate_keys = [os.urandom(SOLUTION_BYTES) for _ in range(tlp_count - 1)]
    solutions = [*intermediate_keys, message]
    base_puzzles = generate_base_puzzles(setup, solutions, workers=workers)

    # Paper Section 4.2: Z* = (Z_1, ct_2, ..., ct_T1), where each
    # ct_{i+1} = encode(Z_{i+1}) XOR H(k_i).
    masked_puzzles = []
    for key, next_puzzle in zip(intermediate_keys, base_puzzles[1:]):
        encoded = pickle.dumps(next_puzzle, protocol=pickle.HIGHEST_PROTOCOL)
        masked_puzzles.append(mask_next_puzzle(key, encoded))

    return CompiledPuzzle(
        first_puzzle=base_puzzles[0],
        masked_puzzles=masked_puzzles,
        root_delay=root_delay,
        tlp_count=tlp_count,
    )


def psolve(puzzle: CompiledPuzzle) -> bytes:
    """Solve each base puzzle in sequence and return the final message."""

    if len(puzzle.masked_puzzles) != puzzle.tlp_count - 1:
        raise ValueError("compiled puzzle has the wrong number of masked puzzles")

    current = puzzle.first_puzzle
    for masked in puzzle.masked_puzzles:
        key = psolve_base(current)
        encoded = unmask_next_puzzle(key, masked)
        current = pickle.loads(encoded)
        if not isinstance(current, BasePuzzle):
            raise ValueError("masked compiler entry is not a base puzzle")
    return psolve_base(current)
