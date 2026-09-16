"""Setup-side construction of the compact seeded CGL delay."""

from sage.all import ZZ

from .backend import params
from .backend import qlapoti as qlpt
from .backend import quaternions as qt
from .curves import (
    cgl_chunk_isogeny,
    cgl_dual_backtrack_subgroup,
    kernel_decomposed_to_ideal_2e,
    point_coords_in_basis_2e,
    seeded_cgl_scalar,
)
from .ideals import short_equivalent_ideal
from .parameters import WALK_SEED_BYTES


def transport_basis(
    composite_ideal,
    target_curve,
    exponent: int,
):
    """Transport the reference torsion basis to ``target_curve``."""

    transported_curve, P_full, Q_full = qlpt.IdealToIsogeny(composite_ideal)
    if not transported_curve.is_isomorphic(target_curve):
        raise RuntimeError("transport ideal has the wrong codomain")

    scale = ZZ(2) ** ZZ(params.f - int(exponent))
    isomorphism = transported_curve.isomorphism_to(target_curve)
    return isomorphism(scale * P_full), isomorphism(scale * Q_full)


def build_seeded_delay(
    seed: bytes,
    delay_length: int,
    *,
    base_ideal,
    base_curve,
    base_transport_basis,
    chunk_size: int,
    vertical_degree: ZZ,
):
    """Build the setup-side ideal for a seeded non-backtracking CGL walk.

    Each chunk is checked as it is constructed. Only the current prefix ideal
    is retained, so setup does not accumulate a transport transcript.
    """

    seed = bytes(seed)
    delay_length = int(delay_length)
    chunk_size = int(chunk_size)
    vertical_degree = ZZ(vertical_degree)
    if len(seed) != WALK_SEED_BYTES:
        raise ValueError(f"walk seed must be exactly {WALK_SEED_BYTES} bytes")
    if not 1 <= chunk_size <= params.f:
        raise ValueError(f"chunk_size must satisfy 1 <= chunk_size <= {params.f}")
    if delay_length <= 0 or delay_length % chunk_size != 0:
        raise ValueError("delay length must be a positive multiple of chunk_size")

    scale = ZZ(2) ** ZZ(params.f - chunk_size)
    base_chunk_basis = tuple(scale * point for point in base_transport_basis)
    curve = base_curve
    prefix_ideal = None
    forbidden_backtrack = None

    for chunk_index in range(delay_length // chunk_size):
        scalar, public_basis = seeded_cgl_scalar(
            seed,
            curve,
            chunk_index,
            chunk_size,
            forbidden_backtrack,
        )

        if prefix_ideal is None:
            transported_basis = base_chunk_basis
            short_prefix = None
        else:
            short_prefix = short_equivalent_ideal(
                prefix_ideal,
                congruence_modulus=vertical_degree,
                degree_coprime_to=ZZ(2) * vertical_degree,
            )
            composite_ideal = base_ideal * short_prefix
            transported_basis = transport_basis(
                composite_ideal,
                curve,
                chunk_size,
            )

        isogeny, _, public_kernel = cgl_chunk_isogeny(
            curve,
            scalar,
            chunk_size,
            basis=public_basis,
        )
        c1, c2 = point_coords_in_basis_2e(
            public_kernel,
            transported_basis[0],
            transported_basis[1],
            chunk_size,
        )
        reference_ideal = kernel_decomposed_to_ideal_2e(c1, c2, chunk_size)

        if short_prefix is None:
            chunk_ideal = qt.Pushforward(reference_ideal, base_ideal)
            prefix_ideal = chunk_ideal
        else:
            # Pushforward through the already-computed composite in one step.
            chunk_ideal = qt.Pushforward(reference_ideal, composite_ideal)
            prefix_ideal = short_prefix * chunk_ideal

        expected_norm = ZZ(2) ** chunk_size
        if ZZ(chunk_ideal.norm()) != expected_norm:
            raise RuntimeError("delay chunk ideal has the wrong norm")
        if prefix_ideal.left_order() != base_ideal.right_order():
            raise RuntimeError("delay prefix escaped the randomized base order")

        forbidden_backtrack = cgl_dual_backtrack_subgroup(
            isogeny,
            public_basis,
            chunk_size,
        )
        curve = isogeny.codomain()

    return prefix_ideal
