"""Secret setup state and serializable public puzzle representations."""

from dataclasses import dataclass

from .parameters import FULL_TORSION_EXPONENT, PSI_TORSION_EXPONENT


@dataclass
class Isogeny2D:
    """Two-dimensional representation of an odd isogeny ``psi: E -> F``.

    The public puzzle carries the action of ``psi`` on the ``2^a`` torsion.
    This stores a full ``2^f`` basis of the domain only so the implementation
    can derive the smaller ``2^a`` basis used by the representation.
    """

    domain: object
    codomain: object
    full_domain_basis: tuple
    images: tuple

    @property
    def P_image(self):
        return self.images[0]

    @property
    def Q_image(self):
        return self.images[1]

    @property
    def torsion_basis(self):
        scale = 2 ** (FULL_TORSION_EXPONENT - PSI_TORSION_EXPONENT)
        P, Q = self.full_domain_basis
        return scale * P, scale * Q


@dataclass
class SetupState:
    """Secret setup state for one reusable base TLP.

    The public puzzle only needs the seeded delay description and a 2D
    representation of ``psi``.  ``base_ideal`` witnesses the endomorphism ring
    of the freshly randomized public base curve.  Only the witness, its
    base-change matrix, and the final shortcut are retained after setup; the
    per-chunk ideals are released as setup advances.
    """

    delay_seed: bytes
    delay_length: int
    base_curve: object
    base_full_basis: tuple
    base_change_matrix: tuple
    base_ideal: object
    shortcut_ideal: object
    shortcut_push_lattice: object


@dataclass
class BasePuzzle:
    """Public base TLP instance generated from ``SetupState``.

    The solver recomputes the CGL delay from ``delay_seed`` and pushes the
    sampled isogeny ``psi`` through that delay.  The final curve's j-invariant
    determines the random-oracle mask used to recover ``masked_solution``.
    """

    delay_seed: bytes
    delay_length: int
    psi: Isogeny2D
    masked_solution: bytes


@dataclass
class CompiledPuzzle:
    """Plain-model puzzle built with the paper's flat compiler.

    ``first_puzzle`` is solved first.  Each entry of ``masked_puzzles`` is the
    encoding of the next base puzzle XORed with a hash of the preceding
    solution, matching the flat compiler in the paper.  The paper's balanced
    setting uses ``tlp_count == root_delay``.
    """

    first_puzzle: BasePuzzle
    masked_puzzles: list[bytes]
    root_delay: int
    tlp_count: int


@dataclass(frozen=True)
class DelayCalibrationChoice:
    """Estimated parameters for a requested wall-clock solve delay."""

    target_seconds: float
    root_delay: int
    tlp_count: int
    predicted_solve_seconds: float
    predicted_pgen_seconds: float
    base_solve_seconds: float
    setup_seconds: float
    tlp_pgen_seconds: float
