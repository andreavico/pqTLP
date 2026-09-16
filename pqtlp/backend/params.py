"""Initialise the fixed Sage objects shared by the mathematical backend."""

from sage.all import GF, EllipticCurve, QuaternionAlgebra, Zmod, matrix, next_prime

from ..parameters import FULL_TORSION_EXPONENT, P
from . import constants

p = P
f = FULL_TORSION_EXPONENT
D_mix = next_prime(p**2)
# Keep field encodings padded to whole 64-bit words.
FP_ENC_BYTES = 8 * ((p.nbits() + 63) // 64)

# Constants used by the Qlapoti ideal-to-isogeny backend.
QUAT_prime_cofactor = next_prime(2 ** p.nbits())

B = O0 = None
Fp = Fp2 = Fp2_i = None
E0 = P0 = Q0 = None
mat_1 = mat_i = mat_ij2 = mat_1k2 = None


def initialize() -> None:
    """Construct the field, reference curve, order, and torsion basis once."""

    global B, O0, Fp, Fp2, Fp2_i, E0, P0, Q0
    global mat_1, mat_i, mat_ij2, mat_1k2

    if E0 is not None:
        return

    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.maximal_order(order_basis=(B(1), i, (i + j) / 2, (1 - k) / 2))

    Fp = GF(p)
    Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p + 1) ** 2)

    P0 = E0(constants.P_X, constants.P_Y)
    Q0 = E0(constants.Q_X, constants.Q_Y)
    ring = Zmod(2**f)
    mat_1 = matrix(ring, 2, [1, 0, 0, 1])
    mat_i = matrix(ring, 2, constants.MAT_I)
    mat_ij2 = matrix(ring, 2, constants.MAT_IJ2)
    mat_1k2 = matrix(ring, 2, constants.MAT_1K2)
