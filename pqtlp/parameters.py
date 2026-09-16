"""The single concrete parameter set used by the implementation."""

from sage.all import ZZ, next_prime

SECURITY_BITS = 128
SOLUTION_BYTES = SECURITY_BITS // 8
WALK_SEED_BYTES = 32

# The supersingular base field and its available power-of-two torsion.
P = ZZ(3) * ZZ(2) ** 324 - 1
FULL_TORSION_EXPONENT = 324

# The vertical isogeny has degree N = q(2^a - q).
PSI_TORSION_EXPONENT = 257
Q = ZZ(next_prime(ZZ(2) ** 256))
VERTICAL_DEGREE = Q * (ZZ(2) ** PSI_TORSION_EXPONENT - Q)

# A chunk must fit both the 2^a representation and its complementary torsion.
CHUNK_SIZE = min(
    PSI_TORSION_EXPONENT,
    FULL_TORSION_EXPONENT - PSI_TORSION_EXPONENT,
)

SHORTCUT_SEARCH_ATTEMPTS = 100
