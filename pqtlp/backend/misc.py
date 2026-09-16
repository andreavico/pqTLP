"""Fixed-width encodings used by the domain-separated hashes."""

from . import params


def encode_fp(x):
    """Encode an Fp element in little-endian, 64-bit-word-padded form."""
    return int(x).to_bytes(params.FP_ENC_BYTES, "little")


def encode_fp2(x):
    """Encode a + b*i as the concatenation of the two Fp encodings."""
    a, b = x
    return encode_fp(a) + encode_fp(b)


def encode_curve_j(E):
    """Encode the curve's j-invariant in Fp2."""
    return encode_fp2(E.j_invariant())
