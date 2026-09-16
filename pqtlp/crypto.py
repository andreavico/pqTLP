"""Domain-separated SHAKE-256 masks for solutions and compiler entries."""

import hashlib

from .backend import misc
from .parameters import SOLUTION_BYTES


def _xor_bytes(a: bytes, b: bytes) -> bytes:
    if len(a) != len(b):
        raise ValueError("XOR inputs must have equal length")
    return bytes(x ^ y for x, y in zip(a, b))


def hpf_curve(E, output_length: int = SOLUTION_BYTES) -> bytes:
    """Instantiate the paper's ``H_pf(j(E))`` with SHAKE-256."""

    output_length = int(output_length)
    if output_length < 0:
        raise ValueError("output_length must be nonnegative")
    return hashlib.shake_256(b"pqtlp-hpf-v1" + misc.encode_curve_j(E)).digest(
        output_length
    )


def mask_solution_with_curve(E, solution: bytes) -> bytes:
    """Return ``H_pf(j(E)) XOR solution`` for a lambda-bit solution."""

    solution = bytes(solution)
    if len(solution) != SOLUTION_BYTES:
        raise ValueError(f"solution must be exactly {SOLUTION_BYTES} bytes")
    return _xor_bytes(solution, hpf_curve(E))


def compiler_mask(key: bytes, length: int) -> bytes:
    """Expand the paper compiler's random-oracle value ``H(key)``."""

    key = bytes(key)
    length = int(length)
    if len(key) != SOLUTION_BYTES:
        raise ValueError(f"compiler key must be exactly {SOLUTION_BYTES} bytes")
    if length < 0:
        raise ValueError("length must be nonnegative")
    return hashlib.shake_256(b"pqtlp-compiler-v1" + key).digest(length)


def mask_next_puzzle(key: bytes, encoded_puzzle: bytes) -> bytes:
    """Return the flat compiler ciphertext ``encode(Z) XOR H(key)``."""

    encoded_puzzle = bytes(encoded_puzzle)
    return _xor_bytes(encoded_puzzle, compiler_mask(key, len(encoded_puzzle)))


def unmask_next_puzzle(key: bytes, masked_puzzle: bytes) -> bytes:
    """Invert ``mask_next_puzzle``."""

    return mask_next_puzzle(key, masked_puzzle)
