"""SageMath implementation of the post-quantum time-lock puzzle."""

from .backend import params as _backend_params
from .calibration import DELAY_PROFILE, choose_parameters, estimate_times
from .parameters import (
    CHUNK_SIZE,
    FULL_TORSION_EXPONENT,
    PSI_TORSION_EXPONENT,
    SECURITY_BITS,
    SOLUTION_BYTES,
    VERTICAL_DEGREE,
    WALK_SEED_BYTES,
    P,
    Q,
)
from .protocol import (
    generate_base_puzzles,
    pgen,
    pgen_base,
    psolve,
    psolve_base,
    setup_delay,
)
from .types import BasePuzzle, CompiledPuzzle, DelayCalibrationChoice, Isogeny2D

_backend_params.initialize()

__all__ = [
    "BasePuzzle",
    "CHUNK_SIZE",
    "CompiledPuzzle",
    "DELAY_PROFILE",
    "DelayCalibrationChoice",
    "FULL_TORSION_EXPONENT",
    "Isogeny2D",
    "P",
    "PSI_TORSION_EXPONENT",
    "Q",
    "SECURITY_BITS",
    "SOLUTION_BYTES",
    "VERTICAL_DEGREE",
    "WALK_SEED_BYTES",
    "choose_parameters",
    "estimate_times",
    "generate_base_puzzles",
    "pgen",
    "pgen_base",
    "psolve",
    "psolve_base",
    "setup_delay",
]
