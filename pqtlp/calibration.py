"""Wall-clock estimates from the recorded full-parameter measurements."""

import math
from typing import Optional

from .types import DelayCalibrationChoice

DELAY_PROFILE = {
    # SageMath 10.8 timings on an Apple M3 for lambda=128,
    # p=3*2^324-1, q=next_prime(2^256), a=257, and 67-step chunks.
    # Entries are medians of three independent full round trips, except for
    # the 4,288-step entry, which is the median of two long validation runs.
    67: {"setup": 2.681016, "tlp_pgen": 1.379246, "base_solve": 0.677777},
    134: {"setup": 7.695152, "tlp_pgen": 1.264468, "base_solve": 3.445756},
    268: {"setup": 18.648962, "tlp_pgen": 1.203975, "base_solve": 11.330542},
    536: {"setup": 40.339866, "tlp_pgen": 1.138889, "base_solve": 18.930294},
    1072: {"setup": 68.336084, "tlp_pgen": 1.125961, "base_solve": 39.172114},
    2144: {"setup": 129.113033, "tlp_pgen": 1.046129, "base_solve": 77.956241},
    4288: {"setup": 259.630579, "tlp_pgen": 1.205141, "base_solve": 164.240522},
}


def estimate_times(
    root_delay: int,
    tlp_count: int,
    profile: Optional[dict] = None,
) -> tuple[float, float]:
    """Estimate generation and solving time for a compiled puzzle.

    ``root_delay`` is the number of degree-two CGL steps. ``tlp_count`` is
    the number of base TLPs that will be solved sequentially.  The model is a
    simple empirical fit: setup is paid once, generation adds one base-puzzle
    cost per TLP, and solving adds one base-solve cost per TLP.
    """

    if profile is None:
        profile = DELAY_PROFILE
    root_delay = int(root_delay)
    tlp_count = int(tlp_count)
    if tlp_count <= 0:
        raise ValueError("tlp_count must be positive")
    if root_delay not in profile:
        raise ValueError(f"no calibration profile for root_delay={root_delay}")
    entry = profile[root_delay]
    pgen = float(entry["setup"]) + tlp_count * float(entry["tlp_pgen"])
    psolve = tlp_count * float(entry["base_solve"])
    return pgen, psolve


def choose_parameters(
    target_seconds: float,
    *,
    profile: Optional[dict] = None,
    min_solve_ratio: float = 1.0,
) -> DelayCalibrationChoice:
    """Choose calibrated ``root_delay`` and ``tlp_count`` values.

    The chooser minimizes estimated puzzle-generation time among calibrated
    candidates whose predicted solving time is at least
    ``min_solve_ratio * target_seconds``. Returned values estimate local
    running time; calibration noise can make an actual solve faster or slower.
    """

    if profile is None:
        profile = DELAY_PROFILE
    target_seconds = float(target_seconds)
    if not math.isfinite(target_seconds) or target_seconds <= 0:
        raise ValueError("target_seconds must be finite and positive")
    min_solve_ratio = float(min_solve_ratio)
    if not math.isfinite(min_solve_ratio) or min_solve_ratio <= 0:
        raise ValueError("min_solve_ratio must be finite and positive")

    candidates = []
    lower_bound = min_solve_ratio * target_seconds
    if not math.isfinite(lower_bound):
        raise ValueError("requested solve time is too large")
    for root_delay, entry in sorted(profile.items()):
        base_solve = float(entry["base_solve"])
        setup = float(entry["setup"])
        base_pgen = float(entry["tlp_pgen"])
        if any(not math.isfinite(v) or v <= 0 for v in (base_solve, setup, base_pgen)):
            raise ValueError("calibration times must be finite and positive")
        # Generation cost increases with the count, so only the smallest
        # count meeting the requested delay can be optimal for this row.
        tlp_count = max(1, math.ceil(lower_bound / base_solve))
        if tlp_count * base_solve < lower_bound:
            tlp_count += 1
        pgen, psolve = estimate_times(root_delay, tlp_count, profile=profile)
        candidates.append(
            DelayCalibrationChoice(
                target_seconds=target_seconds,
                root_delay=int(root_delay),
                tlp_count=tlp_count,
                predicted_solve_seconds=psolve,
                predicted_pgen_seconds=pgen,
                base_solve_seconds=base_solve,
                setup_seconds=setup,
                tlp_pgen_seconds=base_pgen,
            )
        )

    if not candidates:
        raise RuntimeError("no calibrated delay candidate found")
    return min(
        candidates,
        key=lambda c: (
            c.predicted_pgen_seconds,
            abs(c.predicted_solve_seconds - c.target_seconds),
            c.root_delay,
            c.tlp_count,
        ),
    )
