#!/usr/bin/env -S sage -python
"""Choose generation-time-minimizing parameters from measured base timings."""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pqtlp.calibration import DELAY_PROFILE
from pqtlp.parameters import CHUNK_SIZE

TARGETS = (
    ("1 minute", 60),
    ("10 minutes", 10 * 60),
    ("30 minutes", 30 * 60),
    ("1 hour", 60 * 60),
    ("2 hours", 2 * 60 * 60),
    ("12 hours", 12 * 60 * 60),
    ("1 day", 24 * 60 * 60),
    ("1 week", 7 * 24 * 60 * 60),
    ("2 weeks", 14 * 24 * 60 * 60),
    ("1 month", 30 * 24 * 60 * 60),
    ("6 months", 180 * 24 * 60 * 60),
    ("1 year", 365 * 24 * 60 * 60),
)
BASE_PUZZLE_BYTES = 4085
# Medians of three interleaved 128-puzzle batches.
SERIAL_PGEN_SECONDS = 1.231958
PARALLEL_PGEN_SECONDS = {4: 0.374458}


def fit_line(field: str) -> tuple[float, float]:
    points = [
        (root_delay / CHUNK_SIZE, float(row[field]))
        for root_delay, row in sorted(DELAY_PROFILE.items())
    ]
    x_mean = statistics.mean(x for x, _ in points)
    y_mean = statistics.mean(y for _, y in points)
    denominator = sum((x - x_mean) ** 2 for x, _ in points)
    slope = sum((x - x_mean) * (y - y_mean) for x, y in points) / denominator
    return y_mean - slope * x_mean, slope


def estimate_in_range_or_extrapolate(
    chunks: int,
    field: str,
    intercept: float,
    slope: float,
) -> float:
    """Interpolate measured medians, using the fitted line only past them."""

    points = [
        (root_delay // CHUNK_SIZE, float(row[field]))
        for root_delay, row in sorted(DELAY_PROFILE.items())
    ]
    for measured_chunks, value in points:
        if chunks == measured_chunks:
            return value
    for (left_chunks, left_value), (right_chunks, right_value) in zip(
        points,
        points[1:],
    ):
        if left_chunks < chunks < right_chunks:
            fraction = (chunks - left_chunks) / (right_chunks - left_chunks)
            return left_value + fraction * (right_value - left_value)
    return intercept + slope * chunks


def pretty_seconds(seconds: float) -> str:
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 2 * 60 * 60:
        return f"{seconds / 60:.1f} min"
    if seconds < 2 * 24 * 60 * 60:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} d"


def pretty_bytes(size: int) -> str:
    if size < 1024 * 1024:
        return f"{size / 1024:.1f} KiB"
    return f"{size / (1024 * 1024):.2f} MiB"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        default="1,4",
        help="comma-separated worker counts (measured values exist for 1 and 4)",
    )
    args = parser.parse_args()
    worker_counts = [int(value) for value in args.workers.split(",")]

    if any(workers not in {1, *PARALLEL_PGEN_SECONDS} for workers in worker_counts):
        parser.error("measured throughput is available only for 1 or 4 workers")

    setup_intercept, setup_slope = fit_line("setup")
    solve_intercept, solve_slope = fit_line("base_solve")
    serial_per_puzzle = SERIAL_PGEN_SECONDS
    max_measured_chunks = max(DELAY_PROFILE) // CHUNK_SIZE

    print(
        "model: "
        f"setup(k)={setup_intercept:.4f}+{setup_slope:.4f}k s, "
        f"solve(k)={solve_intercept:.4f}+{solve_slope:.4f}k s, "
        f"serial PGen/base={serial_per_puzzle:.4f} s"
    )

    for workers in worker_counts:
        if workers == 1:
            per_puzzle = serial_per_puzzle
        elif workers in PARALLEL_PGEN_SECONDS:
            per_puzzle = PARALLEL_PGEN_SECONDS[workers]
        else:
            raise ValueError(f"no measured parallel throughput for {workers} workers")
        print(f"\nworkers={workers}, effective PGen/base={per_puzzle:.5f} s")
        print(
            "| target | T0 (2-isogenies) | T1 | setup | online PGen | "
            "total generation | projected solve | approximate size | basis |"
        )
        print("|---:|---:|---:|---:|---:|---:|---:|---:|:---|")

        for label, target in TARGETS:
            continuous = math.sqrt(per_puzzle * target / (setup_slope * solve_slope))
            limit = max(32, math.ceil(4 * continuous + 32))
            best = None
            for chunks in range(1, limit + 1):
                base_solve = max(
                    float(DELAY_PROFILE[CHUNK_SIZE]["base_solve"]),
                    estimate_in_range_or_extrapolate(
                        chunks,
                        "base_solve",
                        solve_intercept,
                        solve_slope,
                    ),
                )
                setup = max(
                    float(DELAY_PROFILE[CHUNK_SIZE]["setup"]),
                    estimate_in_range_or_extrapolate(
                        chunks,
                        "setup",
                        setup_intercept,
                        setup_slope,
                    ),
                )
                count = math.ceil(target / base_solve)
                online_pgen = count * per_puzzle
                total_generation = setup + online_pgen
                candidate = (
                    total_generation,
                    chunks,
                    count,
                    setup,
                    online_pgen,
                    base_solve,
                )
                if best is None or candidate < best:
                    best = candidate

            total_generation, chunks, count, setup, online_pgen, base_solve = best
            projected_solve = count * base_solve
            basis = (
                "within measured range"
                if chunks <= max_measured_chunks
                else "extrapolated"
            )
            print(
                f"| {label} | {CHUNK_SIZE * chunks:,} | {count:,} | "
                f"{pretty_seconds(setup)} | {pretty_seconds(online_pgen)} | "
                f"{pretty_seconds(total_generation)} | "
                f"{pretty_seconds(projected_solve)} | "
                f"{pretty_bytes(count * BASE_PUZZLE_BYTES)} | {basis} |"
            )


if __name__ == "__main__":
    main()
