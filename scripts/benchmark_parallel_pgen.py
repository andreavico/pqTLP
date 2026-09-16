#!/usr/bin/env -S sage -python
"""Benchmark independent base-puzzle generation across worker counts."""

from __future__ import annotations

import argparse
import pickle
import statistics
import sys
from pathlib import Path
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pqtlp as tlp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=32)
    parser.add_argument("--workers", default="1,4")
    parser.add_argument("--repeats", type=int, default=1)
    args = parser.parse_args()
    if args.count <= 0:
        raise ValueError("count must be positive")
    if args.repeats <= 0:
        raise ValueError("repeats must be positive")
    worker_counts = [int(value) for value in args.workers.split(",")]
    if any(workers <= 0 for workers in worker_counts):
        parser.error("worker counts must be positive")

    state = tlp.setup_delay(
        tlp.CHUNK_SIZE,
        walk_seed=b"parallel-pgen-benchmark".ljust(tlp.WALK_SEED_BYTES, b"\0"),
    )
    solutions = [
        index.to_bytes(tlp.SOLUTION_BYTES, "big") for index in range(args.count)
    ]
    samples = {workers: {"batch": [], "size": []} for workers in worker_counts}
    for repeat in range(args.repeats):
        # Reverse alternate repetitions to reduce bias from temperature and
        # background-load drift between the serial and parallel measurements.
        order = worker_counts if repeat % 2 == 0 else list(reversed(worker_counts))
        for workers in order:
            start = perf_counter()
            puzzles = tlp.generate_base_puzzles(
                state,
                solutions,
                workers=workers,
            )
            elapsed = perf_counter() - start
            if len(puzzles) != args.count:
                raise RuntimeError("parallel PGen returned the wrong puzzle count")
            if (
                len({puzzle.psi.codomain.j_invariant() for puzzle in puzzles})
                != args.count
            ):
                raise RuntimeError(
                    "parallel workers did not generate independent puzzles"
                )
            encoded_sizes = [
                len(pickle.dumps(puzzle, protocol=pickle.HIGHEST_PROTOCOL))
                for puzzle in puzzles
            ]
            samples[workers]["batch"].append(elapsed)
            samples[workers]["size"].extend(encoded_sizes)
            print(
                f"RUN repeat={repeat + 1} workers={workers}"
                f" batch_s={elapsed:.6f}"
                f" seconds_per_puzzle={elapsed / args.count:.6f}"
                f" median_bytes={statistics.median(encoded_sizes):.1f}",
                flush=True,
            )

    serial_per_puzzle = None
    if 1 in samples:
        serial_per_puzzle = statistics.median(samples[1]["batch"]) / args.count
    print("\nMEDIANS")
    print("workers\tbatch_s\tseconds_per_puzzle\tspeedup\tmedian_bytes")
    for workers in worker_counts:
        batch = statistics.median(samples[workers]["batch"])
        per_puzzle = batch / args.count
        speedup = serial_per_puzzle / per_puzzle if serial_per_puzzle else float("nan")
        size = statistics.median(samples[workers]["size"])
        print(f"{workers}\t{batch:.6f}\t{per_puzzle:.6f}\t{speedup:.3f}\t{size:.1f}")


if __name__ == "__main__":
    main()
