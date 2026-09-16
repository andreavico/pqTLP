#!/usr/bin/env -S sage -python
"""Benchmark one complete compiled puzzle with explicit parameters."""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pqtlp as tlp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-delay", type=int, required=True)
    parser.add_argument("--tlp-count", type=int, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    solution = b"compiled-bench".ljust(tlp.SOLUTION_BYTES, b"\0")
    start = perf_counter()
    puzzle = tlp.pgen(
        args.root_delay,
        solution,
        tlp_count=args.tlp_count,
        workers=args.workers,
    )
    generated = perf_counter()
    recovered = tlp.psolve(puzzle)
    solved = perf_counter()
    if recovered != solution:
        raise RuntimeError("compiled benchmark did not round-trip")
    encoded_size = len(pickle.dumps(puzzle, protocol=pickle.HIGHEST_PROTOCOL))
    print(
        f"root_delay={args.root_delay} tlp_count={args.tlp_count}"
        f" workers={args.workers} pgen_s={generated - start:.6f}"
        f" psolve_s={solved - generated:.6f} bytes={encoded_size}",
        flush=True,
    )


if __name__ == "__main__":
    main()
