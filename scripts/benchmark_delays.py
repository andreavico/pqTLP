#!/usr/bin/env -S sage -python
"""Benchmark full-size base puzzles at several CGL delay lengths."""

from __future__ import annotations

import argparse
import hashlib
import statistics
import sys
from pathlib import Path
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pqtlp as tlp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chunks",
        default="1,2,3,4",
        help="comma-separated numbers of default-size CGL chunks",
    )
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    chunks = [int(value) for value in args.chunks.split(",")]
    if any(value <= 0 for value in chunks):
        raise ValueError("chunk counts must be positive")
    if args.repeats <= 0:
        raise ValueError("repeats must be positive")

    chunk_size = tlp.CHUNK_SIZE
    print(
        "PARAMS"
        f" p={tlp.P}"
        f" p_bits={tlp.P.nbits()}"
        f" a={tlp.PSI_TORSION_EXPONENT}"
        f" q={tlp.Q}"
        f" N_bits={tlp.VERTICAL_DEGREE.nbits()}"
        f" chunk_size={chunk_size}",
        flush=True,
    )
    rows = []
    solution = b"benchmark".ljust(tlp.SOLUTION_BYTES, b"\0")
    for chunk_count in chunks:
        delay = chunk_count * chunk_size
        samples = {"setup": [], "pgen": [], "solve": []}
        for repeat in range(args.repeats):
            seed = hashlib.shake_256(
                f"pqtlp-benchmark:{chunk_count}:{repeat}".encode()
            ).digest(tlp.WALK_SEED_BYTES)
            start = perf_counter()
            state = tlp.setup_delay(delay, walk_seed=seed)
            setup_done = perf_counter()
            puzzle = tlp.pgen_base(state, solution)
            pgen_done = perf_counter()
            recovered = tlp.psolve_base(puzzle)
            solve_done = perf_counter()
            if recovered != solution:
                raise RuntimeError("benchmark puzzle did not round-trip")

            values = {
                "setup": setup_done - start,
                "pgen": pgen_done - setup_done,
                "solve": solve_done - pgen_done,
            }
            for name, value in values.items():
                samples[name].append(value)
            print(
                "RUN"
                f" chunks={chunk_count} repeat={repeat + 1}"
                f" setup={values['setup']:.6f}"
                f" pgen={values['pgen']:.6f}"
                f" solve={values['solve']:.6f}",
                flush=True,
            )

        row = {
            "chunks": chunk_count,
            "delay": delay,
            **{name: statistics.median(values) for name, values in samples.items()},
        }
        rows.append(row)

    print("\nMEDIANS")
    print("chunks\tdelay\tsetup_s\tpgen_s\tsolve_s")
    for row in rows:
        print(
            f"{row['chunks']}\t{row['delay']}\t{row['setup']:.6f}"
            f"\t{row['pgen']:.6f}\t{row['solve']:.6f}"
        )


if __name__ == "__main__":
    main()
