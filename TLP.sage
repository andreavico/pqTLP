#!/usr/bin/env python3
"""Generate or solve a pqTLP puzzle from the command line."""

from __future__ import annotations

import argparse
import base64
import gzip
import os
import pickle
import sys
from time import perf_counter


def ensure_sage_python() -> None:
    """Restart under ``sage -python`` when invoked with ordinary Python."""

    try:
        import sage.all  # noqa: F401

        return
    except ImportError:
        pass

    if os.environ.get("PQTLP_IN_SAGE") == "1":
        raise RuntimeError("SageMath is required")
    environment = dict(os.environ, PQTLP_IN_SAGE="1")
    os.execvpe("sage", ["sage", "-python", __file__, *sys.argv[1:]], environment)


ensure_sage_python()

import pqtlp as tlp  # noqa: E402

MAGIC = "PQTLP"


def encode(envelope: dict) -> str:
    payload = pickle.dumps(envelope, protocol=pickle.HIGHEST_PROTOCOL)
    return f"{MAGIC}\n{base64.b64encode(gzip.compress(payload)).decode()}\n"


def decode(text: str) -> dict:
    fields = text.strip().split(None, 1)
    if len(fields) != 2 or fields[0] != MAGIC:
        raise ValueError("not a pqTLP puzzle")
    payload = gzip.decompress(base64.b64decode(fields[1]))
    envelope = pickle.loads(payload)
    if not isinstance(envelope, dict) or envelope.get("kind") != "pqtlp":
        raise ValueError("invalid pqTLP puzzle")
    return envelope


def read(path: str) -> str:
    if path == "-":
        return sys.stdin.read()
    with open(path, encoding="ascii") as puzzle_file:
        return puzzle_file.read()


def parse_solution(args) -> tuple[bytes, str, int]:
    if args.secret_hex is not None:
        value = bytes.fromhex(args.secret_hex)
        output_format = "hex"
    elif args.secret_text is not None:
        value = args.secret_text.encode()
        output_format = "text"
    else:
        return os.urandom(tlp.SOLUTION_BYTES), "hex", tlp.SOLUTION_BYTES

    if not value or len(value) > tlp.SOLUTION_BYTES:
        raise ValueError(f"the secret must contain 1--{tlp.SOLUTION_BYTES} bytes")
    return value.ljust(tlp.SOLUTION_BYTES, b"\0"), output_format, len(value)


def generate(args) -> None:
    choice = tlp.choose_parameters(args.generate)
    solution, output_format, output_length = parse_solution(args)
    print(
        f"T0={choice.root_delay}, T1={choice.tlp_count}, "
        f"predicted solve={choice.predicted_solve_seconds:.1f}s",
        file=sys.stderr,
    )

    started = perf_counter()
    puzzle = tlp.pgen(
        choice.root_delay,
        solution,
        tlp_count=choice.tlp_count,
        workers=args.workers,
    )
    elapsed = perf_counter() - started
    envelope = {
        "kind": "pqtlp",
        "target_seconds": float(args.generate),
        "predicted_solve_seconds": choice.predicted_solve_seconds,
        "generation_seconds": elapsed,
        "solution_format": output_format,
        "solution_length": output_length,
        "puzzle": puzzle,
    }
    sys.stdout.write(encode(envelope))
    print(f"generation={elapsed:.3f}s", file=sys.stderr)

    if args.reveal_secret:
        revealed = solution[:output_length]
        value = revealed.decode() if output_format == "text" else revealed.hex()
        print(f"secret={value}", file=sys.stderr)


def solve(args) -> None:
    envelope = decode(read(args.solve))
    started = perf_counter()
    solution = tlp.psolve(envelope["puzzle"])
    elapsed = perf_counter() - started

    length = int(envelope["solution_length"])
    solution = solution[:length]
    if envelope["solution_format"] == "text":
        print(solution.decode())
    else:
        print(solution.hex())
    print(f"solve={elapsed:.3f}s", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    action = result.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--generate",
        type=float,
        metavar="SECONDS",
        help="generate a puzzle for an estimated solve duration",
    )
    action.add_argument(
        "--solve",
        metavar="PUZZLE",
        help="solve a trusted puzzle file, or - for standard input",
    )
    result.add_argument(
        "--workers", type=int, default=1, help="generation processes (default: 1)"
    )
    secret = result.add_mutually_exclusive_group()
    secret.add_argument("--secret-hex", help="solution bytes encoded as hex")
    secret.add_argument("--secret-text", help="solution encoded as UTF-8")
    result.add_argument(
        "--reveal-secret",
        action="store_true",
        help="print the generated secret to standard error",
    )
    return result


def main() -> int:
    argument_parser = parser()
    args = argument_parser.parse_args()
    if args.workers <= 0:
        argument_parser.error("--workers must be positive")
    if args.solve and (
        args.secret_hex is not None
        or args.secret_text is not None
        or args.reveal_secret
    ):
        argument_parser.error("secret options only apply to --generate")

    try:
        generate(args) if args.generate is not None else solve(args)
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
