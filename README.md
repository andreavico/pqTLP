# pqTLP

SageMath implementation accompanying the paper. It implements the compact
seeded CGL delay, the ideal-side shortcut, the public two-dimensional
pushforward, and the flat plain-model compiler.

This is research software. It is not constant-time or side-channel resistant,
and its pickle-based encoding must only be loaded from trusted sources.

## Requirements

- SageMath 10.8
- A Unix-like system for multi-process puzzle generation

No separate Python packages are required.

## Parameters

The implementation has one fixed parameter set:

| parameter | value |
|---|---|
| security level | $\lambda=128$ |
| field characteristic | $p=3\cdot2^{324}-1$ |
| vertical prime | $q=\operatorname{nextprime}(2^{256})$ |
| vertical degree | $N=q(2^{257}-q)$ |
| solution length | 128 bits |
| walk seed | 256 bits |
| CGL chunk size | 67 degree-two steps |

The public walk is represented by its seed and length. Setup checks each
chunk as it advances and retains only the base-curve witness, transport data,
and final shortcut.

## Quick check

Run from the extracted repository directory:

```sh
sage -python -m unittest discover -s tests -v
```

The suite checks the concrete parameters, full-parameter base and compiled
puzzle round trips, parallel generation, command-line serialization,
calibration, and arithmetic invariants. It normally takes under a minute on
the reference machine; this is not a security audit.

## Python API

```python
import pqtlp

solution = b"example".ljust(pqtlp.SOLUTION_BYTES, b"\0")
puzzle = pqtlp.pgen(
    root_delay=67,
    message=solution,
    tlp_count=2,
    workers=1,
)
assert pqtlp.psolve(puzzle) == solution
```

Run these examples under `sage -python` or in a SageMath session.
`root_delay` must be a positive multiple of 67. `tlp_count` defaults to
`root_delay`, the paper's balanced setting. The optional `workers` argument only
parallelizes the independent base-puzzle generation calls; setup and solving
remain sequential.

For direct access to one base puzzle:

```python
state = pqtlp.setup_delay(67)
puzzle = pqtlp.pgen_base(state, solution)
assert pqtlp.psolve_base(puzzle) == solution
```

## Command line

`TLP.sage` automatically restarts itself under `sage -python`:

```sh
python3 TLP.sage --generate 60 --workers 4 --secret-text hello > puzzle
python3 TLP.sage --solve puzzle
```

Secrets contain 1–16 bytes and are padded internally; the original length is
restored on solving. Omit both secret options to generate a random solution,
or use `--secret-hex` to supply bytes. `--solve -` reads from standard input.

The requested time is converted to a calibrated pair $(T_0,T_1)$. Generation
progress is written to standard error and the encoded puzzle to standard
output. Calibration uses the recorded serial generation timings, even when
`--workers` is supplied. The requested duration is an estimate on the reference
machine, not a guaranteed minimum solve time.

## Reproducing the results

```sh
sage -python scripts/benchmark_delays.py --chunks 1,2,4,8,16,32 --repeats 3
sage -python scripts/benchmark_delays.py --chunks 64 --repeats 2
sage -python scripts/benchmark_parallel_pgen.py --count 128 --workers 1,4 --repeats 3
sage -python scripts/project_compiled_delays.py --workers 1,4
sage -python scripts/benchmark_compiled.py --root-delay 201 --tlp-count 9
sage -python scripts/benchmark_compiled.py --root-delay 134 --tlp-count 18 --workers 4
```

The projection script interpolates within the measured range and extrapolates
beyond it; the command-line chooser
uses measured delay lengths only. Timings vary with random ideal sampling,
SageMath version, and hardware. Benchmark seeds fix the public walk seed, but
the ideal sampling remains randomized.

## Layout

- `pqtlp/protocol.py`: setup, generation, solving, and flat compilation
- `pqtlp/delay.py`: setup-side seeded CGL walk and transport ideals
- `pqtlp/curves.py`: public two-dimensional pushforward
- `pqtlp/ideals.py`: shortcut and vertical-ideal operations
- `pqtlp/parameters.py`: the concrete parameter set
- `pqtlp/backend/`: bundled Qlapoti and theta arithmetic
- `scripts/`: benchmark and projection scripts
- `tests/`: protocol, arithmetic, calibration, and command-line tests

## Licence

The original pqTLP code is distributed under the [MIT licence](LICENSE),
copyright (c) 2026 pqTLP contributors. Bundled third-party arithmetic retains
its original licences; see [NOTICE](NOTICE) and [LICENSES/](LICENSES/).
