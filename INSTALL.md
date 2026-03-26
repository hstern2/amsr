# Installation

## Python package

```bash
pip install -e .
```

### Dependencies

Defined in `pyproject.toml`:

- rdkit, networkx, anytree, Levenshtein, typer, pandas, torch
- For development: `pip install -e ".[dev]"` (adds pytest, ruff, mypy, pre-commit)

Also requires `scipy` at runtime for conformer generation (`amsr.zmatrix`).

## Optional: C extension

The C extension accelerates the ring geometry optimizer (~10x faster
cost/gradient evaluation).  It is **not required** -- the pure-Python
implementation is used by default.

### Build

```bash
make            # builds amsr/cost_grad.dylib (macOS) or .so (Linux)
make clean      # removes built libraries
```

Requires a C compiler (`cc` / `gcc` / `clang`).  No external libraries
beyond `-lm`.  Tested on:

- macOS ARM64 (Apple M2 Pro, Xcode clang)
- Linux x86_64 (gcc)

### Enable

Set the environment variable before running:

```bash
export AMSR_USE_C=1
```

The C extension produces results that are numerically very close to
Python (differences < 1e-9 per evaluation) but floating-point ordering
differences can cause the L-BFGS-B optimizer to find different local
minima for borderline molecules.  For this reason, **Python is the
default** to ensure reproducible results.

## Batch processing

`roundtrip_sdf.py` supports parallel execution:

```bash
python roundtrip_sdf.py /path/to/sdf/dir -j 8    # 8 parallel workers
python roundtrip_sdf.py /path/to/sdf/dir -j 8 -o results/
```

## Running tests

```bash
pytest tests/
```
