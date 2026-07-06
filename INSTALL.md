# Install and test

## Prerequisites

- Python 3.12+
- [uv](https://docs.astral.sh/uv/)
- C compiler (cc/gcc/clang)

## Install

```bash
uv sync --extra dev
```

## Build C extension

```bash
make
```

Required for conformer generation. Builds `amsr/src/conf_util.dylib` (macOS)
or `amsr/src/conf_util.so` (Linux) from C sources in `amsr/src/`.

## Test

```bash
uv run --extra dev pytest -n 8     # run all tests in parallel (8 workers)
uv run --extra dev pytest          # run all tests sequentially
```

The conformer tests (`test_conf.py`) encode each SDF molecule with 6
different AMSR encodings (1 canonical + 5 randomized) and check that
the round-trip RMSD is below 1.2 Å.

## Batch processing

```bash
uv run python roundtrip_sdf.py /path/to/sdf/dir -j 8     # 8 parallel workers
```
