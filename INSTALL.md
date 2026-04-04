# Install and test

## Install

```bash
pip install -e ".[dev]"
```

## Test

```bash
pytest -n 8             # run all tests in parallel (8 workers)
pytest                  # run all tests sequentially
```

The conformer tests (`test_conf.py`) encode each SDF molecule with 6
different AMSR encodings (1 canonical + 5 randomized) and check that
the round-trip RMSD is below 0.8 Å.

## Optional: C-accelerated conformer generation

```bash
make            # build C extension (~3.5x faster)
make clean      # revert to pure Python
```

Requires a C compiler.  Used automatically when present.

## Batch processing

```bash
python roundtrip_sdf.py /path/to/sdf/dir -j 8     # 8 parallel workers
```
