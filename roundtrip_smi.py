#!/usr/bin/env python
"""Round-trip verification for .smi files: encode to AMSR, decode, compare InChI.

Each molecule is tested with a default encoding plus N_RANDOM_SEEDS randomized
encodings, and the result is recorded as PASSED/FAILED/ERROR based on an
InChI (-FixedH) comparison of the original vs. round-tripped molecule.

Usage:
    python roundtrip_smi.py ~/smi              # sequential
    python roundtrip_smi.py ~/smi -j 10        # 10 parallel workers
    python roundtrip_smi.py ~/smi -n 5         # 5 random seeds per molecule
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated

import typer

from amsr.roundtrip import RoundtripSMI

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _iter_smi_records(path: str):
    """Yield (smiles, name) pairs from a .smi file, skipping blanks and comments."""
    stem = os.path.splitext(os.path.basename(path))[0]
    with open(path) as fh:
        for line_no, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else f"{stem}_{line_no}"
            yield smi, name


def _smi_work(input_dir: Path, n_seeds: int):
    """Yield (smiles, name, seed) by recursively scanning directories for .smi files."""
    seeds = [None] + list(range(n_seeds))
    for dirpath, dirnames, filenames in os.walk(input_dir):
        dirnames.sort()
        for fname in sorted(filenames):
            if fname.endswith(".smi"):
                full = os.path.join(dirpath, fname)
                for smi, name in _iter_smi_records(full):
                    for seed in seeds:
                        yield smi, name, seed


def _report(r, counts):
    counts[3] += 1
    if r["status"] == "PASSED":
        counts[0] += 1
        mark = "\033[32mPASSED\033[0m"
    elif r["status"] == "FAILED":
        counts[1] += 1
        mark = "\033[31mFAILED\033[0m"
    else:
        counts[2] += 1
        mark = "\033[33mERROR\033[0m"
    info = r.get("error", "")
    time_str = f"{r['time']:.2f}s" if r["time"] else ""
    n = counts[3]
    sys.stdout.write(f"[{n}] {r['name']}[seed{r['seed']}] {mark} {info} {time_str}\n")
    if n % 100 == 0:
        sys.stdout.write(
            f"  --- {counts[0]} passed, {counts[1]} failed, {counts[2]} errors / {n} total ---\n"
        )
    sys.stdout.flush()


def _error_result(smi, name, seed):
    return {
        "name": name,
        "smiles": smi,
        "seed": "0" if seed is None else str(seed + 1),
        "amsr": "",
        "status": "ERROR",
        "time": 0.0,
    }


@app.command()
def main(
    input_dir: Annotated[
        Path, typer.Argument(help="Directory (recursively searched) containing .smi files")
    ],
    jobs: Annotated[int, typer.Option("-j", help="Number of parallel workers")] = (
        os.cpu_count() or 1
    ),
    nseeds: Annotated[int, typer.Option("-n", help="Number of random seeds per molecule")] = 10,
):
    """Round-trip verification: encode each SMILES to AMSR, decode, compare InChI."""
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    csv_path = os.path.join(".", "roundtrip_smi_out.csv")
    csv_file = open(csv_path, "w", newline="")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["name", "smiles", "seed", "amsr", "status", "time_s"])

    counts = [0, 0, 0, 0]  # pass, fail, error, total

    def _record(r):
        _report(r, counts)
        csv_writer.writerow(
            [r["name"], r["smiles"], r["seed"], r["amsr"], r["status"], f"{r['time']:.3f}"]
        )
        csv_file.flush()

    if jobs <= 1:
        for smi, name, seed in _smi_work(input_dir, nseeds):
            _record(RoundtripSMI(smi, name, seed))
    else:
        max_pending = jobs * 2
        work = _smi_work(input_dir, nseeds)
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {}
            exhausted = False
            for smi, name, seed in work:
                futures[executor.submit(RoundtripSMI, smi, name, seed)] = (smi, name, seed)
                if len(futures) >= max_pending:
                    break
            else:
                exhausted = True

            while futures:
                done = next(iter(as_completed(futures)))
                smi, name, seed = futures.pop(done)
                try:
                    _record(done.result())
                except Exception:
                    _record(_error_result(smi, name, seed))
                if not exhausted:
                    for smi, name, seed in work:
                        futures[executor.submit(RoundtripSMI, smi, name, seed)] = (
                            smi,
                            name,
                            seed,
                        )
                        if len(futures) >= max_pending:
                            break
                    else:
                        exhausted = True

    csv_file.close()
    total = counts[3]
    print(f"\n{'='*60}")
    print(f"{counts[0]} passed, {counts[1]} failed, {counts[2]} errors / {total} total")
    print(f"Results: {csv_path}")


if __name__ == "__main__":
    app()
