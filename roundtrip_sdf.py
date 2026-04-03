#!/usr/bin/env python
"""Round-trip verification for SDF files: encode to AMSR, decode, compute RMSD.

Each molecule is tested with a default encoding plus N_RANDOM_SEEDS randomized
encodings (same as the pytest suite).

Usage:
    python roundtrip_sdf.py ~/sdf              # sequential
    python roundtrip_sdf.py ~/sdf -j 10        # 10 parallel workers
    python roundtrip_sdf.py ~/sdf -j 10 -t 0.5 # stricter threshold
"""

import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import typer
from rdkit import Chem

from amsr.roundtrip import N_RANDOM_SEEDS, Roundtrip

app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})

_SEEDS = [None] + list(range(N_RANDOM_SEEDS))


def _process_one(sdf_path: str, seed, threshold: float):
    """Round-trip one SDF file with one seed."""
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    seed_str = "0" if seed is None else str(seed + 1)
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    if mol is None:
        return {
            "name": name,
            "seed": seed_str,
            "amsr": "",
            "rmsd": "",
            "status": "ERROR",
            "time": 0.0,
        }

    try:
        t0 = time.time()
        s, rmsd, mol_out = Roundtrip(mol, seed=seed)
        elapsed = time.time() - t0
    except Exception as e:
        return {
            "name": name,
            "seed": seed_str,
            "amsr": "",
            "rmsd": "",
            "status": "ERROR",
            "time": 0.0,
            "error": str(e),
        }

    status = "PASSED" if rmsd < threshold else "FAILED"
    return {
        "name": name,
        "seed": seed_str,
        "amsr": s,
        "rmsd": rmsd,
        "status": status,
        "time": elapsed,
    }


@app.command()
def main(
    input_dir: Path = typer.Argument(..., help="Directory containing SDF files"),
    output_dir: Path = typer.Option(None, "--output", "-o", help="Output directory for CSV"),
    threshold: float = typer.Option(1.0, "--threshold", "-t", help="RMSD threshold for PASS/FAIL"),
    jobs: int = typer.Option(1, "--jobs", "-j", help="Number of parallel workers"),
):
    """Round-trip verification: encode each SDF to AMSR, decode, compute RMSD."""
    if not input_dir.is_dir():
        typer.echo(f"Error: {input_dir} is not a directory", err=True)
        raise typer.Exit(1)

    if output_dir is None:
        output_dir = Path("out")
    output_dir.mkdir(parents=True, exist_ok=True)

    sdf_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".sdf"))
    n_mols = len(sdf_files)
    n_seeds = len(_SEEDS)
    total = n_mols * n_seeds
    typer.echo(
        f"Found {n_mols} SDF files in {input_dir}"
        f" ({n_seeds} seeds each, {total} tests, jobs={jobs})"
    )

    # Build work items: (sdf_path, seed)
    work = []
    for f in sdf_files:
        path = str(input_dir / f)
        for seed in _SEEDS:
            work.append((path, seed))

    results: list[dict] = []
    n_pass = n_fail = n_error = 0

    def _report(r):
        nonlocal n_pass, n_fail, n_error
        if r["status"] == "PASSED":
            n_pass += 1
            mark = "\033[32mPASSED\033[0m"
        elif r["status"] == "FAILED":
            n_fail += 1
            mark = "\033[31mFAILED\033[0m"
        else:
            n_error += 1
            mark = "\033[33mERROR\033[0m"
        rmsd_str = f"rmsd={r['rmsd']:.3f}" if isinstance(r["rmsd"], float) else r.get("error", "")
        time_str = f"{r['time']:.2f}s" if r["time"] else ""
        sys.stdout.write(f"{r['name']}[seed{r['seed']}] {mark} {rmsd_str} {time_str}\n")
        sys.stdout.flush()

    if jobs <= 1:
        for path, seed in work:
            r = _process_one(path, seed, threshold)
            results.append(r)
            _report(r)
    else:
        futures = {}
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            for path, seed in work:
                fut = executor.submit(_process_one, path, seed, threshold)
                futures[fut] = (path, seed)

            for fut in as_completed(futures):
                try:
                    r = fut.result()
                except Exception:
                    path, seed = futures[fut]
                    name = os.path.splitext(os.path.basename(path))[0]
                    seed_str = "0" if seed is None else str(seed + 1)
                    r = {
                        "name": name,
                        "seed": seed_str,
                        "amsr": "",
                        "rmsd": "",
                        "status": "ERROR",
                        "time": 0.0,
                    }
                results.append(r)
                _report(r)

    # Sort by RMSD descending (worst first)
    results.sort(key=lambda r: -r["rmsd"] if isinstance(r["rmsd"], float) else 0.0)

    csv_path = output_dir / "roundtrip_results.csv"
    with open(csv_path, "w", newline="") as csvf:
        w = csv.writer(csvf)
        w.writerow(["name", "seed", "amsr", "rmsd", "status", "time_s"])
        for r in results:
            rmsd_str = f"{r['rmsd']:.3f}" if isinstance(r["rmsd"], float) else ""
            w.writerow([r["name"], r["seed"], r["amsr"], rmsd_str, r["status"], f"{r['time']:.3f}"])

    typer.echo(f"\n{'='*60}")
    typer.echo(f"{n_pass} passed, {n_fail} failed, {n_error} errors / {total} total")
    typer.echo(f"Results: {csv_path}")


if __name__ == "__main__":
    app()
