#!/usr/bin/env python
"""Round-trip verification for SDF files: encode to AMSR, decode, compute RMSD.

Supports parallel processing with --jobs N for multi-core speedup.
"""

import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import typer
from rdkit import Chem

import amsr

app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})


def _process_one(sdf_path: str, out_dir: str, threshold: float):
    """Round-trip one SDF file.  Designed to run in a worker process."""
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    if mol is None:
        return [name, "", "", "parse_error"]

    try:
        s, rmsd, mol_out = amsr.Roundtrip(mol)
    except Exception as e:
        return [name, "", "", f"error: {e}"]

    # Save output SDF files
    match = mol.GetSubstructMatch(mol_out)
    if match:
        mol_reordered = Chem.RenumberAtoms(mol, list(match))
    else:
        mol_reordered = mol
    Chem.MolToMolFile(mol_reordered, os.path.join(out_dir, f"{name}_original.sdf"))
    Chem.MolToMolFile(mol_out, os.path.join(out_dir, f"{name}_out.sdf"))

    status = "OK" if rmsd < threshold else "FAIL"
    return [name, s, f"{rmsd:.3f}", status]


@app.command()
def main(
    input_dir: Path = typer.Argument(..., help="Directory containing SDF files"),
    output_dir: Path = typer.Option(
        None, "--output", "-o", help="Output directory (default: out/)"
    ),
    threshold: float = typer.Option(1.0, "--threshold", "-t", help="RMSD threshold for OK/FAIL"),
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
    total = len(sdf_files)
    typer.echo(f"Found {total} SDF files in {input_dir} (jobs={jobs})")

    sdf_paths = [str(input_dir / f) for f in sdf_files]
    out_str = str(output_dir)
    rows: list[list[str]] = []
    n_ok = n_fail = n_error = 0

    if jobs <= 1:
        # Sequential processing
        for i, path in enumerate(sdf_paths):
            row = _process_one(path, out_str, threshold)
            rows.append(row)
            status = row[3]
            if status == "OK":
                n_ok += 1
            elif status == "FAIL":
                n_fail += 1
            else:
                n_error += 1
            typer.echo(
                f"[{i+1}/{total}] {os.path.basename(path)}: "
                f"{'rmsd=' + row[2] + ' ' if row[2] else ''}{status}"
            )
    else:
        # Parallel processing
        futures = {}
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            for i, path in enumerate(sdf_paths):
                fut = executor.submit(_process_one, path, out_str, threshold)
                futures[fut] = (i, path)

            done = 0
            for fut in as_completed(futures):
                done += 1
                i, path = futures[fut]
                try:
                    row = fut.result()
                except Exception as e:
                    name = os.path.splitext(os.path.basename(path))[0]
                    row = [name, "", "", f"error: {e}"]
                rows.append(row)
                status = row[3]
                if status == "OK":
                    n_ok += 1
                elif status == "FAIL":
                    n_fail += 1
                else:
                    n_error += 1
                if done % 100 == 0 or done == total:
                    typer.echo(f"  [{done}/{total}] OK={n_ok} FAIL={n_fail} ERROR={n_error}")

    # Sort by rmsd descending (worst first)
    def sort_key(row):
        try:
            return -float(row[2])
        except (ValueError, IndexError):
            return 0.0

    rows.sort(key=sort_key)

    csv_path = output_dir / "roundtrip_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "amsr", "rmsd", "status"])
        w.writerows(rows)

    typer.echo(f"\nDone. OK={n_ok} FAIL={n_fail} ERROR={n_error} / {total} total")
    typer.echo(f"Results: {csv_path}")


if __name__ == "__main__":
    app()
