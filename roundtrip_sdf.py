#!/usr/bin/env python
"""Round-trip verification for SDF files: encode to AMSR, decode, compute RMSD.

Uses the same round-trip logic as test_zmatrix.py but runs on an arbitrary
directory of SDF files.
"""

import csv
import os
import traceback
from pathlib import Path

import typer
from rdkit import Chem

import amsr

app = typer.Typer(context_settings={"help_option_names": ["-h", "--help"]})


def _process_one(mol, name, out_dir: Path):
    """Round-trip one molecule. Same logic as test_zmatrix.py::test_roundtrip."""
    s, rmsd_raw, mol_raw, rmsd_refined, mol_refined = amsr.Roundtrip(mol)

    match = mol.GetSubstructMatch(mol_raw)
    if match:
        mol_reordered = Chem.RenumberAtoms(mol, list(match))
    else:
        mol_reordered = mol
    Chem.MolToMolFile(mol_reordered, str(out_dir / f"{name}_original.sdf"))
    Chem.MolToMolFile(mol_raw, str(out_dir / f"{name}_raw.sdf"))
    Chem.MolToMolFile(mol_refined, str(out_dir / f"{name}_refined.sdf"))

    return s, rmsd_raw, rmsd_refined


@app.command()
def main(
    input_dir: Path = typer.Argument(..., help="Directory containing SDF files"),
    output_dir: Path = typer.Option(
        None, "--output", "-o", help="Output directory (default: out/)"
    ),
    threshold: float = typer.Option(0.8, "--threshold", "-t", help="RMSD threshold for OK/FAIL"),
):
    """Round-trip verification: encode each SDF to AMSR, decode, compute RMSD."""
    if not input_dir.is_dir():
        typer.echo(f"Error: {input_dir} is not a directory", err=True)
        raise typer.Exit(1)

    if output_dir is None:
        output_dir = Path("out")
    output_dir.mkdir(parents=True, exist_ok=True)

    sdf_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".sdf"))
    typer.echo(f"Found {len(sdf_files)} SDF files in {input_dir}")

    rows: list[list[str]] = []
    n_ok = 0
    n_fail = 0
    n_error = 0

    for i, fname in enumerate(sdf_files):
        name = os.path.splitext(fname)[0]
        mol = Chem.MolFromMolFile(str(input_dir / fname), removeHs=True)
        if mol is None:
            typer.echo(f"[{i+1}/{len(sdf_files)}] {fname}: SKIP (could not parse)")
            rows.append([name, "", "", "", "parse_error"])
            n_error += 1
            continue

        try:
            s, rmsd_raw, rmsd_refined = _process_one(mol, name, output_dir)
            status = "OK" if rmsd_refined < threshold else "FAIL"
            if status == "OK":
                n_ok += 1
            else:
                n_fail += 1
            typer.echo(
                f"[{i+1}/{len(sdf_files)}] {fname}:"
                f" raw={rmsd_raw:.3f} refined={rmsd_refined:.3f} {status}"
            )
            rows.append([name, s, f"{rmsd_raw:.3f}", f"{rmsd_refined:.3f}", status])
        except Exception as e:
            n_error += 1
            typer.echo(f"[{i+1}/{len(sdf_files)}] {fname}: ERROR ({e})")
            traceback.print_exc()
            rows.append([name, "", "", "", f"error: {e}"])

    # Sort by rmsd_refined descending (worst first), matching test_zmatrix.py
    def sort_key(row):
        try:
            return -float(row[3])
        except (ValueError, IndexError):
            return 0.0

    rows.sort(key=sort_key)

    csv_path = output_dir / "roundtrip_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["name", "amsr", "rmsd_raw", "rmsd_refined", "status"])
        w.writerows(rows)

    typer.echo(f"\nDone. OK={n_ok} FAIL={n_fail} ERROR={n_error} / {len(sdf_files)} total")
    typer.echo(f"Results: {csv_path}")


if __name__ == "__main__":
    app()
