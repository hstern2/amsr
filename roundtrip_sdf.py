#!/usr/bin/env python
"""Round-trip verification for SDF files: encode to AMSR, decode, compute RMSD."""

import csv
import os
import traceback
from pathlib import Path

import typer
from rdkit import Chem
from rdkit.Chem import rdMolAlign

import amsr

app = typer.Typer()


def _roundtrip(mol, name, out_dir: Path):
    """Encode 3D mol to AMSR, decode, generate z-matrix conformer.
    Align to original, write both SDFs, return RMSD."""
    dihedral: dict[tuple[int, int, int, int], int] = {}
    s = amsr.FromMol(mol)
    mol2 = amsr.ToMol(s, dihedral=dihedral)
    mol3 = amsr.GetConformer(mol2, dihedral=dihedral)

    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    match = mol.GetSubstructMatch(mol3)
    if match:
        atom_map = [(i, match[i]) for i in range(mol3.GetNumAtoms())]
        rdMolAlign.AlignMol(mol3, mol, atomMap=atom_map)

    Chem.MolToMolFile(mol, str(out_dir / f"{name}_original.sdf"))
    Chem.MolToMolFile(mol3, str(out_dir / f"{name}_zmatrix.sdf"))

    return s, rmsd


@app.command()
def main(
    input_dir: Path = typer.Argument(..., help="Directory containing SDF files"),
    output_dir: Path = typer.Option(
        None, "--output", "-o", help="Output directory (default: out/)"
    ),
    threshold: float = typer.Option(1.0, "--threshold", "-t", help="RMSD threshold for OK/FAIL"),
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

    results = []
    n_ok = 0
    n_fail = 0
    n_error = 0

    for i, fname in enumerate(sdf_files):
        name = os.path.splitext(fname)[0]
        mol = Chem.MolFromMolFile(str(input_dir / fname), removeHs=True)
        if mol is None:
            typer.echo(f"[{i+1}/{len(sdf_files)}] {fname}: SKIP (could not parse)")
            results.append(
                {"name": name, "file": fname, "amsr": "", "rmsd": "", "status": "parse_error"}
            )
            n_error += 1
            continue

        try:
            s, rmsd = _roundtrip(mol, name, output_dir)
            status = "OK" if rmsd < threshold else "FAIL"
            if status == "OK":
                n_ok += 1
            else:
                n_fail += 1
            typer.echo(f"[{i+1}/{len(sdf_files)}] {fname}: RMSD={rmsd:.3f} {status}")
            results.append(
                {"name": name, "file": fname, "amsr": s, "rmsd": f"{rmsd:.3f}", "status": status}
            )
        except Exception as e:
            n_error += 1
            typer.echo(f"[{i+1}/{len(sdf_files)}] {fname}: ERROR ({e})")
            traceback.print_exc()
            results.append(
                {"name": name, "file": fname, "amsr": "", "rmsd": "", "status": f"error: {e}"}
            )

    csv_path = output_dir / "roundtrip_results.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["name", "file", "amsr", "rmsd", "status"])
        writer.writeheader()
        writer.writerows(results)

    typer.echo(f"\nDone. OK={n_ok} FAIL={n_fail} ERROR={n_error} / {len(sdf_files)} total")
    typer.echo(f"Results: {csv_path}")
    typer.echo(f"SDF output: {output_dir}")


if __name__ == "__main__":
    app()
