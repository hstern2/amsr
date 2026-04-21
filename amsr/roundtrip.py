"""Shared round-trip logic: encode to AMSR, decode, generate conformer, compute RMSD."""

import os
import random
import time
from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdMolAlign

from .conf import GetConformer
from .decode import ToMol
from .encode import FromMol

N_RANDOM_SEEDS = 5  # number of randomized encodings per molecule


def Roundtrip(mol: Chem.Mol, seed=None) -> tuple[str, float, Chem.Mol]:
    """Encode mol to AMSR, decode, generate conformer.

    Args:
        mol: RDKit molecule with 3D conformer (reference geometry).
        seed: If not None, randomize the AMSR encoding with this seed.

    Returns (amsr_str, rmsd, mol_conformer).
    """
    if seed is not None:
        random.seed(seed)
    dihedral: dict[tuple[int, int, int, int], int] = {}
    s = FromMol(mol, randomize=(seed is not None))
    mol2 = ToMol(s, dihedral=dihedral)
    mol3 = GetConformer(mol2, dihedral=dihedral)
    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    return s, rmsd, mol3


def RoundtripSDF(
    sdf_path: str, seed, threshold: float = 1.0, output_dir: Optional[str] = None
) -> dict:
    """Round-trip one SDF file with one seed.

    Returns a dict with keys: name, seed, amsr, rmsd, status, time.
    If output_dir is set, writes orig and roundtrip SDF files there.
    """
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

    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)
        match = mol.GetSubstructMatch(mol_out)
        mol_orig = Chem.RenumberAtoms(mol, list(match)) if match else mol
        Chem.MolToMolFile(mol_orig, os.path.join(output_dir, f"{name}_{seed_str}_orig.sdf"))
        Chem.MolToMolFile(mol_out, os.path.join(output_dir, f"{name}_{seed_str}_out.sdf"))

    return {
        "name": name,
        "seed": seed_str,
        "amsr": s,
        "rmsd": rmsd,
        "status": status,
        "time": elapsed,
    }


def RoundtripSMI(smiles: str, name: str, seed) -> dict:
    """Round-trip one SMILES with one seed, using an InChI (-FixedH) comparison.

    Returns a dict with keys: name, smiles, seed, amsr, status, time.
    Status is PASSED when the post-roundtrip InChI matches the original,
    FAILED when it differs, or ERROR on exception.
    """
    seed_str = "0" if seed is None else str(seed + 1)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "name": name,
            "smiles": smiles,
            "seed": seed_str,
            "amsr": "",
            "status": "ERROR",
            "time": 0.0,
            "error": "MolFromSmiles returned None",
        }

    try:
        t0 = time.time()
        i1 = Chem.MolToInchi(mol, options="-FixedH")
        if seed is not None:
            random.seed(seed)
        s = FromMol(mol, randomize=(seed is not None))
        mol2 = ToMol(s)
        i2 = Chem.MolToInchi(mol2, options="-FixedH")
        elapsed = time.time() - t0
    except Exception as e:
        return {
            "name": name,
            "smiles": smiles,
            "seed": seed_str,
            "amsr": "",
            "status": "ERROR",
            "time": 0.0,
            "error": str(e),
        }

    status = "PASSED" if i1 == i2 else "FAILED"
    return {
        "name": name,
        "smiles": smiles,
        "seed": seed_str,
        "amsr": s,
        "status": status,
        "time": elapsed,
    }
