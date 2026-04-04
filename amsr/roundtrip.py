"""Shared round-trip logic: encode to AMSR, decode, generate conformer, compute RMSD."""

import random

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
