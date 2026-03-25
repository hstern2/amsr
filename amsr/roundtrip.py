"""Shared round-trip logic: encode to AMSR, decode, generate conformer, compute RMSD."""

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def Roundtrip(mol: Chem.Mol) -> tuple[str, float, Chem.Mol]:
    """Encode mol to AMSR, decode, generate conformer.

    Returns (amsr_str, rmsd, mol_conformer).
    """
    from .decode import ToMol
    from .encode import FromMol
    from .zmatrix import GetConformer

    dihedral: dict[tuple[int, int, int, int], int] = {}
    s = FromMol(mol)
    mol2 = ToMol(s, dihedral=dihedral)
    mol3 = GetConformer(mol2, dihedral=dihedral)
    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    return s, rmsd, mol3
