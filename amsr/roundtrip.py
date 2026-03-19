"""Shared round-trip logic: encode to AMSR, decode, generate conformer, compute RMSD."""

from rdkit import Chem
from rdkit.Chem import rdMolAlign


def Roundtrip(mol: Chem.Mol) -> tuple[str, float, Chem.Mol, float, Chem.Mol]:
    """Encode mol to AMSR, decode, generate conformer.

    Returns (amsr_str, rmsd_raw, mol_raw, rmsd_refined, mol_refined).
    rmsd_raw/mol_raw: z-matrix placement only (no ring refinement).
    rmsd_refined/mol_refined: with ring closure refinement.
    """
    from .decode import ToMol
    from .encode import FromMol
    from .zmatrix import GetConformer

    dihedral: dict[tuple[int, int, int, int], int] = {}
    s = FromMol(mol)
    mol2 = ToMol(s, dihedral=dihedral)
    mol_raw = GetConformer(mol2, dihedral=dihedral, refine_rings=False)
    rmsd_raw = rdMolAlign.GetBestRMS(mol_raw, mol)
    mol_refined = GetConformer(mol2, dihedral=dihedral, refine_rings=True)
    rmsd_refined = rdMolAlign.GetBestRMS(mol_refined, mol)
    return s, rmsd_raw, mol_raw, rmsd_refined, mol_refined
