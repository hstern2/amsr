"""Shared round-trip logic: encode to AMSR, decode, generate conformer, compute RMSD."""

import logging

from rdkit import Chem
from rdkit.Chem import rdMolAlign

log = logging.getLogger(__name__)


def Roundtrip(mol: Chem.Mol) -> tuple[str, float, Chem.Mol]:
    """Encode mol to AMSR, decode, generate conformer. Returns (amsr_str, rmsd, mol_zmatrix)."""
    from .decode import ToMol
    from .encode import FromMol
    from .zmatrix import GetConformer

    dihedral: dict[tuple[int, int, int, int], int] = {}
    s = FromMol(mol)
    mol2 = ToMol(s, dihedral=dihedral)
    mol3 = GetConformer(mol2, dihedral=dihedral)
    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    # Log RMSD before and after optimization (for debugging)
    if mol3.HasProp("_embed_coords"):
        import numpy as np

        embed_coords = np.frombuffer(
            mol3.GetProp("_embed_coords").encode("latin-1"), dtype=np.float64
        ).reshape(-1, 3)
        embed_mol = Chem.RWMol(mol2)
        embed_conf = Chem.Conformer(mol2.GetNumAtoms())
        embed_conf.Set3D(True)
        for i in range(mol2.GetNumAtoms()):
            embed_conf.SetAtomPosition(i, embed_coords[i].tolist())
        embed_mol.RemoveAllConformers()
        embed_mol.AddConformer(embed_conf, assignId=True)
        embed_rmsd = rdMolAlign.GetBestRMS(embed_mol.GetMol(), mol)
        log.info("RMSD embed=%.3f opt=%.3f", embed_rmsd, rmsd)
    return s, rmsd, mol3
