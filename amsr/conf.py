from rdkit import Chem
from typing import Tuple, Optional, Dict
import numpy as np


def GetConformerAndEnergy(
    mol: Chem.Mol, dihedral: Optional[Dict[Tuple[int, int, int, int], int]] = None
) -> Chem.Mol:
    """Return a conformer

    :param mol: RDKit Mol
    :param dihedral: optional dictionary of dihedral angle constraints
    :return: RDKit Mol
    """
    mol = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol)
    mp = Chem.AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = Chem.AllChem.MMFFGetMoleculeForceField(mol, mp)
    if dihedral is not None:
        for (i, j, k, l), v in dihedral.items():
            if not mol.GetBondBetweenAtoms(j, k).IsInRing():
                Chem.rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), i, j, k, l, v)
                ff.MMFFAddTorsionConstraint(i, j, k, l, False, v, v, 1e3)
    ff.Minimize(maxIts=100000)
    Chem.rdMolTransforms.CanonicalizeConformer(mol.GetConformer())
    ener = ff.CalcEnergy()  # kcal/mol
    return Chem.RemoveHs(mol), ener


def GetRoundedDihedral(
    mol: Chem.Mol, dihedral: Tuple[int, int, int, int], ndeg: int
) -> int:
    """Return dihedral angle rounded to nearest ndeg degrees

    :param mol: RDKit Mol
    :param dihedral: tuple of four atom indices
    :param ndeg: int
    :return: rounded dihedral angle
    """
    return (
        round(
            Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(0), *dihedral) / ndeg
        )
        * ndeg
    )
