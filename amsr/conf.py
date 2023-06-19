from rdkit import Chem
import rdkit.Chem.AllChem
from typing import Tuple, Optional, Dict


def GetConformer(
    mol: Chem.Mol, dihedral: Optional[Dict[Tuple[int, int, int, int], int]] = None
) -> Chem.Mol:
    """Return a conformer

    :param mol: RDKit Mol
    :param dihedral: optional dictionary of dihedral angle constraints
    :return: RDKit Mol
    """
    mol = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol)
    if dihedral is not None:
        for (i, j, k, l), v in dihedral.items():
            if not mol.GetBondBetweenAtoms(j, k).IsInRing():
                Chem.rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), i, j, k, l, v)
    mp = Chem.AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    if mp is not None:
        ff = Chem.AllChem.MMFFGetMoleculeForceField(mol, mp)
        if ff is not None:
            if dihedral is not None:
                for (i, j, k, l), v in dihedral.items():
                    ff.MMFFAddTorsionConstraint(i, j, k, l, False, v, v, 4e5)
            ff.Minimize(maxIts=10000)
    return Chem.RemoveHs(mol)


def GetRoundedDihedral(mol: Chem.Mol, dihedral: Tuple[int, int, int, int]) -> int:
    """Return dihedral angle rounded to nearest 60 degrees

    :param mol: RDKit Mol
    :param dihedral: tuple of four atom indices
    :return: rounded dihedral angle
    """
    return (
        round(Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(0), *dihedral) / 60)
        * 60
    )
