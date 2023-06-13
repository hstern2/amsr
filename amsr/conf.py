from rdkit import Chem


def GetConformer(mol, dihedral=None):
    mol = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol)
    if dihedral is not None:
        for (i, j, k, l), v in dihedral.items():
            if not mol.GetBondBetweenAtoms(j, k).IsInRing():
                Chem.rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), i, j, k, l, v)
    mp = Chem.AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = Chem.AllChem.MMFFGetMoleculeForceField(mol, mp)
    if dihedral is not None:
        for (i, j, k, l), v in dihedral.items():
            ff.MMFFAddTorsionConstraint(i, j, k, l, False, v, v, 4e5)
    ff.Minimize(maxIts=10000)
    return mol
