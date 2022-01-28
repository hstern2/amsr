from rdkit import Chem
from .atom import GetSeenIndex
from .tokens import E, Z


class Bond:
    def __init__(self, sym=None):
        self.sym = sym

    def rdStereo(self):
        if self.sym == E:
            return Chem.BondStereo.STEREOE
        elif self.sym == Z:
            return Chem.BondStereo.STEREOZ
        else:
            return None

    def asToken(self, b, mol):
        a1 = b.GetBeginAtom()
        i1 = a1.GetIdx()
        a2 = b.GetEndAtom()
        i2 = a2.GetIdx()
        n1 = [GetSeenIndex(c) for c in a1.GetNeighbors() if c.GetIdx() != i2]
        n2 = [GetSeenIndex(c) for c in a2.GetNeighbors() if c.GetIdx() != i1]
        flip = False
        for a in b.GetStereoAtoms():
            i = GetSeenIndex(mol.GetAtomWithIdx(a))
            if i in n1 and i != min(n1):
                flip = not flip
            elif i in n2 and i != min(n2):
                flip = not flip
        if self.sym == E:
            return Z if flip else E
        if self.sym == Z:
            return E if flip else Z

    @classmethod
    def fromRD(cls, b):
        s = b.GetStereo()
        if s == Chem.BondStereo.STEREOE:
            return cls(E)
        if s == Chem.BondStereo.STEREOZ:
            return cls(Z)
        return None
