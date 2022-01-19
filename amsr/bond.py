from rdkit import Chem
from .atom import visitedIndex


class Bond:
    def __init__(self, sym=None):
        self.sym = sym

    def rdStereo(self):
        if self.sym == "_":
            return Chem.BondStereo.STEREOE
        elif self.sym == "^":
            return Chem.BondStereo.STEREOZ
        else:
            return None

    def asToken(self, b, atom):
        if self.rdStereo() is None:
            return
        a1 = b.GetBeginAtom()
        i1 = a1.GetIdx()
        a2 = b.GetEndAtom()
        i2 = a2.GetIdx()
        n1 = [visitedIndex(c, atom) for c in a1.GetNeighbors() if c.GetIdx() != i2]
        n2 = [visitedIndex(c, atom) for c in a2.GetNeighbors() if c.GetIdx() != i1]
        flip = False
        for a in b.GetStereoAtoms():
            i = atom[a].visitedIndex
            if i in n1 and i != min(n1):
                flip = not flip
            elif i in n2 and i != min(n2):
                flip = not flip
        if self.sym == "_":
            yield "^" if flip else "_"
        if self.sym == "^":
            yield "_" if flip else "^"

    @classmethod
    def fromRD(cls, b):
        t = b.GetBondType()
        s = b.GetStereo()
        if s == Chem.BondStereo.STEREOE:
            return cls("_")
        if s == Chem.BondStereo.STEREOZ:
            return cls("^")
        return cls()


def bondStereo(bond):
    if bond == "_":
        return Chem.BondStereo.STEREOE
    elif bond == "^":
        return Chem.BondStereo.STEREOZ
    else:
        return None
