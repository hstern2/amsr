from rdkit import Chem

from .atom import GetSeenIndex
from .conf import GetRoundedDihedral
from .tokens import BOND_SYMBOL_FOR_DIHEDRAL, E, Z


def _is_rotatable(b):
    if b.GetBondType() != Chem.rdchem.BondType.SINGLE:
        return False
    if b.GetIsAromatic():
        return False
    if b.GetBeginAtom().GetDegree() == 1 or b.GetEndAtom().GetDegree() == 1:
        return False
    if b.IsInRing():
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()
        return (
            a1.GetHybridization() == Chem.HybridizationType.SP3
            and a2.GetHybridization() == Chem.HybridizationType.SP3
        )
    return True


def _earliestSeenNotIncluding(a, bi, avoid_equiv_terminals=False):
    """Pick dihedral reference neighbor: earliest seen.

    If avoid_equiv_terminals is True and the earliest-seen candidate is a
    terminal (degree 1) with other terminals present, prefer the earliest-
    seen non-terminal instead.  Must match decode.py:_dihedral_ref.
    """
    nbrs = []
    for c in a.GetNeighbors():
        ci = c.GetIdx()
        if ci != bi:
            nbrs.append((GetSeenIndex(c), c.GetDegree(), ci))
    if not nbrs:
        return None
    nbrs.sort()  # by seenIndex
    pick_seen, pick_deg, pick_ci = nbrs[0]
    if avoid_equiv_terminals and pick_deg == 1:
        if sum(1 for _, d, _ in nbrs if d == 1) > 1:
            for _, d, ci in nbrs:
                if d > 1:
                    return ci
    return pick_ci


class Bond:
    def __init__(self, sym=None, isRotatable=False):
        self.sym = sym
        self.isRotatable = isRotatable

    def rdStereo(self):
        if self.sym == E:
            return Chem.BondStereo.STEREOE
        elif self.sym == Z:
            return Chem.BondStereo.STEREOZ
        else:
            return None

    def asToken(self, b):
        m = b.GetOwningMol()
        a1 = b.GetBeginAtom()
        i1 = a1.GetIdx()
        a2 = b.GetEndAtom()
        i2 = a2.GetIdx()
        if self.isRotatable:
            j1 = _earliestSeenNotIncluding(a1, i2, avoid_equiv_terminals=True)
            j2 = _earliestSeenNotIncluding(a2, i1, avoid_equiv_terminals=True)
            return BOND_SYMBOL_FOR_DIHEDRAL[GetRoundedDihedral(m, (j1, i1, i2, j2), 30)]
        n1 = [GetSeenIndex(c) for c in a1.GetNeighbors() if c.GetIdx() != i2]
        n2 = [GetSeenIndex(c) for c in a2.GetNeighbors() if c.GetIdx() != i1]
        flip = False
        for a in b.GetStereoAtoms():
            i = GetSeenIndex(m.GetAtomWithIdx(a))
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
        m = b.GetOwningMol()
        if m.GetNumConformers() > 0 and m.GetConformer().Is3D() and _is_rotatable(b):
            return cls(isRotatable=True)
        return None
