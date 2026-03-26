from rdkit import Chem

from .atom import GetSeenIndex
from .conf import GetRoundedDihedral
from .tokens import BOND_SYMBOL_FOR_DIHEDRAL, E, Z


def _kekulized_bond_type(b):
    """Return the Kekulized bond type, caching the Kekulized mol on first call."""
    mol = b.GetOwningMol()
    if not hasattr(mol, "_kekulized"):
        mol_k = Chem.RWMol(mol)
        try:
            Chem.Kekulize(mol_k, clearAromaticFlags=False)
        except Exception:
            mol_k = mol
        mol._kekulized = mol_k
    return mol._kekulized.GetBondWithIdx(b.GetIdx()).GetBondType()


def _is_rotatable(b):
    if b.GetBondType() != Chem.rdchem.BondType.SINGLE:
        # Aromatic bonds in large (>6) rings may be Kekulized single bonds
        # that can rotate (e.g. tropone).
        if b.GetIsAromatic() and b.IsInRing():
            ri = b.GetOwningMol().GetRingInfo()
            sizes = set(ri.BondRingSizes(b.GetIdx()))
            if all(s > 6 for s in sizes) and _kekulized_bond_type(b) == Chem.BondType.SINGLE:
                return True
        return False
    if b.GetIsAromatic():
        return False
    if b.GetBeginAtom().GetDegree() == 1 or b.GetEndAtom().GetDegree() == 1:
        return False
    if b.IsInRing():
        a1 = b.GetBeginAtom()
        a2 = b.GetEndAtom()
        h1, h2 = a1.GetHybridization(), a2.GetHybridization()
        # At least one SP3 endpoint — excludes rigid SP2-SP2 bonds
        # but allows SP3-SP2 ring bonds (e.g. C-N in piperidine lactams).
        return Chem.HybridizationType.SP3 in (h1, h2)
    return True


def _earliestSeenNotIncluding(a, bi, avoid_equiv_terminals=False):
    """Pick dihedral reference neighbor: earliest seen.

    If avoid_equiv_terminals is True and the earliest-seen candidate is a
    terminal (degree 1) with other terminals present, prefer the earliest-
    seen non-terminal instead.  When all candidates are terminals, break
    ties by highest atomic number to stay consistent with decode.py
    regardless of group-expansion atom ordering.
    Must match decode.py:_dihedral_ref.
    """
    mol = a.GetOwningMol()
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
        n_terminals = sum(1 for _, d, _ in nbrs if d == 1)
        if n_terminals > 1:
            for _, d, ci in nbrs:
                if d > 1:
                    return ci
            # All terminals: pick by highest atomic number, then highest
            # bond order to parent (stable across group-expansion
            # reorderings where atom index may differ).
            ai = a.GetIdx()
            terminals = [
                (
                    -mol.GetAtomWithIdx(ci).GetAtomicNum(),
                    -int(mol.GetBondBetweenAtoms(ai, ci).GetBondTypeAsDouble() * 10),
                    seen,
                    ci,
                )
                for seen, d, ci in nbrs
                if d == 1
            ]
            terminals.sort()
            return terminals[0][3]
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
            # Order atoms so a1 is the earlier-seen (parent in DFS),
            # matching the decoder's convention for dihedral references.
            if GetSeenIndex(a1) > GetSeenIndex(a2):
                a1, a2 = a2, a1
                i1, i2 = i2, i1
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
