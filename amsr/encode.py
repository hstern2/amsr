from rdkit import Chem
from typing import Optional, List, Set, FrozenSet
from .atom import Atom, GetSeenIndex, SetSeenIndex, IsSeen
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups
from .tokens import DOT, NOP, RING_DIGITS


def _partition3456(n):
    if n <= 6:
        yield n
    elif n <= 8:
        yield from _partition3456(n - 3)
        yield 3
    else:
        yield 6
        yield from _partition3456(n - 6)


def _ringTokens(n, nSkip):
    yield from map(str, _partition3456(n))
    yield from iter(DOT * nSkip)
    yield NOP


def _searchOrder(b, a):
    # 1. seen atoms before unseen atoms (i.e. rings)
    # 2. aromatic bonds first
    # 3. small rings before larger (for seen) .. otherwise atom index (for unseen)
    c = b.GetOtherAtom(a)
    isSeen = IsSeen(c)
    return (
        not isSeen,
        not b.GetIsAromatic(),
        GetSeenIndex(a) - GetSeenIndex(c) if isSeen else c.GetIdx(),
    )


def _isDot(t):
    return t == DOT


def _isRingDigit(t):
    return t in RING_DIGITS


def _isDotOrRingDigit(t):
    return _isDot(t) or _isRingDigit(t)


def _needsNOP(prv, nxt):
    return (_isRingDigit(prv) and _isDotOrRingDigit(nxt)) or (
        _isDot(prv) and _isDot(nxt)
    )


def _removeTrailingDots(t):
    i = len(t) - 1
    while i >= 0:
        if not _isDot(t[i]):
            if not _isRingDigit(t[i]):
                del t[i + 1 :]
            break
        i -= 1
    return t


def _omitUnneededNOPs(t):
    n = len(t)
    while n > 0 and t[n - 1] == NOP:
        n -= 1
    return [
        t[i]
        for i in range(n)
        if t[i] != NOP or (0 < i < n - 1 and _needsNOP(t[i - 1], t[i + 1]))
    ]


def FromMolToTokens(mol: Chem.Mol, useGroups: Optional[bool] = True) -> List[str]:
    """Convert RDKit Mol to list of AMSR tokens

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations
    :return: list of AMSR tokens
    """

    atom = [Atom.fromRDAtom(a) for a in mol.GetAtoms()]
    seenBonds: Set[FrozenSet[int]] = set()
    nSeenAtoms = 0

    def _search(a):
        nonlocal nSeenAtoms
        i = a.GetIdx()
        ai = atom[i]
        SetSeenIndex(a, nSeenAtoms)
        nSeenAtoms += 1
        for b in sorted(a.GetBonds(), key=lambda b: _searchOrder(b, a)):
            c = b.GetOtherAtom(a)
            j = c.GetIdx()
            ij = frozenset([i, j])
            if ij in seenBonds:
                continue
            aj = atom[j]
            bond = Bond.fromRD(b)
            if IsSeen(c):  # ring
                nSkip = 0
                for k, depth in BFSFind(a, j, seenBonds):
                    if k == j:
                        if bond is not None:
                            yield (b, bond)
                        yield from _ringTokens(depth + 1, nSkip)
                        seenBonds.add(ij)
                        ai.nNeighbors += 1
                        aj.nNeighbors += 1
                        break
                    elif atom[k].canBond():
                        nSkip += 1
            else:  # new atom
                seenBonds.add(ij)
                ai.nNeighbors += 1
                aj.nNeighbors += 1
                if bond is not None:
                    yield (b, bond)
                yield (c, aj)
                yield from _search(c)
        # end loop over bonds
        if ai.canBond():
            ai.isSaturated = True
            yield DOT

    def _getPreTokens():
        for i, a in enumerate(mol.GetAtoms()):
            if not IsSeen(a):
                yield a, atom[i]
                yield from _search(a)

    t = _omitUnneededNOPs(
        [
            t[1].asToken(t[0], mol) if type(t) == tuple else t
            for t in list(_getPreTokens())
        ]
    )
    if useGroups:
        t = EncodeGroups(t)
    return _removeTrailingDots(t)


def FromMol(mol: Chem.Mol, useGroups: Optional[bool] = True) -> str:
    """Convert RDKit Mol to AMSR

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations
    :return: AMSR
    """
    return "".join(FromMolToTokens(mol, useGroups=useGroups))


def FromSmiles(s: str, useGroups: Optional[bool] = True) -> str:
    """Convert SMILES to AMSR

    :param s: SMILES
    :param useGroups: use group symbols/abbreviations
    :return: AMSR
    """
    return FromMol(Chem.MolFromSmiles(s), useGroups=useGroups)


def FromSmilesToTokens(s, useGroups: Optional[bool] = True) -> List[str]:
    """Convert SMILES to list of AMSR tokens

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations
    :return: AMSR
    """
    return FromMolToTokens(Chem.MolFromSmiles(s), useGroups=useGroups)
