from rdkit import Chem
from typing import Optional, List, Set, FrozenSet
from .atom import Atom, GetSeenIndex, SetSeenIndex, IsSeen
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups
from .tokens import DOT, SKIP


def _ringTokens(n, nSkip):
    yield f"[{n}]" if n > 9 else f"{n}"
    yield from iter(SKIP * nSkip)


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


def _removeTrailingDots(t):
    i = len(t) - 1
    while i >= 0:
        if not _isDot(t[i]):
            del t[i + 1 :]
            break
        i -= 1
    return t


def _bondTokens(b, bond):
    if bond is not None:
        yield (b, bond)


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
                        yield from _bondTokens(b, bond)
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
                yield from _bondTokens(b, bond)
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

    t = [
        t[1].asToken(t[0], mol) if type(t) == tuple else t
        for t in list(_getPreTokens())
    ]

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
