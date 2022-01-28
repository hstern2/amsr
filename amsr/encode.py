from rdkit import Chem
from .atom import Atom, GetSeenIndex, SetSeenIndex, IsSeen
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups
from .tokens import DOT, NOP, RemoveTrailingDots, OmitUnneededNOPs


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


def FromMolToTokens(mol):

    atom = [Atom.fromRDAtom(a) for a in mol.GetAtoms()]
    seenBonds = set()
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

    def _getTokens():
        return [
            t[1].asToken(t[0], mol) if type(t) == tuple else t
            for t in list(_getPreTokens())
        ]

    return RemoveTrailingDots(EncodeGroups(OmitUnneededNOPs(_getTokens())))


def FromMol(m):
    return "".join(FromMolToTokens(m))


def FromSmiles(s):
    return FromMol(Chem.MolFromSmiles(s))


def FromSmilesToTokens(s):
    return FromMolToTokens(Chem.MolFromSmiles(s))
