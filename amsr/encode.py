from rdkit import Chem
from itertools import repeat
from .atom import Atom, visitedIndex
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups


def Partition3456(n):
    if n <= 6:
        yield str(n)
    elif n <= 8:
        yield from Partition3456(n - 3)
        yield 3
    else:
        yield 6
        yield from Partition3456(n - 6)


def RingTokens(n, nSkip):
    for _ in map(str, Partition3456(n)):
        yield _
    for _ in repeat(None, nSkip):
        yield "."
    yield " "


def FromMolToTokens(mol):

    atom = [Atom.fromRD(a) for a in mol.GetAtoms()]
    seenBonds = set()
    nSeenAtoms = 0

    def byVisitedIndex(n):
        return visitedIndex(n.name, atom)

    def seen(a):
        return visitedIndex(a, atom) is not None

    def searchOrder(b, a):
        # 1. seen atoms before unseen atoms (i.e. rings)
        # 2. small rings before larger (for seen) .. otherwise atom index (for unseen)
        c = b.GetOtherAtom(a)
        isSeen = seen(c)
        return (
            not isSeen,
            visitedIndex(a, atom) - visitedIndex(c, atom) if isSeen else c.GetIdx(),
        )

    def search(a):
        nonlocal nSeenAtoms
        i = a.GetIdx()
        si = atom[i]
        si.visitedIndex = nSeenAtoms
        nSeenAtoms += 1
        for b in sorted(a.GetBonds(), key=lambda b: searchOrder(b, a)):
            c = b.GetOtherAtom(a)
            j = c.GetIdx()
            ij = frozenset([i, j])
            if ij in seenBonds:
                continue
            sj = atom[j]
            bond = Bond.fromRD(b)
            if seen(c):  # ring
                nSkip = 0
                for n in sorted(
                    BFSFind(a, j, seenBonds), key=byVisitedIndex, reverse=True
                ):
                    k = n.name.GetIdx()
                    if k == j:
                        if bond is not None:
                            yield (b, bond)
                        yield from RingTokens(n.depth + 1, nSkip)
                        seenBonds.add(ij)
                        si.nNeighbors += 1
                        sj.nNeighbors += 1
                        break
                    elif atom[k].canBond():
                        nSkip += 1
            else:  # new atom
                seenBonds.add(ij)
                si.nNeighbors += 1
                sj.nNeighbors += 1
                if bond is not None:
                    yield (b, bond)
                yield (c, sj)
                yield from search(c)
        # end loop over bonds
        if si.canBond():
            si.isSaturated = True
            yield "."

    def getPreTokens():
        for a in mol.GetAtoms():
            if not seen(a):
                yield a, atom[a.GetIdx()]
                yield from search(a)

    def getTokens():
        return [
            t[1].asToken(t[0], atom) if type(t) == tuple else t
            for t in list(getPreTokens())
        ]

    def isDot(t):
        return t == "."

    def isDigit(t):
        return t in ("3", "4", "5", "6")

    def isDotOrDigit(t):
        return isDot(t) or isDigit(t)

    def needsSpace(prv, nxt):
        return (isDigit(prv) and isDotOrDigit(nxt)) or (isDot(prv) and isDot(nxt))

    def omitUnneededSpaces(tok):
        ntok = len(tok)
        return [
            t
            for i, t in enumerate(tok)
            if t != " " or (0 < i < ntok - 1 and needsSpace(tok[i - 1], tok[i + 1]))
        ]

    def removeTrailingDots(tok):
        i = len(tok) - 1
        while i >= 0:
            if tok[i] != ".":
                if not isDigit(tok[i]):
                    del tok[i + 1 :]
                break
            i -= 1
        return tok

    return removeTrailingDots(EncodeGroups(omitUnneededSpaces(getTokens())))


def FromMol(m):
    return "".join(FromMolToTokens(m))


def FromSmiles(s):
    return FromMol(Chem.MolFromSmiles(s))


def FromSmilesToTokens(s):
    return FromMolToTokens(Chem.MolFromSmiles(s))
