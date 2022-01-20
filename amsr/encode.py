from rdkit import Chem
from itertools import repeat
from .atom import Atom, visitedIndex
from .bfs import BFSFind
from .bond import Bond


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
    nAtoms = mol.GetNumAtoms()
    if nAtoms == 0:
        return
    atom = [Atom.fromRD(a) for a in mol.GetAtoms()]
    seenBonds = set()
    nSeenAtoms = 0
    nBonds = mol.GetNumBonds()

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
                yield (b, bond)
                yield (c, sj)
                yield from search(c)
        # end loop over bonds
        if si.canBond() and (len(seenBonds) < nBonds or nSeenAtoms < nAtoms):
            si.isSaturated = True
            yield "."

    def getTokens1():
        for a in mol.GetAtoms():
            if not seen(a):
                yield a, atom[a.GetIdx()]
                yield from search(a)

    def getTokens2():
        for t in list(getTokens1()):
            yield from t[1].asToken(t[0], atom) if type(t) == tuple else t

    def isDot(t):
        return t == "."

    def isDigit(t):
        return t in ("3", "4", "5", "6")

    def isDotOrDigit(t):
        return isDot(t) or isDigit(t)

    def needsSpace(prv, nxt):
        return (isDigit(prv) and isDotOrDigit(nxt)) or (isDot(prv) and isDot(nxt))

    def getTokens3():
        tok = list(getTokens2())
        ntok = len(tok)
        for i, t in enumerate(tok):
            if t != " " or (0 < i < ntok - 1 and needsSpace(tok[i - 1], tok[i + 1])):
                yield t

    tok = list(getTokens3())
    ntok = len(tok)
    PHENYL = list("cccccc6 .....")
    N_PHENYL = len(PHENYL)
    i = 0
    while i < ntok:
        j = max(7, min(N_PHENYL, ntok - i))
        if tok[i : i + j] == PHENYL[:j]:
            yield "[Ph]"
            i += N_PHENYL
        else:
            yield tok[i]
            i += 1


def FromMol(m):
    return "".join(FromMolToTokens(m))


def FromSmiles(s):
    return FromMol(Chem.MolFromSmiles(s))


def FromSmilesToTokens(s):
    return FromMolToTokens(Chem.MolFromSmiles(s))
