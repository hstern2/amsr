from rdkit import Chem
from typing import Optional, List, Set, FrozenSet
from .atom import Atom, GetSeenIndex, SetSeenIndex, IsSeen, UnSee
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups
from .tokens import DOT, SKIP, MOLSEP
from random import shuffle


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


def FromMolToTokens(
    mol: Chem.Mol,
    useGroups: Optional[bool] = True,
    stringent: Optional[bool] = True,
    randomize: Optional[bool] = False,
    canonical: Optional[bool] = False,
    useStereo: Optional[bool] = True,
) -> List[str]:
    """Convert RDKit Mol to list of AMSR tokens

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations (default: True)
    :param stringent: try to exclude unstable or synthetically inaccessible molecules (default: True)
    :param randomize: randomize order of graph traversal (default: False)
    :param canonical: canonical order of graph traversal (default: False)
    :param useStereo: encode stereochemistry (default: True)
    :return: list of AMSR tokens
    """

    assert not (randomize and canonical)

    if not useStereo:
        Chem.RemoveStereochemistry(mol)

    if randomize:
        i = list(range(mol.GetNumAtoms()))
        shuffle(i)
        mol = Chem.RenumberAtoms(mol, i)
    elif canonical:
        ranks = list(Chem.CanonicalRankAtoms(mol, includeChirality=useStereo))
        i = sorted(range(len(ranks)), key=lambda x: ranks[x])
        mol = Chem.RenumberAtoms(mol, i)

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
                        ai.addBondTo(aj)
                        break
                    elif atom[k].canBond() and ai.canBondWith(
                        atom[k], stringent=stringent
                    ):
                        nSkip += 1
            else:  # new atom
                seenBonds.add(ij)
                ai.addBondTo(aj)
                yield from _bondTokens(b, bond)
                yield (c, aj)
                yield from _search(c)
        # end loop over bonds
        if ai.canBond():
            ai.isSaturated = True
            yield DOT

    def _getPreTokens():
        _ = False
        for i, a in enumerate(mol.GetAtoms()):
            if not IsSeen(a):
                if _:
                    yield MOLSEP
                yield a, atom[i]
                yield from _search(a)
                _ = True

    t = [t[1].asToken(t[0]) if type(t) == tuple else t for t in list(_getPreTokens())]

    for a in mol.GetAtoms():
        UnSee(a)

    if useGroups:
        t = EncodeGroups(t)

    return _removeTrailingDots(t)


def FromMol(
    mol: Chem.Mol,
    useGroups: Optional[bool] = True,
    stringent: Optional[bool] = True,
    randomize: Optional[bool] = False,
    canonical: Optional[bool] = False,
    useStereo: Optional[bool] = True,
) -> str:
    """Convert RDKit Mol to AMSR

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations (default: True)
    :param stringent: try to exclude unstable or synthetically inaccessible molecules (default: True)
    :param randomize: randomize order of graph traversal (default: False)
    :param canonical: canonical order of graph traversal (default: False)
    :param useStereo: encode stereochemistry (default: True)
    :return: list of AMSR tokens
    """
    return "".join(
        FromMolToTokens(
            mol,
            useGroups=useGroups,
            stringent=stringent,
            randomize=randomize,
            canonical=canonical,
            useStereo=useStereo,
        )
    )


def FromSmiles(
    s: str,
    useGroups: Optional[bool] = True,
    stringent: Optional[bool] = True,
    randomize: Optional[bool] = False,
    canonical: Optional[bool] = False,
    useStereo: Optional[bool] = True,
) -> str:
    """Convert SMILES to AMSR

    :param s: SMILES
    :param useGroups: use group symbols/abbreviations
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :param randomize: randomize order of graph traversal
    :param canonical: canonical order of graph traversal (default: False)
    :param useStereo: encode stereochemistry (default: True)
    :return: AMSR
    """
    return FromMol(
        Chem.MolFromSmiles(s),
        useGroups=useGroups,
        stringent=stringent,
        randomize=randomize,
        canonical=canonical,
        useStereo=useStereo,
    )


def FromSmilesToTokens(
    s,
    useGroups: Optional[bool] = True,
    stringent: Optional[bool] = True,
    randomize: Optional[bool] = False,
    canonical: Optional[bool] = False,
    useStereo: Optional[bool] = True,
) -> List[str]:
    """Convert SMILES to list of AMSR tokens

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :param randomize: randomize order of graph traversal
    :param canonical: canonical order of graph traversal (default: False)
    :param useStereo: encode stereochemistry (default: True)
    :return: AMSR
    """
    return FromMolToTokens(
        Chem.MolFromSmiles(s),
        useGroups=useGroups,
        stringent=stringent,
        randomize=randomize,
        canonical=canonical,
        useStereo=useStereo,
    )
