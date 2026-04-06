import math
from random import shuffle
from typing import Optional

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

from .atom import Atom, GetSeenIndex, IsSeen, SetSeenIndex, UnSee
from .bfs import BFSFind
from .bond import Bond
from .groups import EncodeGroups
from .tokens import DOT, MOLSEP, SKIP


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
) -> list[str]:
    """Convert RDKit Mol to list of AMSR tokens

    :param mol: RDKit Mol
    :param useGroups: use group symbols/abbreviations (default: True)
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
        (default: True)
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

    # Assign pseudo-chirality from 3D geometry for centers whose
    # substituents are graph-equivalent.  RDKit won't mark these as
    # chiral, but the 3D arrangement matters for conformer reconstruction.
    # Skip aromatic atoms (always planar) and atoms whose neighbors are
    # all equivalent (e.g. neopentane).  Use out-of-plane distance to
    # distinguish genuinely pyramidal centers from nearly-planar ones.
    if useStereo and mol.GetNumConformers() > 0 and mol.GetConformer().Is3D():
        conf = mol.GetConformer()
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        for a in mol.GetAtoms():
            if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            if a.GetIsAromatic() or a.GetDegree() < 3:
                continue
            hyb = a.GetHybridization()
            if hyb not in (Chem.HybridizationType.SP3, Chem.HybridizationType.SP2):
                continue
            nbrs = [n.GetIdx() for n in a.GetNeighbors()]
            if len(set(ranks[n] for n in nbrs)) < 2:
                continue
            pc = conf.GetAtomPosition(a.GetIdx())
            p0 = conf.GetAtomPosition(nbrs[0])
            p1 = conf.GetAtomPosition(nbrs[1])
            p2 = conf.GetAtomPosition(nbrs[2])
            v0 = (p0.x - pc.x, p0.y - pc.y, p0.z - pc.z)
            v1 = (p1.x - pc.x, p1.y - pc.y, p1.z - pc.z)
            v2 = (p2.x - pc.x, p2.y - pc.y, p2.z - pc.z)
            # Scalar triple product = signed volume of parallelepiped
            vol = (
                v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
                + v0[1] * (v1[2] * v2[0] - v1[0] * v2[2])
                + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0])
            )
            # Out-of-plane distance: |vol| / |cross(v1-v0, v2-v0)|
            cx = (v1[1] - v0[1]) * (v2[2] - v0[2]) - (v1[2] - v0[2]) * (v2[1] - v0[1])
            cy = (v1[2] - v0[2]) * (v2[0] - v0[0]) - (v1[0] - v0[0]) * (v2[2] - v0[2])
            cz = (v1[0] - v0[0]) * (v2[1] - v0[1]) - (v1[1] - v0[1]) * (v2[0] - v0[0])
            cross_mag = math.sqrt(cx * cx + cy * cy + cz * cz)
            if cross_mag < 1e-10:
                continue
            oop_dist = abs(vol) / cross_mag
            if oop_dist < 0.15:
                continue
            # Non-junction amide/carbamate N is planar by resonance;
            # only ring-junction N can be forced pyramidal by strain.
            if a.GetSymbol() == "N" and hyb == Chem.HybridizationType.SP2:
                ri = mol.GetRingInfo()
                is_junction = sum(1 for r in ri.AtomRings() if a.GetIdx() in r) > 1
                if not is_junction:
                    is_amide = any(
                        nb.GetSymbol() == "C"
                        and any(
                            mol.GetBondBetweenAtoms(nb.GetIdx(), nb2.GetIdx()).GetBondTypeAsDouble()
                            == 2
                            and nb2.GetSymbol() == "O"
                            for nb2 in nb.GetNeighbors()
                            if nb2.GetIdx() != a.GetIdx()
                        )
                        for nb in a.GetNeighbors()
                    )
                    if is_amide:
                        continue
            if vol > 0:
                a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif vol < 0:
                a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        # Assign pseudo-E/Z on double bonds whose substituents are
        # graph-equivalent (e.g. guanidine C=N with two identical
        # NH-cyclohexyl groups).  RDKit won't assign E/Z here, but
        # the 3D arrangement matters for conformer reconstruction.
        for b in mol.GetBonds():
            if b.GetBondTypeAsDouble() < 1.5:
                continue
            if b.GetStereo() != Chem.BondStereo.STEREONONE:
                continue
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            ni = [c.GetIdx() for c in b.GetBeginAtom().GetNeighbors() if c.GetIdx() != j]
            nj = [c.GetIdx() for c in b.GetEndAtom().GetNeighbors() if c.GetIdx() != i]
            if not ni or not nj:
                continue
            has_equiv = (len(ni) >= 2 and len(set(ranks[n] for n in ni)) < len(ni)) or (
                len(nj) >= 2 and len(set(ranks[n] for n in nj)) < len(nj)
            )
            if not has_equiv:
                continue
            ri, rj = min(ni), min(nj)
            dih = rdMolTransforms.GetDihedralDeg(conf, ri, i, j, rj)
            if abs(dih) < 90:
                b.SetStereo(Chem.BondStereo.STEREOZ)
            else:
                b.SetStereo(Chem.BondStereo.STEREOE)
            b.SetStereoAtoms(ri, rj)

    atom = [Atom.fromRDAtom(a) for a in mol.GetAtoms()]
    seenBonds: set[frozenset[int]] = set()
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
                    elif atom[k].canBond() and ai.canBondWith(atom[k], stringent=stringent):
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

    t = [t[1].asToken(t[0]) if isinstance(t, tuple) else t for t in list(_getPreTokens())]

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
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
        (default: True)
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
) -> list[str]:
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
