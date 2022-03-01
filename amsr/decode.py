from rdkit import Chem
from typing import List
from re import match, escape
from .atom import Atom
from .bond import Bond
from .pibonds import PiBonds
from .bfs import BFSTree
from .parity import IsEvenParity
from .groups import DecodeGroups
from .tokens import RegExp, SKIP, L_BRACKET, R_BRACKET


def _addBond(mol, atom, i, j, bond):
    atom[i].nNeighbors += 1
    atom[j].nNeighbors += 1
    n = mol.AddBond(i, j, Chem.BondType.SINGLE)
    s = bond.rdStereo()
    if s is not None:
        mol.GetBondWithIdx(n - 1).SetStereo(s)


def _addAtom(mol, atom, a, bond, contiguous=False):
    atom.append(a)
    mol.AddAtom(a.asRDAtom())
    if not a.canBond():
        return
    j = len(atom) - 1
    for i in reversed(range(j)):
        if atom[i].canBond():
            _addBond(mol, atom, i, j, bond)
            return
    if not contiguous:
        return
    for i in reversed(range(j)):
        if atom[i].nNeighbors < atom[i].maxNeighbors:
            _addBond(mol, atom, i, j, bond)
            return


def _saturate(atom):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].isSaturated = True
            return


def _ring(mol, atom, ringStr, bond):
    m = match(
        f"{escape(L_BRACKET)}?([0-9]+){escape(R_BRACKET)}?({escape(SKIP)}*)", ringStr
    )
    n = int(m.group(1))
    if n < 3:
        return
    nSkip = len(m.group(2))
    for i in reversed(range(len(atom))):
        if not atom[i].canBond():
            continue
        for node in BFSTree(mol.GetAtomWithIdx(i), n - 1):
            j = node.name.GetIdx()
            if not atom[j].canBond():
                continue
            if nSkip == 0:
                _addBond(mol, atom, i, j, bond)
                return
            else:
                nSkip -= 1


def ToMol(s: str, contiguous: bool = False, useFilters: bool = True) -> Chem.Mol:
    """Convert AMSR to an RDKit Mol

    :param s: AMSR
    :param contiguous: make second pass (considering "saturated" atoms) to form bonds, rather than start new molecule
    :param useFilters: use a subset of filters from J. Chem. Inf. Model. 52, 2864 (2012) to exclude unstable or synthetically inaccessible molecules
    :return: RDKit Mol
    """
    mol = Chem.RWMol()
    atom: List[Atom] = []
    for m in RegExp.finditer(DecodeGroups(s)):
        if m.group("ring"):
            _ring(mol, atom, m.group("ring"), Bond(m.group("bond")))
        elif m.group("atom"):
            _addAtom(
                mol, atom, Atom(m.group("atom")), Bond(m.group("bond")), contiguous
            )
        elif m.group("saturate"):
            _saturate(atom)
    for a in mol.GetAtoms():
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            n = len(a.GetBonds())
            if n < 3:
                a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            elif not IsEvenParity([b.GetIdx() for b in a.GetNeighbors()]):
                a.InvertChirality()
    for b in mol.GetBonds():
        if b.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
            ai, aj = b.GetBeginAtom(), b.GetEndAtom()
            i, j = ai.GetIdx(), aj.GetIdx()
            ni = [c.GetIdx() for c in ai.GetNeighbors() if c.GetIdx() != j]
            nj = [c.GetIdx() for c in aj.GetNeighbors() if c.GetIdx() != i]
            if len(ni) == 0 or len(nj) == 0:
                b.SetStereo(Chem.BondStereo.STEREONONE)
            else:
                b.SetStereoAtoms(min(ni), min(nj))
    PiBonds(mol, atom, useFilters)
    for i, a in enumerate(atom):
        if a.bangs > 0 and a.canBond():
            mol.GetAtomWithIdx(i).SetNumExplicitHs(a.maxNeighbors - a.nNeighbors)
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol)
    for a in reversed(mol.GetAtoms()):
        if atom[a.GetIdx()].canBond():
            a.SetBoolProp("_active", True)
            break

    return mol


def ToSmiles(s: str, contiguous: bool = False, useFilters: bool = True) -> str:
    """Convert AMSR to SMILES

    :param s: AMSR
    :param contiguous: make second pass (considering "saturated" atoms) to form bonds, rather than start new molecule
    :param useFilters: use a subset of filters from J. Chem. Inf. Model. 52, 2864 (2012) to exclude unstable or synthetically inaccessible molecules
    :return: SMILES
    """
    return Chem.MolToSmiles(ToMol(s, contiguous, useFilters))
