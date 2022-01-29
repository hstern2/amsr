from rdkit import Chem
from typing import Optional, List
from re import match
from .atom import Atom
from .bond import Bond
from .pibonds import PiBonds
from .bfs import BFSTree
from .parity import IsEvenParity
from .groups import DecodeGroups
from .tokens import RegExp


def _addBond(mol, atom, i, j, bond):
    atom[i].nNeighbors += 1
    atom[j].nNeighbors += 1
    n = mol.AddBond(i, j, Chem.BondType.SINGLE)
    s = bond.rdStereo()
    if s is not None:
        mol.GetBondWithIdx(n - 1).SetStereo(s)


def _addAtom(mol, atom, a, bond):
    atom.append(a)
    mol.AddAtom(a.asRDAtom())
    if a.canBond():
        j = len(atom) - 1
        for i in reversed(range(j)):
            if atom[i].canBond():
                _addBond(mol, atom, i, j, bond)
                return


def _saturate(atom):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].isSaturated = True
            return


def _addDanglingBond(atom, dangling):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].nNeighbors += 1
            dangling.append(i)
            return


def _connectToDanglingBond(atom, dangling):
    if not dangling:
        return
    i = dangling.pop()
    atom[i].nNeighbors -= 1
    for a in atom[i + 1 :]:
        a.isSaturated = True


def _ring(mol, atom, ring, bond):
    m = match(r"(\d+)(\.*)", ring)
    n = sum(map(int, m.group(1)))
    nSkip = len(m.group(2))
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            for j in BFSTree(mol.GetAtomWithIdx(i), n - 1):
                if atom[j].canBond():
                    if nSkip == 0:
                        _addBond(mol, atom, i, j, bond)
                        return
                    else:
                        nSkip -= 1


def ToMol(s: str) -> Chem.Mol:
    """Convert AMSR to an RDKit Mol

    :param s: AMSR
    :return: RDKit Mol
    """
    mol = Chem.RWMol()
    atom = []
    dangling = []
    for m in RegExp.finditer(DecodeGroups(s)):
        if m.group("ring"):
            _ring(mol, atom, m.group("ring"), Bond(m.group("bond")))
        elif m.group("atom"):
            _addAtom(mol, atom, Atom(m.group("atom")), Bond(m.group("bond")))
        elif m.group("saturate"):
            _saturate(atom)
        elif m.group("dangling"):
            _addDanglingBond(atom, dangling)
        elif m.group("ampersand"):
            _connectToDanglingBond(atom, dangling)
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
    PiBonds(mol, atom)
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


def ToSmiles(s: str) -> str:
    """Convert AMSR to SMILES

    :param s: AMSR
    :return: SMILES
    """
    return Chem.MolToSmiles(ToMol(s))
