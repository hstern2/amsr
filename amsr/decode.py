from rdkit import Chem
from typing import Tuple, Optional, List, Dict
from re import match, escape
from .atom import Atom
from .bond import Bond
from .pibonds import PiBonds
from .bfs import BFSTree
from .parity import IsEvenParity
from .groups import DecodeGroups
from .tokens import RegExp, SKIP, L_BRACKET, R_BRACKET, DIHEDRAL_FOR_BOND_SYMBOL


def _addBond(mol, atom, i, j, bond, dihedral_for_bond):
    atom[i].addBondTo(atom[j])
    n = mol.AddBond(i, j, Chem.BondType.SINGLE) - 1
    s = bond.rdStereo()
    if s is not None:
        mol.GetBondWithIdx(n).SetStereo(s)
    if bond.sym:
        dihedral_for_bond[n] = DIHEDRAL_FOR_BOND_SYMBOL[bond.sym]


def _addAtom(mol, atom, a):
    atom.append(a)
    mol.AddAtom(a.asRDAtom())


def _bondAtom(mol, atom, a, bond, makeBond, stringent, dihedral_for_bond):
    if not makeBond or not a.canBond():
        _addAtom(mol, atom, a)
        return
    for i in reversed(range(len(atom))):
        if atom[i].canBond() and atom[i].canBondWith(a, stringent):
            _addAtom(mol, atom, a)
            _addBond(mol, atom, i, len(atom) - 1, bond, dihedral_for_bond)
            return
    for i in reversed(range(len(atom))):
        if atom[i].nNeighbors < atom[i].maxNeighbors and atom[i].canBondWith(
            a, stringent
        ):
            _addAtom(mol, atom, a)
            _addBond(mol, atom, i, len(atom) - 1, bond, dihedral_for_bond)
            return


def _saturate(atom):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].isSaturated = True
            return


def _ring(mol, atom, ringStr, bond, stringent, dihedral_for_bond):
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
            if not atom[j].canBond() or not atom[i].canBondWith(atom[j], stringent):
                continue
            if nSkip == 0:
                _addBond(mol, atom, i, j, bond, dihedral_for_bond)
                return
            else:
                nSkip -= 1


def ToMol(
    s: str,
    stringent: Optional[bool] = True,
    dihedral: Optional[Dict[Tuple[int, int, int, int], int]] = None,
) -> Chem.Mol:
    """Convert AMSR to an RDKit Mol

    :param s: AMSR
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :param dihedral: return dictionary of dihedral angles, where keys are indices and values are angles in degrees
    :return: RDKit Mol
    """
    mol = Chem.RWMol()
    atom: List[Atom] = []
    dihedral_for_bond: Dict[int] = {}
    makeBond = False
    for m in RegExp.finditer(DecodeGroups(s)):
        if m.group("ring"):
            _ring(
                mol,
                atom,
                m.group("ring"),
                Bond(m.group("bond")),
                stringent,
                dihedral_for_bond,
            )
        elif m.group("atom"):
            _bondAtom(
                mol,
                atom,
                Atom(m.group("atom")),
                Bond(m.group("bond")),
                makeBond,
                stringent,
                dihedral_for_bond,
            )
            makeBond = True
        elif m.group("saturate"):
            _saturate(atom)
        elif m.group("molsep"):
            makeBond = False
    for a in mol.GetAtoms():
        if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            n = len(a.GetBonds())
            if n < 3:
                a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            elif not IsEvenParity([b.GetIdx() for b in a.GetNeighbors()]):
                a.InvertChirality()
    for b in mol.GetBonds():
        k = b.GetIdx()
        is_EZ = b.GetStereo() in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE)
        is_dihedral = dihedral is not None and k in dihedral_for_bond
        if is_EZ or is_dihedral:
            ai, aj = b.GetBeginAtom(), b.GetEndAtom()
            i, j = ai.GetIdx(), aj.GetIdx()
            ni = [c.GetIdx() for c in ai.GetNeighbors() if c.GetIdx() != j]
            nj = [c.GetIdx() for c in aj.GetNeighbors() if c.GetIdx() != i]
            if len(ni) == 0 or len(nj) == 0:
                b.SetStereo(Chem.BondStereo.STEREONONE)
            else:
                mi = min(ni)
                mj = min(nj)
                if is_EZ:
                    b.SetStereoAtoms(mi, mj)
                if is_dihedral:
                    dihedral[mi, i, j, mj] = dihedral_for_bond[k]
    PiBonds(mol, atom, stringent)
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


def ToSmiles(s: str, stringent: Optional[bool] = True) -> str:
    """Convert AMSR to SMILES

    :param s: AMSR
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :return: SMILES
    """
    return Chem.MolToSmiles(ToMol(s, stringent))
