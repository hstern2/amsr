from rdkit import Chem
from re import finditer
from .atom import Atom
from .bond import Bond
from .pibonds import PiBonds
from .bfs import BFSTree
from .parity import IsEvenParity
from .groups import DecodeGroups


def AddBond(mol, atom, i, j, bond):
    atom[i].nNeighbors += 1
    atom[j].nNeighbors += 1
    n = mol.AddBond(i, j, Chem.BondType.SINGLE)
    s = bond.rdStereo()
    if s is not None:
        mol.GetBondWithIdx(n - 1).SetStereo(s)


def AddAtom(mol, atom, a, bond):
    atom.append(a)
    mol.AddAtom(a.asRDAtom())
    if a.canBond():
        j = len(atom) - 1
        for i in reversed(range(j)):
            if atom[i].canBond():
                AddBond(mol, atom, i, j, bond)
                return


def Saturate(atom):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].isSaturated = True
            return


def AddDanglingBond(atom, dangling):
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            atom[i].nNeighbors += 1
            dangling.append(i)
            return


def ConnectToDanglingBond(atom, dangling):
    if not dangling:
        return
    i = dangling.pop()
    atom[i].nNeighbors -= 1
    for a in atom[i + 1 :]:
        a.isSaturated = True


def Ring(mol, atom, ring, skip, bond):
    n = sum(map(int, ring))
    nSkip = len(skip) if skip else 0
    for i in reversed(range(len(atom))):
        if atom[i].canBond():
            for j in BFSTree(mol.GetAtomWithIdx(i), n - 1):
                if atom[j].canBond():
                    if nSkip == 0:
                        AddBond(mol, atom, i, j, bond)
                        return
                    else:
                        nSkip -= 1


def ToMol(s, activeAtom=None):
    pbond = r"(?P<bond>[\^\_])"
    c = r"[\+\-\:\!\*\(\)]*"
    patom = r"(?P<atom>\[\d*[A-Za-z][a-z]?" + c + r"\]|[A-Za-z]" + c + r")"
    pring = r"(?P<ring>[3-6]+)"
    pskip = r"(?P<skip>\.*)"
    pdot = r"(?P<dot>\.)"
    pdangling = r"(?P<dangling>\[\])"
    pampersand = r"(?P<ampersand>\&)"
    mol = Chem.RWMol()
    atom = []
    dangling = []
    s = DecodeGroups(s)
    for m in finditer(
        f"({pbond}?({patom}|({pring}{pskip})))|{pdot}|{pdangling}|{pampersand}", s
    ):
        if m.group("ring"):
            Ring(mol, atom, m.group("ring"), m.group("skip"), Bond(m.group("bond")))
        elif m.group("atom"):
            AddAtom(mol, atom, Atom(m.group("atom")), Bond(m.group("bond")))
        elif m.group("dot"):
            Saturate(atom)
        elif m.group("dangling"):
            AddDanglingBond(atom, dangling)
        elif m.group("ampersand"):
            ConnectToDanglingBond(atom, dangling)
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
    if activeAtom is not None:
        for i in reversed(range(len(atom))):
            if atom[i].canBond():
                activeAtom.append(i)
                break
    return mol


def ToSmiles(s):
    return Chem.MolToSmiles(ToMol(s))
