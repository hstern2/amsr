from rdkit import Chem
from networkx import Graph
from networkx.algorithms import max_weight_matching


def piBond(mol, atom, i, j, order):
    mol.GetBondBetweenAtoms(i, j).SetBondType(
        Chem.BondType.TRIPLE if order == 3 else Chem.BondType.DOUBLE
    )
    atom[i].nPiBonds += order - 1
    atom[j].nPiBonds += order - 1


def incrementBond(mol, atom, i, j):
    b = mol.GetBondBetweenAtoms(i, j)
    t = b.GetBondType()
    if t == Chem.BondType.SINGLE:
        b.SetBondType(Chem.BondType.DOUBLE)
    elif t == Chem.BondType.DOUBLE:
        b.SetBondType(Chem.BondType.TRIPLE)
    atom[i].nPiBonds += 1
    atom[j].nPiBonds += 1


def subgraph(d, condition):
    return {i: {j for j in n if condition(j)} for i, n in d.items() if condition(i)}


def PiBonds(mol, atom):
    def canPiBond(i):
        return atom[i].nAvailablePiBonds() >= 1

    def canMultiplePiBond(i):
        return atom[i].nAvailablePiBonds() >= 2

    def canBePi(i, j):
        return canPiBond(i) and canPiBond(j)

    def canBeTriple(i, j):
        return canMultiplePiBond(i) and canMultiplePiBond(j)

    def singleCoordinate(g, i):
        return g.degree[i] == 1

    def possiblePiBonds():
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if canBePi(i, j):
                yield i, j

    # subgraph of atoms that can make pi bonds
    g = Graph(possiblePiBonds())
    # single coordinate
    done = False
    while not done:
        done = True
        for i in g.nodes:
            if singleCoordinate(g, i):
                for j in g.neighbors(i):
                    if canBeTriple(i, j):
                        piBond(mol, atom, i, j, 3)
                    elif canBePi(i, j):
                        piBond(mol, atom, i, j, 2)
                    else:
                        continue
                    done = False
        if not done:
            g = g.subgraph(filter(canPiBond, g.nodes))
    # isolated sp centers
    for i in g.nodes:
        if canMultiplePiBond(i):
            for j in g.neighbors(i):
                if canMultiplePiBond(j):
                    break
            else:
                for j in g.neighbors(i):
                    if canBePi(i, j):
                        piBond(mol, atom, i, j, 2)
    # remaining pi bonds
    done = False
    while not done:
        done = True
        g = g.subgraph(filter(canPiBond, g.nodes))
        for i, j in max_weight_matching(g):
            incrementBond(mol, atom, i, j)
            done = False
