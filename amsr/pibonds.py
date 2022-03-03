from rdkit import Chem
from networkx import Graph
from networkx.algorithms import max_weight_matching
from .ringinfo import BridgeheadAtoms, RingAtoms


class PiBonds:
    def piBond(self, i, j, order):
        self.mol.GetBondBetweenAtoms(i, j).SetBondType(
            Chem.BondType.TRIPLE if order == 3 else Chem.BondType.DOUBLE
        )
        self.atom[i].nPiBonds += order - 1
        self.atom[j].nPiBonds += order - 1

    def incrementBond(self, i, j):
        b = self.mol.GetBondBetweenAtoms(i, j)
        t = b.GetBondType()
        if t == Chem.BondType.SINGLE:
            b.SetBondType(Chem.BondType.DOUBLE)
        elif t == Chem.BondType.DOUBLE:
            b.SetBondType(Chem.BondType.TRIPLE)
        self.atom[i].nPiBonds += 1
        self.atom[j].nPiBonds += 1

    def canPiBond(self, i):
        return self.atom[i].nAvailablePiBonds() >= 1 and i not in self.excluded

    def canMultiplePiBond(self, i):
        return self.atom[i].nAvailablePiBonds() >= 2 and i not in self.excluded

    def canBePi(self, i, j):
        return self.canPiBond(i) and self.canPiBond(j)

    def canBeTriple(self, i, j):
        return self.canMultiplePiBond(i) and self.canMultiplePiBond(j)

    def possiblePiBonds(self):
        for b in self.mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if self.canBePi(i, j):
                yield i, j

    def isCarbon(self, i):
        return self.atom[i].isCarbon()

    def reduceGraph(self):
        self.graph = self.graph.subgraph(
            [i for i in self.graph.nodes if self.canPiBond(i)]
        )

    def singleCoordinate(self, i):
        return self.graph.degree[i] == 1

    def __init__(self, mol, atom, useFilters=True):

        self.mol = mol
        self.atom = atom
        self.excluded = set()

        if useFilters:
            # use a subset of filters from GDB17
            # Ruddigkeit et al., J. Chem. Inf. Model. 52, 2864 (2012)
            Chem.GetSSSR(self.mol)
            ringInfo = self.mol.GetRingInfo()
            self.excluded.update(BridgeheadAtoms(mol, ringInfo, 7))
            self.excluded.update(RingAtoms(ringInfo, 3))

        # subgraph of atoms that can make pi bonds
        self.graph = Graph(self.possiblePiBonds())

        # single coordinate - heteroatoms first
        done = False
        while not done:
            done = True
            for i in sorted(self.graph.nodes, key=lambda k: self.isCarbon(k)):
                if self.singleCoordinate(i):
                    for j in self.graph.neighbors(i):
                        if self.canBeTriple(i, j):
                            self.piBond(i, j, 3)
                        elif self.canBePi(i, j):
                            self.piBond(i, j, 2)
                        else:
                            continue
                        done = False
            if not done:
                self.reduceGraph()

        # isolated sp centers
        for i in self.graph.nodes:
            if self.canMultiplePiBond(i):
                for j in self.graph.neighbors(i):
                    if self.canMultiplePiBond(j):
                        break
                else:
                    for j in self.graph.neighbors(i):
                        if self.canBePi(i, j):
                            self.piBond(i, j, 2)

        # remaining pi bonds
        done = False
        while not done:
            done = True
            self.reduceGraph()
            for i, j in max_weight_matching(self.graph):
                self.incrementBond(i, j)
                done = False
