from .encode import FromMol, FromMolToTokens
from .decode import ToMol
from random import choices, gauss, shuffle
from .count import Count
import numpy
from rdkit import Chem


class Sampler:
    def __init__(self, mols):
        d = {}
        n = []
        for m in mols:
            tok = FromMolToTokens(m, useGroups=True)
            n.append(len(tok))
            for t in tok:
                Count(d, t)
        self.tokens = list(d.keys())
        self.weights = list(d.values())
        self.navg = numpy.mean(n)
        self.nstddev = numpy.std(n)

    def sampleTokens(self, navg=None, nstddev=None):
        if navg is None:
            navg = self.navg
        if nstddev is None:
            nstddev = self.nstddev
        return choices(
            population=self.tokens,
            weights=self.weights,
            k=max(round(gauss(navg, nstddev)), 1),
        )

    def sample(self, navg=None, nstddev=None):
        return "".join(self.sampleTokens(navg, nstddev))

    def decorate(self, mol, navg=10, nstddev=2):
        n = mol.GetNumAtoms()
        i = list(range(n))
        shuffle(i)
        return ToMol(FromMol(Chem.RenumberAtoms(mol, i)) + self.sample(navg, nstddev))
