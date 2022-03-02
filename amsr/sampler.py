from .encode import FromMolToTokens
from numpy import mean, std
from random import choices, gauss
from .count import Count


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
        self.navg = mean(n)
        self.nstddev = std(n)

    def sampleTokens(self):
        return choices(
            population=self.tokens,
            weights=self.weights,
            k=max(round(gauss(self.navg, self.nstddev)), 1),
        )

    def sample(self):
        return "".join(self.sampleTokens())
