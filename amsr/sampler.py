from .encode import FromMolToTokens
from random import choices
from .count import Count


class Sampler:
    def __init__(self, mols):
        d = {}
        for m in mols:
            tok = FromMolToTokens(m, useGroups=True)
            for t in tok:
                Count(d, t)
        self.tokens = list(d.keys())
        self.weights = list(d.values())

    def sampleTokens(self, k):
        return choices(population=self.tokens, weights=self.weights, k=k)

    def sample(self, k):
        return "".join(self.sampleTokens(k))
