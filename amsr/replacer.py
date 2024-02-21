from .count import Count
from .encode import FromMolToTokens
from typing import List
from rdkit import Chem
from random import randrange, choices


class Replacer:
    """replace tokens in an AMSR string

    :param mols: rdkit Mols, from which to draw token frequencies
    """

    def __init__(self, mols: List[Chem.Mol]):
        p = dict()
        for m in mols:
            for t in FromMolToTokens(m):
                Count(p, t)
        self.tokens = list(p.keys())
        self.counts = list(p.values())

    def replaceToken(self, t: List[str]):
        if t:
            t[randrange(len(t))] = choices(self.tokens, self.counts)[0]
