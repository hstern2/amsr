from .encode import FromMolToTokens
from rdkit import Chem
from typing import Optional, List
from random import choices
from .count import Count

START_TOKEN = "[START]"
END_TOKEN = "[END]"


class Markov:
    """generate molecules using a simple Markov model

    :param mols: rdkit Mols, from which to draw token frequencies
    """

    def __init__(self, mols: List[Chem.Mol]):

        self.p = dict()
        self.p[START_TOKEN] = dict()
        for m in mols:
            t = FromMolToTokens(m, useGroups=True)
            if len(t) == 0:
                continue
            t0 = t[0]
            Count(self.p[START_TOKEN], t0)
            for i in range(1, len(t)):
                t0 = t[i - 1]
                t1 = t[i]
                if t0 not in self.p:
                    self.p[t0] = dict()
                Count(self.p[t0], t1)
            tn = t[-1]
            if tn not in self.p:
                self.p[tn] = dict()
            Count(self.p[tn], END_TOKEN)
        for k, v in self.p.items():
            self.p[k] = (list(v.keys()), list(v.values()))

    def _nextToken(self, t0):
        return choices(*self.p[t0])[0]

    def generateTokens(self, nmax: Optional[int] = -1):
        """generate sequence of tokens

        :param nmax: maximum number of tokens to generate
        :return: sequence of tokens
        """
        t = self._nextToken(START_TOKEN)
        n = 0
        while t != END_TOKEN and (nmax < 0 or n < nmax):
            yield t
            n += 1
            t = self._nextToken(t)

    def generate(self, nmax: Optional[int] = -1):
        """generate an AMSR string

        :param nmax: maximum length of string
        :return: AMSR string
        """
        return "".join(self.generateTokens(nmax=nmax))
