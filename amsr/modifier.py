from rdkit import Chem
from typing import List, Optional
from .encode import FromMolToTokens
from .decode import ToMol
from .markov import Markov
from .replacer import Replacer
from random import expovariate, shuffle


class Modifier:
    """modify a molecule by shuffling atom order, deleting tokens,
    then adding tokens

    :param mols: rdkit Mols, from which to draw token frequencies
    :param nDeleteAvg: average number of tokens to delete
    :param nAddMax: maximum number of tokens to add
    :param nReplaceAvg: number of token replacements
    """

    def __init__(
        self,
        mols: List[Chem.Mol],
        nDeleteAvg: Optional[int] = 2,
        nAddMax: Optional[int] = 10,
        nReplaceAvg: Optional[int] = 3,
    ):
        self.markov = Markov(mols)
        self.replacer = Replacer(mols)
        self.nDeleteAvg = nDeleteAvg
        self.nAddMax = nAddMax
        self.nReplaceAvg = nReplaceAvg

    def modify(self, mol: Chem.Mol) -> Chem.Mol:
        """modify given molecule

        :param mol: molecule to modify
        :return: modified molecule
        """
        n = mol.GetNumAtoms()
        i = list(range(n))
        shuffle(i)
        a = FromMolToTokens(Chem.RenumberAtoms(mol, i))
        if self.nDeleteAvg > 0:  # delete tokens
            k = min(round(expovariate(1 / self.nDeleteAvg)), len(a) - 1)
            if k > 0:
                a = a[:-k]
        if self.nReplaceAvg > 0:  # replace tokens
            for j in range(round(expovariate(1 / self.nReplaceAvg))):
                self.replacer.replaceToken(a)
        if self.nAddMax > 0:
            a += self.markov.generateTokens(nmax=self.nAddMax)
        return ToMol("".join(a))
