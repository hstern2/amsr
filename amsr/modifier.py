from rdkit import Chem
from typing import Iterable, Optional
from .encode import FromMolToTokens
from .decode import ToMol
from .markov import Markov
from random import expovariate, shuffle


class Modifier:
    """modify a molecule by shuffling atom order, deleting tokens,
    then adding tokens

    :param mols: rdkit Mols, from which to draw token frequencies
    :param nDeleteAvg: average number of tokens to delete
    :param nAddMax: maximum number of tokens ato add
    """

    def __init__(
        self,
        mols: Iterable[Chem.Mol],
        nDeleteAvg: Optional[int] = 2,
        nAddMax: Optional[int] = 10,
    ):
        self.markov = Markov(mols)
        self.nDeleteAvg = nDeleteAvg
        self.nAddMax = nAddMax

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
        if self.nAddMax > 0:
            a += self.markov.generateTokens(nmax=self.nAddMax)
        return ToMol("".join(a))
