from rdkit import Chem
from typing import List, Optional
from .encode import FromMolToTokens
from .decode import ToMol
from .sampler import Sampler
from random import expovariate, shuffle


class Modifier:
    def __init__(
        self,
        mols: List[Chem.Mol],
        nDeleteAvg: Optional[int] = 2,
        nAddAvg: Optional[int] = 5,
    ):
        self.sampler = Sampler(mols)
        self.nDeleteAvg = nDeleteAvg
        self.nAddAvg = nAddAvg

    def modify(self, mol: Chem.Mol):
        n = mol.GetNumAtoms()
        i = list(range(n))
        shuffle(i)
        a = FromMolToTokens(Chem.RenumberAtoms(mol, i))
        if self.nDeleteAvg > 0:  # delete tokens
            k = min(round(expovariate(1 / self.nDeleteAvg)), len(a) - 1)
            if k > 0:
                a = a[:-k]
        if self.nAddAvg > 0:
            k = max(round(expovariate(1 / self.nAddAvg)), 1)
            a += self.sampler.sampleTokens(k)
        return ToMol("".join(a))
