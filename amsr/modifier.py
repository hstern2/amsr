from rdkit import Chem
from typing import List, Optional
from .encode import FromMolToTokens
from .decode import ToMol
from .lstm import LSTMModel
from random import expovariate, randrange


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
        model_path: str,
        nDeleteAvg: Optional[int] = 2,
        nAddMax: Optional[int] = 10,
        nReplaceAvg: Optional[int] = 3,
    ):
        self.model = LSTMModel.from_saved_model(model_path)
        self.nDeleteAvg = nDeleteAvg
        self.nAddMax = nAddMax
        self.nReplaceAvg = nReplaceAvg

    def modify(self, mol: Chem.Mol) -> Chem.Mol:
        """modify given molecule

        :param mol: molecule to modify
        :return: modified molecule
        """
        a = FromMolToTokens(mol, randomize=True)
        if self.nDeleteAvg > 0:  # delete tokens
            k = min(round(expovariate(1 / self.nDeleteAvg)), len(a) - 1)
            if k > 0:
                a = a[:-k]
        if self.nReplaceAvg > 0:  # replace tokens
            for _ in range(round(expovariate(1 / self.nReplaceAvg))):
                a = self.model.replace_token(a, randrange(len(a)))
        if self.nAddMax > 0:
            a = self.model.generate_tokens(a, self.nAddMax)

        return ToMol("".join(a))

    def modifySmiles(self, s: str) -> str:
        """modify given SMILES

        :param s: SMILES to modify
        :return: modified SMILES
        """
        return Chem.MolToSmiles(self.modify(Chem.MolFromSmiles(s)))
