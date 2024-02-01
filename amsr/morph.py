from Levenshtein import editops
from .encode import FromSmilesToTokens
from .decode import ToMol
from typing import List
from rdkit import Chem


class Morph:
    """morph between two molecules, by taking the minimum-edit pathway
    between their string representations
    """

    def __init__(self, s: List[str], t: List[str]):
        """initialize from two lists of tokens"""
        s = list(s)[:]
        self.amsr = ["".join(s)]
        k = 0
        for op, i, j in editops(s, t):
            if op == "insert":
                s.insert(i + k, t[j])
                k += 1
            elif op == "delete":
                del s[i + k]
                k -= 1
            elif op == "replace":
                s[i + k] = t[j]
            self.amsr.append("".join(s))
        self.mol = []
        i0 = None
        for s in self.amsr:
            m = ToMol(s)
            i = Chem.MolToInchi(m, options="-FixedH")
            if i0 is None or i != i0:
                self.mol.append(m)
                i0 = i

    @classmethod
    def fromSmiles(cls, s, t):
        """initialize from two SMILES strings"""
        return cls(FromSmilesToTokens(s), FromSmilesToTokens(t))

    def showAsSmiles(self):
        """display each mol in the morph as SMILES"""
        for m in self.mol:
            print(Chem.MolToSmiles(m))
