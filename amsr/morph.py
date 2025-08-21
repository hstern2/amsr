from Levenshtein import editops
from rdkit import Chem

from .decode import ToMol
from .encode import FromSmilesToTokens


class Morph:
    """morph between two molecules, by taking the minimum-edit pathway
    between their string representations

    :param s: list of AMSR tokens
    :param t: list of AMSR tokens
    """

    def __init__(self, s: list[str], t: list[str]):
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
        for amsr_str in self.amsr:
            m = ToMol(amsr_str)
            i = Chem.MolToInchi(m, options="-FixedH")
            if i0 is None or i != i0:
                self.mol.append(m)
                i0 = i

    @classmethod
    def fromSmiles(cls, s: str, t: str):
        """create morph from two SMILES strings

        :param s: SMILES
        :param t: SMILES
        :return: Morph object
        """
        return cls(FromSmilesToTokens(s), FromSmilesToTokens(t))

    def showAsSmiles(self):
        """display each mol in the morph as SMILES"""
        for m in self.mol:
            print(Chem.MolToSmiles(m))

    def asGridImage(self):
        from rdkit.Chem.Draw import MolsToGridImage

        """return mols as a grid image"""
        return MolsToGridImage(self.mol)
