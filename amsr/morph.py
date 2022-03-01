from Levenshtein import editops
from .encode import FromSmilesToTokens
from .decode import ToMol
from rdkit import Chem


class Morph:
    def __init__(self, s, t):
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
            m = ToMol(s, contiguous=True, useFilters=True)
            i = Chem.MolToInchi(m, options="-FixedH")
            if i0 is None or i != i0:
                self.mol.append(m)
                i0 = i

    @classmethod
    def fromSmiles(cls, s, t):
        return cls(FromSmilesToTokens(s), FromSmilesToTokens(t))

    def showAsSmiles(self):
        for m in self.mol:
            print(Chem.MolToSmiles(m))
