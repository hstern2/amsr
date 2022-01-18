from rdkit import Chem
from .encode import FromMol
from .decode import ToMol

def CheckSmiles(s):
    return CheckMol(Chem.MolFromSmiles(s))

def CheckMol(m):
    i1 = Chem.MolToInchi(m, options='-FixedH')
    i2 = Chem.MolToInchi(ToMol(FromMol(m)), options='-FixedH')
    if i1 == i2:
        return True
    else:
        print(i1)
        print(i2)
        return False
