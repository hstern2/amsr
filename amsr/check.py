from rdkit import Chem
from typing import Optional
from .encode import FromMol
from .decode import ToMol


def CheckMol(m1: Chem.Mol, stringent: Optional[bool] = False) -> bool:
    """Do round trip from RDKit mol m to AMSR.
    Compare InChI strings (with -FixedH) before and after round trip.

    :param m: RDKit Mol
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :return: do InChI strings match?
    """
    i1 = Chem.MolToInchi(m1, options="-FixedH")
    a = FromMol(m1, stringent=stringent)
    m2 = ToMol(a, stringent=stringent)
    i2 = Chem.MolToInchi(m2, options="-FixedH")
    if i1 == i2:
        return True
    else:
        print("CheckMol failed.")
        print(Chem.MolToSmiles(m1))
        print(i1)
        print(a)
        print(Chem.MolToSmiles(m2))
        print(i2)
        return False


def CheckSmiles(s: str, stringent: Optional[bool] = False) -> bool:
    """Convert SMILES s to RDKit Mol, then do round trip to AMSR.
    Compare InChI strings (with -FixedH) before and after round trip.

    :param s: SMILES
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :return: do InChI strings match?
    """
    return CheckMol(Chem.MolFromSmiles(s), stringent=stringent)
