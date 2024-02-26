from rdkit import Chem
from typing import Optional
from .encode import FromMol
from .decode import ToMol


def CheckMol(m1: Chem.Mol, stringent: Optional[bool] = True) -> bool:
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
        print(f"CheckMol failed. stringent: {stringent}")
        print(f"AMSR: {a}")
        print(f"original InChI: {i1}")
        print(f"final    InChI: {i2}")
        print(f"original SMILES: {Chem.MolToSmiles(m1)}")
        print(f"final    SMILES: {Chem.MolToSmiles(m2)}")
        return False


def CheckSmiles(s: str, stringent: Optional[bool] = True) -> bool:
    """Convert SMILES s to RDKit Mol, then do round trip to AMSR.
    Compare InChI strings (with -FixedH) before and after round trip.

    :param s: SMILES
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :return: do InChI strings match?
    """
    m = Chem.MolFromSmiles(s)
    if m is None:
        print(f"rdkit returned None for {s}")
        return False
    if CheckMol(m, stringent=stringent):
        return True
    print(f"CheckSmiles (stringent: {stringent}) failed for {s}")
    return False


def CheckAMSR(s: str, stringent: Optional[bool] = True) -> bool:
    """Decode AMSR and check for valid molecule

    :param s: AMSR
    :param stringent: try to exclude unstable or synthetically inaccessible molecules
    :return: valid molecule?
    """
    m = ToMol(s, stringent=stringent)
    if CheckMol(m, stringent=stringent):
        return True
    print(f"CheckAMSR failed for {s}")
    return False
