from rdkit import Chem
from .encode import FromMol
from .decode import ToMol


def CheckMol(m: Chem.Mol) -> bool:
    """Do round trip from RDKit mol m to AMSR.
    Compare InChI strings (with -FixedH) before and after round trip.

    :param m: RDKit Mol
    :return: do InChI strings match?
    """
    i1 = Chem.MolToInchi(m, options="-FixedH")
    i2 = Chem.MolToInchi(ToMol(FromMol(m)), options="-FixedH")
    if i1 == i2:
        return True
    else:
        print(i1)
        print(i2)
        return False


def CheckSmiles(s: str) -> bool:
    """Convert SMILES s to RDKit Mol, then do round trip to AMSR.
    Compare InChI strings (with -FixedH) before and after round trip.

    :param s: SMILES
    :return: do InChI strings match?
    """
    return CheckMol(Chem.MolFromSmiles(s))
