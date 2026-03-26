import random

from rdkit import Chem

import amsr

from .conftest import read_csv

_N_RANDOM = 5


def _check_smiles_randomized(smiles, stringent=True):
    """Check canonical + N_RANDOM randomized round trips."""
    if not amsr.CheckSmiles(smiles, stringent):
        return False
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return True
    i1 = Chem.MolToInchi(m, options="-FixedH")
    for seed in range(_N_RANDOM):
        random.seed(seed)
        try:
            a = amsr.FromMol(m, stringent=stringent, randomize=True)
            m2 = amsr.ToMol(a, stringent=stringent)
        except Exception as e:
            print(f"Randomized round trip error (seed={seed}) for {smiles}: {e}")
            return False
        try:
            i2 = Chem.MolToInchi(m2, options="-FixedH")
        except Exception as e:
            print(f"Randomized round trip error (seed={seed}) for {smiles}: {e}")
            return False
        if i1 != i2:
            print(f"Randomized round trip failed (seed={seed}) for {smiles}")
            print(f"  AMSR: {a}")
            print(f"  original InChI: {i1}")
            print(f"  final    InChI: {i2}")
            return False
    return True


def _test_csv(csv_file: str, stringent: bool = True) -> None:
    df = read_csv(csv_file)
    results = df.apply(lambda m: _check_smiles_randomized(m.SMILES, stringent), axis=1)
    n_fail = (~results).sum()
    assert n_fail == 0, f"{n_fail}/{len(df)} molecules failed randomized round trip"


def test_Mg_compounds() -> None:
    _test_csv("Mg-compounds.csv", stringent=False)


def test_NP() -> None:
    _test_csv("natural_products.csv")


def test_FDA() -> None:
    _test_csv("some_FDA_approved_structures.csv", stringent=False)


def test_ertl() -> None:
    _test_csv("some_ertl_npsubs.csv")


def test_DEL() -> None:
    _test_csv("DEL_compounds.csv")


def test_chembl() -> None:
    _test_csv("chembl_35_5000.csv", stringent=False)
