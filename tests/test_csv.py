import amsr

from .conftest import read_csv


def _test_csv(csv_file: str, stringent: bool = True) -> None:
    assert read_csv(csv_file).apply(lambda m: amsr.CheckSmiles(m.SMILES, stringent), axis=1).all()


def test_Mg_compounds() -> None:
    _test_csv("Mg-compounds.csv", stringent=False)


def test_NP() -> None:
    _test_csv("natural_products.csv")


def test_FDA() -> None:
    _test_csv("some_FDA_approved_structures.csv")


def test_ertl() -> None:
    _test_csv("some_ertl_npsubs.csv")


def test_DEL() -> None:
    _test_csv("DEL_compounds.csv")
