import amsr, os
import pandas as pd


def test_version():
    assert amsr.__version__


def test_methane():
    assert amsr.CheckSmiles("C")


def test_caffeine():
    assert amsr.CheckSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")


def test_taxol():
    assert amsr.CheckSmiles(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )


def _test_csv(csv_file):
    assert (
        pd.read_csv(os.path.join(os.path.dirname(__file__), csv_file))
        .apply(lambda m: amsr.CheckSmiles(m.SMILES), axis=1)
        .all()
    )


def test_FDA():
    _test_csv("some_FDA_approved_structures.csv")


def test_NP():
    _test_csv("natural_products.csv")
