import amsr, os
import pandas as pd


def test_version():
    assert amsr.__version__


def _test_csv(csv_file):
    assert (
        pd.read_csv(os.path.join(os.path.dirname(__file__), csv_file))
        .apply(lambda m: amsr.CheckSmiles(m.SMILES), axis=1)
        .all()
    )


def test_Mg():
    _test_csv("Mg.csv")


def test_NP():
    _test_csv("natural_products.csv")


def test_FDA():
    _test_csv("some_FDA_approved_structures.csv")
