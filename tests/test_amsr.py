import amsr, os
import pandas as pd


def test_version():
    assert amsr.__version__


def test_methane():
    assert amsr.CheckSmiles("C")


def test_taxol():
    assert amsr.CheckSmiles(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )


def test_some_FDA_Approved_structures():
    pd.read_csv(
        os.path.join(os.path.dirname(__file__), "some_FDA_Approved_structures.csv")
    ).apply(lambda m: amsr.CheckSmiles(m.SMILES), axis=1)
