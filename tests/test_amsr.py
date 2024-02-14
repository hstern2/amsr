import amsr, os, pandas
from random import seed
from rdkit import Chem
from itertools import combinations


def test_version():
    assert amsr.__version__


def _read_csv(csv_file):
    return pandas.read_csv(os.path.join(os.path.dirname(__file__), csv_file))


def _test_csv(csv_file):
    assert _read_csv(csv_file).apply(lambda m: amsr.CheckSmiles(m.SMILES), axis=1).all()


def test_Mg_compounds():
    _test_csv("Mg-compounds.csv")


def test_NP():
    _test_csv("natural_products.csv")


def test_FDA():
    _test_csv("some_FDA_approved_structures.csv")


def test_ertl():
    _test_csv("ertl_npsubstituents.csv")


def test_markov():
    seed(0)
    fda = _read_csv("some_FDA_approved_structures.csv")
    markov = amsr.Markov((Chem.MolFromSmiles(s) for s in fda["SMILES"]))
    for _ in range(100):
        assert amsr.CheckMol(amsr.ToMol(markov.generate()))


def test_modify():
    seed(0)
    fda = _read_csv("some_FDA_approved_structures.csv")
    np = _read_csv("natural_products.csv")
    modifier = amsr.Modifier((Chem.MolFromSmiles(s) for s in fda["SMILES"]))
    for _ in range(10):
        for m in (Chem.MolFromSmiles(s) for s in np["SMILES"]):
            assert amsr.CheckMol(modifier.modify(m))


def test_morph():
    seed(0)
    np = _read_csv("natural_products.csv")
    for s, t in combinations(np["SMILES"], 2):
        for m in amsr.morph.Morph.fromSmiles(s, t).mol:
            assert amsr.CheckMol(m)
