import amsr, os, pandas
from random import seed
from rdkit import Chem


def test_version():
    assert amsr.__version__


def _read_csv(csv_file):
    return pandas.read_csv(os.path.join(os.path.dirname(__file__), csv_file))


def _test_csv(csv_file):
    assert _read_csv(csv_file).apply(lambda m: amsr.CheckSmiles(m.SMILES), axis=1).all()


def test_elements():
    _test_csv("elements.csv")


def test_Mg_compounds():
    _test_csv("Mg-compounds.csv")


def test_NP():
    _test_csv("natural_products.csv")


def test_FDA():
    _test_csv("some_FDA_approved_structures.csv")


# def test_chembl_33():
#    df = pandas.read_csv(os.path.join(os.path.dirname(__file__), "chembl_33.smi"), delimiter=" ")
#    assert(df.apply(lambda m: amsr.CheckSmiles(m.canonical_smiles), axis=1).all())


def test_sampler():
    seed(0)
    df = _read_csv("some_FDA_approved_structures.csv")
    sampler = amsr.Sampler((Chem.MolFromSmiles(s) for s in df["SMILES"]))
    for _ in range(1000):
        assert amsr.CheckMol(amsr.ToMol(sampler.sample(20)))
