import amsr, os, pandas
from random import seed
from rdkit import Chem
from itertools import combinations

caffeine_smi = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
taxol_smi = (
    "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]"
    "([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)N"
    "C(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
)


def test_version():
    assert amsr.__version__


def test_caffeine():
    assert amsr.CheckSmiles(caffeine_smi)


def test_taxol():
    assert amsr.CheckSmiles(taxol_smi)


def test_cage():
    assert amsr.CheckAMSR("CCccCccc6oC..CCCC6C6.6")


def _read_csv(csv_file):
    return pandas.read_csv(
        os.path.join(os.path.dirname(__file__), "..", "data", csv_file)
    )


def _test_csv(csv_file, stringent=True):
    assert (
        _read_csv(csv_file)
        .apply(lambda m: amsr.CheckSmiles(m.SMILES, stringent), axis=1)
        .all()
    )


def test_Mg_compounds():
    _test_csv("Mg-compounds.csv", stringent=False)


def test_NP():
    _test_csv("natural_products.csv")


def test_FDA():
    _test_csv("some_FDA_approved_structures.csv")


def test_ertl():
    _test_csv("some_ertl_npsubs.csv")


def test_DEL():
    _test_csv("DEL_compounds.csv")


_model_path = os.path.join(os.path.dirname(__file__), "..", "models", "model.pth")


def test_lstm():
    lstm = amsr.LSTMModel.from_saved_model(_model_path)
    for _ in range(20):
        assert amsr.CheckAMSR(lstm.generate(["C"]))


def test_markov():
    seed(0)
    fda = _read_csv("some_FDA_approved_structures.csv")
    markov = amsr.Markov([Chem.MolFromSmiles(s) for s in fda["SMILES"]])
    for _ in range(20):
        assert amsr.CheckAMSR(markov.generate())


def test_modify():
    seed(0)
    fda = _read_csv("some_FDA_approved_structures.csv")
    np = _read_csv("natural_products.csv")
    modifier = amsr.Modifier(_model_path)
    for _ in range(10):
        for m in (Chem.MolFromSmiles(s) for s in np["SMILES"]):
            assert amsr.CheckMol(modifier.modify(m))


def test_morph():
    seed(0)
    np = _read_csv("natural_products.csv")
    for s, t in combinations(np["SMILES"][:10], 2):
        for a in amsr.morph.Morph.fromSmiles(s, t).amsr:
            assert amsr.CheckAMSR(a)


def test_no_stereo():
    assert amsr.CheckAMSR(amsr.FromSmiles(taxol_smi, useStereo=False))


def test_canonical():
    m = Chem.MolFromSmiles(caffeine_smi)
    s = amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True)
    for _ in range(20):
        assert amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True) == s
