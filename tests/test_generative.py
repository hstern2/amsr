from random import seed

from rdkit import Chem

import amsr

from .conftest import MODEL_PATH, read_csv


def test_lstm() -> None:
    lstm = amsr.LSTMModel.from_saved_model(MODEL_PATH)
    for _ in range(20):
        assert amsr.CheckAMSR(lstm.generate(["C"]))


def test_markov() -> None:
    seed(0)
    fda = read_csv("some_FDA_approved_structures.csv")
    markov = amsr.Markov([Chem.MolFromSmiles(s) for s in fda["SMILES"]])
    for _ in range(20):
        assert amsr.CheckAMSR(markov.generate())


def test_modify() -> None:
    seed(0)
    np = read_csv("natural_products.csv")
    modifier = amsr.Modifier(MODEL_PATH)
    for _ in range(10):
        for m in (Chem.MolFromSmiles(s) for s in np["SMILES"]):
            assert amsr.CheckMol(modifier.modify(m))
