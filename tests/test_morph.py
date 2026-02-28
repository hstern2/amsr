from itertools import combinations
from random import seed

import amsr

from .conftest import read_csv


def test_morph() -> None:
    seed(0)
    np = read_csv("natural_products.csv")
    for s, t in combinations(np["SMILES"][:10], 2):
        for a in amsr.morph.Morph.fromSmiles(s, t).amsr:
            assert amsr.CheckAMSR(a)
