from rdkit import Chem

import amsr

caffeine_smi = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
taxol_smi = (
    "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]"
    "([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)N"
    "C(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
)


def test_version() -> None:
    assert amsr.__version__


def test_caffeine() -> None:
    assert amsr.CheckSmiles(caffeine_smi)


def test_taxol() -> None:
    assert amsr.CheckSmiles(taxol_smi)


def test_cage() -> None:
    assert amsr.CheckAMSR("CCccCccc6oC..CCCC6C6.6")


def test_no_stereo() -> None:
    assert amsr.CheckAMSR(amsr.FromSmiles(taxol_smi, useStereo=False))


def test_canonical() -> None:
    m = Chem.MolFromSmiles(caffeine_smi)
    s = amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True)
    for _ in range(20):
        assert amsr.FromSmiles(Chem.MolToSmiles(m, doRandom=True), canonical=True) == s
