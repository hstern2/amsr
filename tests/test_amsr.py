import amsr
from rdkit.Chem import MolFromSmiles


def test_C():
    assert amsr.CheckMol(MolFromSmiles("C"))


def test_sane():
    assert amsr.__version__
