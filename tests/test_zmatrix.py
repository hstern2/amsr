import csv
import glob
import os

import pytest
from rdkit import Chem

import amsr

_data_dir = os.path.join(os.path.dirname(__file__), "data")
_out_dir = os.path.join(os.path.dirname(__file__), "out")
_csv_path = os.path.join(_out_dir, "test_out.csv")

_sdf_files = sorted(glob.glob(os.path.join(_data_dir, "*.sdf")))


@pytest.fixture(scope="session", autouse=True)
def _csv_header():
    os.makedirs(_out_dir, exist_ok=True)
    with open(_csv_path, "w", newline="") as f:
        csv.writer(f).writerow(["name", "amsr", "rmsd"])
    yield


@pytest.mark.parametrize("sdf_path", _sdf_files, ids=[os.path.basename(f) for f in _sdf_files])
def test_roundtrip(sdf_path):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    assert mol is not None, f"Could not parse {sdf_path}"

    s, rmsd, mol3 = amsr.Roundtrip(mol)

    os.makedirs(_out_dir, exist_ok=True)
    Chem.MolToMolFile(mol, os.path.join(_out_dir, f"{name}_original.sdf"))
    Chem.MolToMolFile(mol3, os.path.join(_out_dir, f"{name}_zmatrix.sdf"))

    with open(_csv_path, "a", newline="") as f:
        csv.writer(f).writerow([name, s, f"{rmsd:.3f}"])

    assert mol3.GetConformer().Is3D()
    assert rmsd < 1.0, f"RMSD {rmsd:.3f} Å too large for {name}"
