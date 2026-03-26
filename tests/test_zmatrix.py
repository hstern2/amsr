import csv
import glob
import os

import pytest
from rdkit import Chem

import amsr

from .conftest import SDF_DIR

_out_dir = os.path.join(os.path.dirname(__file__), "out")
_csv_path = os.path.join(_out_dir, "out.csv")

_sdf_files = sorted(glob.glob(os.path.join(SDF_DIR, "*.sdf")))


@pytest.fixture(scope="session", autouse=True)
def _csv_output():
    os.makedirs(_out_dir, exist_ok=True)
    with open(_csv_path, "w", newline="") as f:
        csv.writer(f).writerow(["name", "amsr", "rmsd"])
    yield
    with open(_csv_path, newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    data.sort(key=lambda r: -float(r[2]))
    with open(_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)


@pytest.mark.parametrize("sdf_path", _sdf_files, ids=[os.path.basename(f) for f in _sdf_files])
def test_roundtrip_sdf(sdf_path):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    assert mol is not None, f"Could not parse {sdf_path}"

    s, rmsd, mol_out = amsr.Roundtrip(mol)

    os.makedirs(_out_dir, exist_ok=True)
    match = mol.GetSubstructMatch(mol_out)
    if match:
        mol_reordered = Chem.RenumberAtoms(mol, list(match))
    else:
        mol_reordered = mol
    Chem.MolToMolFile(mol_reordered, os.path.join(_out_dir, f"{name}_original.sdf"))
    Chem.MolToMolFile(mol_out, os.path.join(_out_dir, f"{name}_out.sdf"))

    with open(_csv_path, "a", newline="") as f:
        csv.writer(f).writerow([name, s, f"{rmsd:.3f}"])

    assert mol_out.GetConformer().Is3D()
    assert rmsd < 0.8, f"RMSD {rmsd:.3f} Å too large for {name}"
