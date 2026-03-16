import glob
import os

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolAlign

import amsr

_data_dir = os.path.join(os.path.dirname(__file__), "data")
_out_dir = os.path.join(os.path.dirname(__file__), "out")

_sdf_files = sorted(glob.glob(os.path.join(_data_dir, "*.sdf")))


@pytest.mark.parametrize("sdf_path", _sdf_files, ids=[os.path.basename(f) for f in _sdf_files])
def test_roundtrip(sdf_path):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    assert mol is not None, f"Could not parse {sdf_path}"

    dihedral = {}
    s = amsr.FromMol(mol)
    mol2 = amsr.ToMol(s, dihedral=dihedral)
    mol3 = amsr.GetConformer(mol2, dihedral=dihedral)

    rmsd = rdMolAlign.GetBestRMS(mol3, mol)

    os.makedirs(_out_dir, exist_ok=True)
    Chem.MolToMolFile(mol3, os.path.join(_out_dir, f"{name}_zmatrix.sdf"))

    assert mol3.GetConformer().Is3D()
    assert rmsd < 1.0, f"RMSD {rmsd:.3f} Å too large for {name}"
