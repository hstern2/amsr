import csv
import glob
import os
import random
import time

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolAlign

import amsr.cost_grad as _cg
from amsr.decode import ToMol
from amsr.encode import FromMol
from amsr.zmatrix import GetConformer

from .conftest import SDF_DIR

_out_dir = os.path.join(os.path.dirname(__file__), "out")
_csv_path = os.path.join(_out_dir, "out.csv")

_sdf_files = sorted(glob.glob(os.path.join(SDF_DIR, "*.sdf")))
_N_RANDOM = 5  # number of randomized encodings per molecule


@pytest.fixture(scope="session", autouse=True)
def _csv_output():
    os.makedirs(_out_dir, exist_ok=True)
    with open(_csv_path, "w", newline="") as f:
        csv.writer(f).writerow(["name", "backend", "seed", "amsr", "rmsd", "time_s"])
    yield
    with open(_csv_path, newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    data.sort(key=lambda r: (-float(r[4]) if r[4] else 0))
    with open(_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)


def _roundtrip(mol, seed=None):
    """Roundtrip with optional randomized encoding.  Returns (amsr, rmsd, mol_out, elapsed)."""
    if seed is not None:
        random.seed(seed)
    t0 = time.time()
    dihedral = {}
    s = FromMol(mol, randomize=(seed is not None))
    mol2 = ToMol(s, dihedral=dihedral)
    mol3 = GetConformer(mol2, dihedral=dihedral)
    elapsed = time.time() - t0
    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    return s, rmsd, mol3, elapsed


def _backend_name():
    return "C" if _cg.is_available() else "Python"


# Seeds: None = default (canonical), then 5 deterministic random seeds
_seeds = [None] + list(range(1000, 1000 + _N_RANDOM))


@pytest.mark.parametrize("sdf_path", _sdf_files, ids=[os.path.basename(f) for f in _sdf_files])
@pytest.mark.parametrize(
    "seed", _seeds, ids=["canonical"] + [f"seed{s}" for s in range(1000, 1000 + _N_RANDOM)]
)
def test_roundtrip_sdf(sdf_path, seed):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    assert mol is not None, f"Could not parse {sdf_path}"

    s, rmsd, mol_out, elapsed = _roundtrip(mol, seed=seed)
    backend = _backend_name()

    os.makedirs(_out_dir, exist_ok=True)
    seed_str = "canonical" if seed is None else str(seed)
    if seed is None:
        # Save SDF output for the canonical encoding only
        match = mol.GetSubstructMatch(mol_out)
        if match:
            mol_reordered = Chem.RenumberAtoms(mol, list(match))
        else:
            mol_reordered = mol
        Chem.MolToMolFile(mol_reordered, os.path.join(_out_dir, f"{name}_original.sdf"))
        Chem.MolToMolFile(mol_out, os.path.join(_out_dir, f"{name}_out.sdf"))

    with open(_csv_path, "a", newline="") as f:
        csv.writer(f).writerow([name, backend, seed_str, s, f"{rmsd:.3f}", f"{elapsed:.3f}"])

    assert mol_out.GetConformer().Is3D()
    assert (
        rmsd < 0.8
    ), f"RMSD {rmsd:.3f} Å too large for {name} (seed={seed_str}, backend={backend})"
