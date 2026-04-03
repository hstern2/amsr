import csv
import glob
import os
import time

import pytest
from rdkit import Chem

import amsr.cost_grad as _cg
from amsr.roundtrip import N_RANDOM_SEEDS, Roundtrip

from .conftest import SDF_DIR

_out_dir = os.path.join(os.path.dirname(__file__), "out")
_csv_path = os.path.join(_out_dir, "out.csv")


def pytest_generate_tests(metafunc):
    if "sdf_path" in metafunc.fixturenames:
        sdf_files = sorted(glob.glob(os.path.join(SDF_DIR, "*.sdf")))
        metafunc.parametrize("sdf_path", sdf_files, ids=[os.path.basename(f) for f in sdf_files])
    if "seed" in metafunc.fixturenames:
        seeds = [None] + list(range(N_RANDOM_SEEDS))
        metafunc.parametrize("seed", seeds, ids=[f"seed{s}" for s in range(len(seeds))])


@pytest.fixture(scope="session", autouse=True)
def _csv_output():
    os.makedirs(_out_dir, exist_ok=True)
    with open(_csv_path, "w", newline="") as f:
        csv.writer(f).writerow(["name", "backend", "seed", "amsr", "rmsd", "time_s"])
    yield
    with open(_csv_path, newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    data.sort(key=lambda r: -float(r[4]) if r[4] else 0)
    with open(_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)


def test_roundtrip_sdf(sdf_path, seed):
    name = os.path.splitext(os.path.basename(sdf_path))[0]
    mol = Chem.MolFromMolFile(sdf_path, removeHs=True)
    assert mol is not None, f"Could not parse {sdf_path}"

    backend = "C" if _cg.is_available() else "Python"
    seed_str = "0" if seed is None else str(seed + 1)

    t0 = time.time()
    s, rmsd, mol_out = Roundtrip(mol, seed=seed)
    elapsed = time.time() - t0

    os.makedirs(_out_dir, exist_ok=True)
    match = mol.GetSubstructMatch(mol_out)
    if match:
        mol_reordered = Chem.RenumberAtoms(mol, list(match))
    else:
        mol_reordered = mol
    Chem.MolToMolFile(mol_reordered, os.path.join(_out_dir, f"{name}_{seed_str}_orig.sdf"))
    Chem.MolToMolFile(mol_out, os.path.join(_out_dir, f"{name}_{seed_str}_out.sdf"))

    with open(_csv_path, "a", newline="") as f:
        csv.writer(f).writerow([name, backend, seed_str, s, f"{rmsd:.3f}", f"{elapsed:.3f}"])

    assert mol_out.GetConformer().Is3D()
    assert (
        rmsd < 1.0
    ), f"RMSD {rmsd:.3f} Å too large for {name} (seed={seed_str}, backend={backend})"
