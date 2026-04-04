import csv
import glob
import os

import pytest

import amsr.cost_grad as _cg
from amsr.roundtrip import N_RANDOM_SEEDS, RoundtripSDF

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
        csv.writer(f).writerow(["name", "seed", "amsr", "rmsd", "time_s"])
    yield
    with open(_csv_path, newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    data.sort(key=lambda r: -float(r[3]) if r[3] else 0)
    with open(_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(data)


def test_roundtrip_sdf(sdf_path, seed):
    r = RoundtripSDF(sdf_path, seed, output_dir=_out_dir)

    assert r["status"] != "ERROR", f"Error: {r.get('error', 'could not parse')}"

    backend = "C" if _cg.is_available() else "Python"
    with open(_csv_path, "a", newline="") as f:
        csv.writer(f).writerow(
            [r["name"], r["seed"], r["amsr"], f"{r['rmsd']:.3f}", f"{r['time']:.3f}"]
        )

    assert r["rmsd"] < 1.0, (
        f"RMSD {r['rmsd']:.3f} Å too large for {r['name']}"
        f" (seed={r['seed']}, backend={backend})"
    )
