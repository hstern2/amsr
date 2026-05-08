import csv
import glob
import os

import pytest
from filelock import FileLock

from amsr.roundtrip import N_RANDOM_SEEDS, SDF_RMSD_THRESHOLD, RoundtripSDF

from .conftest import SDF_DIR

_out_dir = os.path.join(os.path.dirname(__file__), "out")
_csv_path = os.path.join(_out_dir, "out.csv")
_lock_path = _csv_path + ".lock"


def pytest_generate_tests(metafunc):
    if "sdf_path" in metafunc.fixturenames:
        sdf_files = sorted(glob.glob(os.path.join(SDF_DIR, "*.sdf")))
        metafunc.parametrize("sdf_path", sdf_files, ids=[os.path.basename(f) for f in sdf_files])
    if "seed" in metafunc.fixturenames:
        seeds = [None] + list(range(N_RANDOM_SEEDS))
        metafunc.parametrize("seed", seeds, ids=[f"seed{s}" for s in range(len(seeds))])


@pytest.fixture(scope="session", autouse=True)
def _csv_output(worker_id, tmp_path_factory):
    if worker_id == "master":
        # Not running under xdist — no coordination needed
        os.makedirs(_out_dir, exist_ok=True)
        with open(_csv_path, "w", newline="") as f:
            csv.writer(f).writerow(["name", "seed", "amsr", "rmsd", "time_s"])
    else:
        # xdist workers: use a shared lock so only the first worker writes the header
        root_tmp = tmp_path_factory.getbasetemp().parent
        lock = FileLock(str(root_tmp / "csv_header.lock"))
        os.makedirs(_out_dir, exist_ok=True)
        with lock:
            if not os.path.exists(_csv_path) or os.path.getsize(_csv_path) == 0:
                with open(_csv_path, "w", newline="") as f:
                    csv.writer(f).writerow(["name", "seed", "amsr", "rmsd", "time_s"])
    yield
    # Sort only once: in master mode always, in xdist only on gw0
    if worker_id in ("master", "gw0"):
        with FileLock(_lock_path):
            with open(_csv_path, newline="") as f:
                rows = list(csv.reader(f))
            if rows:
                header, data = rows[0], rows[1:]
                data.sort(key=lambda r: -float(r[3]) if len(r) > 3 and r[3] else 0)
                with open(_csv_path, "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(header)
                    w.writerows(data)


def test_roundtrip_sdf(sdf_path, seed):
    r = RoundtripSDF(sdf_path, seed, output_dir=_out_dir)

    assert r["status"] != "ERROR", f"Error: {r.get('error', 'could not parse')}"

    with FileLock(_lock_path):
        with open(_csv_path, "a", newline="") as f:
            csv.writer(f).writerow(
                [r["name"], r["seed"], r["amsr"], f"{r['rmsd']:.3f}", f"{r['time']:.3f}"]
            )

    assert r["rmsd"] < SDF_RMSD_THRESHOLD, (
        f"RMSD {r['rmsd']:.3f} Å too large for {r['name']}" f" (seed={r['seed']})"
    )
