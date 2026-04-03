import os

import pandas
import pytest

_ROOT = os.path.join(os.path.dirname(__file__), "..")
SDF_DIR = os.path.join(_ROOT, "data", "sdf")
MODEL_PATH = os.path.join(_ROOT, "models", "model.pth")


def pytest_addoption(parser):
    parser.addoption(
        "--sdf-dir",
        default=None,
        help="Override SDF directory for roundtrip tests (default: data/sdf/)",
    )


@pytest.fixture(scope="session")
def sdf_dir(request):
    custom = request.config.getoption("--sdf-dir")
    return custom if custom else SDF_DIR


def read_csv(csv_file: str):
    return pandas.read_csv(os.path.join(_ROOT, "data", "csv", csv_file))
