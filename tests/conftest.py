import os

import pandas

_ROOT = os.path.join(os.path.dirname(__file__), "..")
SDF_DIR = os.path.join(_ROOT, "data", "sdf")
MODEL_PATH = os.path.join(_ROOT, "models", "model.pth")


def read_csv(csv_file: str):
    return pandas.read_csv(os.path.join(_ROOT, "data", "csv", csv_file))
