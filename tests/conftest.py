import os

import pandas

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
MODEL_PATH = os.path.join(os.path.dirname(__file__), "..", "models", "model.pth")


def read_csv(csv_file: str):
    return pandas.read_csv(os.path.join(os.path.dirname(__file__), "..", "data", csv_file))
