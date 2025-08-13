import re
from pathlib import Path
from sys import stdout

import pandas as pd
from typer import Argument, Option, Typer

from amsr import FromSmilesToTokens, LSTMModel, ToSmiles

DEFAULT_T = 1.0
DEFAULT_MAX_LENGTH = 20
DEFAULT_N = 10
DEFAULT_MODEL = Path(__file__).resolve().parent / "models" / "model.pth"
DEFAULT_ATTACHMENT_POINT = "[1C]"

app = Typer(add_completion=False, context_settings={"help_option_names": ["-h", "--help"]})


def _grow_from_seed(
    smi: str,
    name: str,
    attachment_point: str,
    temperature: float,
    max_length: int,
    model: LSTMModel,
    n: int,
) -> pd.DataFrame:
    a = FromSmilesToTokens(smi)
    i_attach = a.index(attachment_point)
    stem = a[:i_attach] + [re.sub(r"\[[0-9]+", "", attachment_point)]
    g = []
    for _ in range(n):
        b = model.generate_tokens(stem, max_length=max_length, temperature=temperature)
        g.append(ToSmiles("".join(stem + ["&"] + b[i_attach + 1 :] + ["&"] + a[i_attach + 1 :])))
    return pd.DataFrame({"SMILES": g, "NAME": [f"{name}_gen{_}" for _ in range(n)]})


def _grow(
    seeds: pd.DataFrame,
    attachment_point: str,
    temperature: float,
    max_length: int,
    num_seqs: int,
    model_path: str,
):
    # Load the model
    model = LSTMModel.from_saved_model(model_path)
    df = pd.DataFrame()
    for _, row in seeds.iterrows():
        df = pd.concat(
            [
                df,
                _grow_from_seed(
                    row.SMILES,
                    row.NAME,
                    attachment_point,
                    temperature,
                    max_length,
                    model,
                    num_seqs,
                ),
            ]
        )
    return df


@app.command()
def grow(
    seeds: str = Argument(help=".smi file containing seeds"),
    attachment_point=Option(
        DEFAULT_ATTACHMENT_POINT, "--attachment_point", "-a", help="attachment point"
    ),
    temperature: float = Option(DEFAULT_T, "--temperature", "-t", help="temperature"),
    max_length: int = Option(DEFAULT_MAX_LENGTH, "--max_length", "-l", help="maximum length"),
    num_seqs: int = Option(
        DEFAULT_N, "--num_seqs", "-n", help="number of sequences to generate per seed"
    ),
    model_path: str = Option(DEFAULT_MODEL, "--model_path", "-m", help="path to model file"),
):
    """Generate AMSR sequences from a seed SMILES string with an attachment point."""
    s = pd.read_csv(seeds, sep=r"\s+", header=None, names=["SMILES", "NAME"])
    _grow(s, attachment_point, temperature, max_length, num_seqs, model_path).to_csv(
        stdout, sep=" ", index=False, header=False
    )


if __name__ == "__main__":
    app()
