from amsr import ToTokens, ToSmiles, LSTMModel
from argparse import ArgumentParser
import os

DEFAULT_T = 1.0
DEFAULT_MAX_LENGTH = 20
DEFAULT_N = 10
DEFAULT_MODEL = os.path.join(os.path.dirname(__file__), "models", "model.pth")

ap = ArgumentParser(description="generate AMSR sequences")
ap.add_argument("s", help="initial AMSR sequence")
ap.add_argument(
    "-t",
    "--temperature",
    type=float,
    help=f"temperature. Default: {DEFAULT_T}",
    default=DEFAULT_T,
)
ap.add_argument(
    "-l",
    "--max_length",
    type=int,
    help=f"maximum length. Default: {DEFAULT_MAX_LENGTH}",
    default=DEFAULT_MAX_LENGTH,
)
ap.add_argument(
    "-n",
    "--num_seqs",
    type=int,
    help=f"number of sequences. Default: {DEFAULT_N}",
    default=DEFAULT_N,
)
ap.add_argument(
    "-m",
    "--model",
    help=f"path to model file. Default: {DEFAULT_MODEL}",
    default=DEFAULT_MODEL,
)
a = ap.parse_args()

loaded_model = LSTMModel.from_saved_model(a.model)
toks = ToTokens(a.s)
print(f"Initial tokens: {toks}")
print("Generated sequences:")
for _ in range(a.num_seqs):
    print(
        ToSmiles(
            loaded_model.generate(
                toks, temperature=a.temperature, max_length=a.max_length
            )
        )
    )
