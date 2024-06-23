import pandas
from amsr import FromMolToTokens
from lstm import LSTMModel
from rdkit.Chem import MolFromSmiles, RenumberAtoms
from random import shuffle

model = LSTMModel()
df = pandas.read_csv("chembl_33_filtered.csv")
smi = df.SMILES  # .sample(1000)
a = []
for s in smi:
    m = MolFromSmiles(s)
    k = list(range(m.GetNumAtoms()))
    for _ in range(10):
        shuffle(k)
        a.append(FromMolToTokens(RenumberAtoms(m, k)))
print(f"Training on {len(a)} AMSR strings")
pth = "model.pth"
model.train_and_save_model(a, pth)
print("Done.")
