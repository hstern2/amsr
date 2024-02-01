from flask import Flask, render_template, request
from rdkit import Chem
import rdkit.Chem.Draw, rdkit.Chem.AllChem
import pandas, json, amsr, math, sys, os
from rdkit.Chem.QED import qed
from rdkit.Chem.Descriptors import TPSA
from random import expovariate

# synthetic accessibility score; smaller means more accessible
# Ertl & Schuffenhauer, J. Cheminf. 2009, 1 (8)
sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
import sascorer


def flip_mol(m, axis):
    conf = m.GetConformer()
    for a in m.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        if axis == "X":
            conf.SetAtomPosition(i, (-pos.x, pos.y, pos.z))
        elif axis == "Y":
            conf.SetAtomPosition(i, (pos.x, -pos.y, pos.z))
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)


def rotate_mol(m, rotationValue: int):
    # rotationValue in degrees
    conf = m.GetConformer()
    radians = math.radians(rotationValue)
    for a in m.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        newX = pos.x * math.cos(radians) - pos.y * math.sin(radians)
        newY = pos.x * math.sin(radians) + pos.y * math.cos(radians)
        conf.SetAtomPosition(i, (newX, newY, pos.z))


def get_svg(mol, flipX: bool, flipY: bool, rotationValue: int):
    rdkit.Chem.AllChem.Compute2DCoords(mol)
    if flipX:
        flip_mol(mol, "X")
    if flipY:
        flip_mol(mol, "Y")
    rotate_mol(mol, rotationValue)
    d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(396, 396)
    actives = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp("_active")]
    d.DrawMolecule(mol, highlightAtoms=actives)
    d.FinishDrawing()
    return d.GetDrawingText()


def mol_isOK(mol):
    return mol and mol.GetNumAtoms() > 0


app = Flask(__name__)
methods = ["GET", "POST"]
url = "https://raw.githubusercontent.com/hstern2/amsr/main/tests/some_FDA_approved_structures.csv"
df = pandas.read_csv(url)
markov = amsr.Markov((Chem.MolFromSmiles(s) for s in df["SMILES"]))


@app.route("/", methods=methods)
def index():
    return render_template("index.html")


@app.route("/random_mol", methods=methods)
def random_mol():
    k = max(round(expovariate(1 / 20)), 1)
    return json.dumps({"amsr": markov.sample(nmax=25)})


@app.route("/mol_changed", methods=methods)
def mol_changed():
    inString = request.form.get("inString")
    smiles_to_amsr = request.form.get("smiles_to_amsr") == "true"
    threeD = request.form.get("threeD") == "true"
    stringent = request.form.get("stringent") == "true"
    flipX = request.form.get("flipX") == "true"
    flipY = request.form.get("flipY") == "true"
    rotationValue = int(request.form.get("rotationValue"))  # degrees
    svg, sdf, outString, ener, QED, tpsa, sa = [""] * 7
    dih = {}
    mol = (
        Chem.MolFromSmiles(inString)
        if smiles_to_amsr
        else amsr.ToMol(inString, dihedral=dih)
    )
    if mol_isOK(mol):
        svg = get_svg(mol, flipX, flipY, rotationValue)
        QED = f"QED score: {qed(mol):.3f}"
        tpsa = f"TPSA: {TPSA(mol):.3f} &#8491;<sup>2</sup>"
        sa = f"SA score: {sascorer.calculateScore(mol):.3f}"
        if threeD:
            mol, ener = amsr.GetConformerAndEnergy(mol, dihedral=dih)
            sdf = Chem.MolToMolBlock(mol)
            ener = f"Energy: {ener:.3f} kcal/mol"
        if smiles_to_amsr:
            outString = amsr.FromMol(mol, stringent=stringent)
        else:
            outString = Chem.MolToSmiles(mol)
    return json.dumps(
        {
            "svg": svg,
            "sdf": sdf,
            "outString": outString,
            "ener": ener,
            "qed": QED,
            "tpsa": tpsa,
            "sa": sa,
        }
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    p = ArgumentParser(description="flask app to demo AMSR")
    default_port = 5000
    p.add_argument(
        "-p",
        "--port",
        type=int,
        help=f"server port (default={default_port})",
        default=default_port,
    )
    a = p.parse_args()
    app.run(debug=True, port=a.port)
