import json
import math
import os
import sys

import rdkit.Chem.AllChem
import rdkit.Chem.Draw
from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem.Descriptors import TPSA
from rdkit.Chem.QED import qed

import amsr

# synthetic accessibility score; smaller means more accessible
# Ertl & Schuffenhauer, J. Cheminf. 2009, 1 (8)
sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
import sascorer


def flip_mol(m):
    "flip a molecule i.e. 180 degree rotation around vertical axis"
    conf = m.GetConformer()
    for a in m.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (-pos.x, pos.y, pos.z))
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)


def rotate_mol(m, rotationValue: int):
    "rotate molecule in XY plane"
    conf = m.GetConformer()
    radians = math.radians(rotationValue)
    for a in m.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        newX = pos.x * math.cos(radians) - pos.y * math.sin(radians)
        newY = pos.x * math.sin(radians) + pos.y * math.cos(radians)
        conf.SetAtomPosition(i, (newX, newY, pos.z))


def get_svg(mol, flip: bool, rotationValue: int):
    # Minimal stereochemistry assignment without adding atoms
    Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)

    # Set atom notes for stereocenters (smaller annotations)
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            cip_code = atom.GetProp("_CIPCode")
            atom.SetProp("atomNote", f"({cip_code})")

    rdkit.Chem.AllChem.Compute2DCoords(mol)
    if flip:
        flip_mol(mol)
    rotate_mol(mol, rotationValue)
    d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(396, 396)
    opts = d.drawOptions()
    opts.annotationFontScale = 0.75  # Make annotations smaller
    opts.multipleBondOffset = 0.15  # Bring annotations closer
    actives = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp("_active")]
    d.DrawMolecule(mol, highlightAtoms=actives)
    d.FinishDrawing()
    return d.GetDrawingText().replace("svg:", "")


def mol_isOK(mol):
    return mol and mol.GetNumAtoms() > 0


app = Flask(__name__)
methods = ["GET", "POST"]
model_path = os.path.join(os.path.dirname("__file__"), "..", "models", "model.pth")
lstm = amsr.LSTMModel.from_saved_model(model_path)


@app.route("/", methods=methods)
def index():
    return render_template("index.html")


@app.route("/generate", methods=methods)
def generate():
    seed = request.form.get("seed")
    if len(seed) == 0:
        seed = "C"
    return json.dumps({"amsr": lstm.generate(amsr.ToTokens(seed))})


@app.route("/mol_changed", methods=methods)
def mol_changed():
    inString = request.form.get("inString")
    smiles_to_amsr = request.form.get("smiles_to_amsr") == "true"
    threeD = request.form.get("threeD") == "true"
    stringent = request.form.get("stringent") == "true"
    flipMol = request.form.get("flipMol") == "true"
    rotationValue = int(request.form.get("rotationValue"))  # degrees
    svg, sdf, outString, ener, QED, tpsa, sa = [""] * 7
    dih = {}
    mol = (
        Chem.MolFromSmiles(inString)
        if smiles_to_amsr
        else amsr.ToMol(inString, stringent=stringent, dihedral=dih)
    )
    if mol_isOK(mol):
        svg = get_svg(mol, flipMol, rotationValue)
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
            outString = amsr.ToSmiles(inString, stringent=stringent, isomericSmiles=True)
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
