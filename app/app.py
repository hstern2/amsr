import json
import math
import os
import sys

import rdkit.Chem.AllChem
import rdkit.Chem.Draw
from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import TPSA, MolWt
from rdkit.Chem.Lipinski import (
    HeavyAtomCount,
    NumHAcceptors,
    NumHDonors,
    NumRotatableBonds,
)
from rdkit.Chem.QED import qed

import amsr

# synthetic accessibility score; smaller means more accessible
# Ertl & Schuffenhauer, J. Cheminf. 2009, 1 (8)
sys.path.append(os.path.join(Chem.RDConfig.RDContribDir, "SA_Score"))
import sascorer

# store the last rendered molecule (with 2D coords) for alignment
prev_mol = None


def _align_2d(mol, ref):
    """Align mol's 2D coords to ref's using identity atom mapping.

    Uses least-squares rigid alignment on shared atoms (by index),
    and picks the better of no-flip vs X-flip to handle mirror images.
    """
    n = min(mol.GetNumAtoms(), ref.GetNumAtoms())
    if n < 2:
        return
    ref_conf = ref.GetConformer()
    conf = mol.GetConformer()
    px, py, qx, qy = [], [], [], []
    for i in range(n):
        p = ref_conf.GetAtomPosition(i)
        q = conf.GetAtomPosition(i)
        px.append(p.x)
        py.append(p.y)
        qx.append(q.x)
        qy.append(q.y)
    cpx = sum(px) / n
    cpy = sum(py) / n
    cqx = sum(qx) / n
    cqy = sum(qy) / n
    dpx = [x - cpx for x in px]
    dpy = [y - cpy for y in py]
    dqx = [x - cqx for x in qx]
    dqy = [y - cqy for y in qy]

    def _best_rotation(fqx, fqy):
        num = sum(fqx[i] * dpy[i] - fqy[i] * dpx[i] for i in range(n))
        den = sum(fqx[i] * dpx[i] + fqy[i] * dpy[i] for i in range(n))
        theta = math.atan2(num, den)
        c, s = math.cos(theta), math.sin(theta)
        ssd = sum(
            (c * fqx[i] - s * fqy[i] - dpx[i]) ** 2 + (s * fqx[i] + c * fqy[i] - dpy[i]) ** 2
            for i in range(n)
        )
        return c, s, ssd

    c0, s0, ssd0 = _best_rotation(dqx, dqy)
    c1, s1, ssd1 = _best_rotation([-x for x in dqx], dqy)
    use_flip = ssd1 < ssd0
    cos_r, sin_r = (c1, s1) if use_flip else (c0, s0)
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x, y = pos.x - cqx, pos.y - cqy
        if use_flip:
            x = -x
        conf.SetAtomPosition(i, (cos_r * x - sin_r * y + cpx, sin_r * x + cos_r * y + cpy, 0))


def flip_mol(m):
    "flip a molecule i.e. 180 degree rotation around vertical axis"
    conf = m.GetConformer()
    for a in m.GetAtoms():
        i = a.GetIdx()
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (-pos.x, pos.y, pos.z))
    Chem.AssignStereochemistry(m, force=True, cleanIt=True)


def rotate_mol(m, rotationValue: float):
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

    # Set bond notes for E/Z double bonds (smaller annotations)
    for bond in mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREOE:
            bond.SetProp("bondNote", "(E)")
        elif bond.GetStereo() == Chem.BondStereo.STEREOZ:
            bond.SetProp("bondNote", "(Z)")

    global prev_mol
    rdkit.Chem.AllChem.Compute2DCoords(mol)
    if prev_mol is not None:
        _align_2d(mol, prev_mol)
    prev_mol = Chem.RWMol(mol)
    if flip:
        flip_mol(mol)
    rotate_mol(mol, rotationValue)
    d = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(396, 396)
    opts = d.drawOptions()
    opts.annotationFontScale = 0.75  # Make annotations smaller
    opts.multipleBondOffset = 0.15  # Bring annotations closer
    opts.setBondNoteColour((0, 0, 0))  # Set bond note color to black
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
    try:
        rotationValue = float(request.form.get("rotationValue", "0"))
    except (ValueError, TypeError):
        rotationValue = 0.0
    svg, sdf, outString, ener, QED, tpsa, sa = [""] * 7
    hac, mw, clogp, hbd, hba, n_rot_bonds, passes_ro5 = [""] * 7
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
        # Additional properties from props.py
        hac_val = HeavyAtomCount(mol)
        mw_val = MolWt(mol)
        clogp_val = MolLogP(mol)
        hbd_val = NumHDonors(mol)
        hba_val = NumHAcceptors(mol)
        n_rot_bonds_val = NumRotatableBonds(mol)
        passes_ro5_val = mw_val <= 500 and clogp_val <= 5 and hbd_val <= 5 and hba_val <= 10

        hac = f"Heavy atom count: {hac_val}"
        mw = f"Molecular weight: {mw_val:.2f} Da"
        clogp = f"LogP: {clogp_val:.3f}"
        hbd = f"H-bond donors: {hbd_val}"
        hba = f"H-bond acceptors: {hba_val}"
        n_rot_bonds = f"Rotatable bonds: {n_rot_bonds_val}"
        passes_ro5 = f"Passes Rule of 5: {'Yes' if passes_ro5_val else 'No'}"

        if threeD:
            mol, ener = amsr.GetConformerAndEnergy(mol, dihedral=dih)
            sdf = Chem.MolToMolBlock(mol)
            ener = f"Energy: {ener:.3f} kcal/mol"
        if smiles_to_amsr:
            outString = amsr.FromMol(mol, stringent=stringent)
        else:
            outString = amsr.ToSmiles(inString, stringent=stringent)
    return json.dumps(
        {
            "svg": svg,
            "sdf": sdf,
            "outString": outString,
            "ener": ener,
            "qed": QED,
            "tpsa": tpsa,
            "sa": sa,
            "hac": hac,
            "mw": mw,
            "clogp": clogp,
            "hbd": hbd,
            "hba": hba,
            "n_rot_bonds": n_rot_bonds,
            "passes_ro5": passes_ro5,
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
