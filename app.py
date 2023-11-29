from flask import Flask, render_template_string, request
from rdkit import Chem
import rdkit.Chem.Draw, rdkit.Chem.AllChem, json, amsr, math


template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>AMSR demo</title>
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"></script>
    <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css">
    <style>
        body {
            font-family: Arial;
        }
        #viewer2d {
            position: relative;
            height: 396px;
            width: 396px;
            border: solid black 1px;
            padding: 1px;
        }
        #viewer3d {
            position: relative;
            height: 400px;
            width: 400px;
        }
        #amsr, #smiles {
            font-family: Courier;
            font-size: 16px;
        }
    </style>
</head>
<body>
    <table>
    <tr><td>AMSR:

    </td>
    <td><input
            placeholder="enter AMSR"
            id="amsr"
            type="text"
            name="text"
            size="76"
            autocapitalize="off"
            autocorrect="off"
            autocomplete="on"
            onkeyup="amsr_changed()"
        ></input></td>
    </tr>
    <tr><td>SMILES:</td>
    <td><input
            placeholder="enter SMILES"
            id="smiles"
            type="text"
            name="text"
            size="76"
            autocapitalize="off"
            autocorrect="off"
            autocomplete="on"
            onkeyup="smiles_changed()"
        ></input>
    </td></tr>
    <tr>
    <td><input
            type="checkbox"
            id="flipX"
            name="flipX"
            onclick="amsr_changed()"
            >flipX
    <input
            type="checkbox"
            id="flipY"
            name="flipY"
            onclick="amsr_changed()"
            >flipY
    <button
            id="counterClockwise"
            name="counterClockwise"
            onclick="rotateCounterClockwise()"
            <i class="fas fa-rotate-left"></i>
    </button>
    <button
            id="clockwise"
            name="clockwise"
            onclick="rotateClockwise()"
            <i class="fas fa-rotate-right"></i>
    </button>
    <input
            type="checkbox"
            id="threeD"
            name="threeD"
            checked
            onclick="smiles_changed()"
            >3D
    </td></tr>
    </table>

    <table><tr>
        <td><div id="viewer2d"></div></td>
        <td><div id="viewer3d"></div></td>
    </tr></table>

    <script>

    function element(id) {
        return document.getElementById(id);
    }

    let viewer3d = $3Dmol.createViewer("viewer3d", {border: 'solid', backgroundColor: 'gainsboro'});
    let viewer2d = element("viewer2d");
    let pending = null;
    let rotationValue = 0;
    const rotationIncrement = 30; // degrees

    function refresh_viewer3d(sdf) {
        viewer3d.removeAllModels();
        viewer3d.addModel(sdf, "sdf")
        viewer3d.setStyle({}, {stick: {}});
        viewer3d.zoomTo();
        viewer3d.zoom(0.9);
        viewer3d.render();
    }

    function refresh_viewer2d(svg) {
        viewer2d.innerHTML = svg;
    }

    function rotateClockwise() {
        rotationValue -= rotationIncrement
        amsr_changed();
    }

    function rotateCounterClockwise() {
        rotationValue += rotationIncrement
        amsr_changed();
    }

    function smiles_changed() {
        if (element('threeD').checked) {
           element('viewer3d').style.display = 'block';
        } else {
           element('viewer3d').style.display = 'none';
        }
        if (pending)
            pending.abort();
        pending = $.post('/smiles_changed',
            {'smiles': element('smiles').value,
            'flipX': element('flipX').checked,
            'flipY': element('flipY').checked,
            'threeD': element('threeD').checked,
            'rotationValue': rotationValue},
            function(response) {
                var data = JSON.parse(response);
                refresh_viewer2d(data.svg);
                refresh_viewer3d(data.sdf);
                element('amsr').value = data.amsr;
            }
        );
    }

    function amsr_changed() {
        if (pending)
            pending.abort();
        pending = $.post('/amsr_changed',
            {'amsr': element('amsr').value,
            'flipX': element('flipX').checked,
            'flipY': element('flipY').checked,
            'threeD': element('threeD').checked,
            'rotationValue': rotationValue},
            function(response) {
                var data = JSON.parse(response);
                refresh_viewer2d(data.svg);
                refresh_viewer3d(data.sdf);
                element('smiles').value = data.smiles;
            }
        );
    }

    </script>

</body>
</html>
"""


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


@app.route("/", methods=["GET", "POST"])
def index():
    return render_template_string(template)


@app.route("/smiles_changed", methods=["GET", "POST"])
def smiles_changed():
    smiles = request.form.get("smiles")
    threeD = request.form.get("threeD") == "true"
    flipX = request.form.get("flipX") == "true"
    flipY = request.form.get("flipY") == "true"
    rotationValue = int(request.form.get("rotationValue"))  # degrees
    svg, sdf, a = "", "", ""
    mol = Chem.MolFromSmiles(smiles)
    if mol_isOK(mol):
        svg = get_svg(mol, flipX, flipY, rotationValue)
        if threeD:
            mol = amsr.GetConformer(mol)
            sdf = Chem.MolToMolBlock(mol)
        a = amsr.FromMol(mol)
    return json.dumps({"svg": svg, "sdf": sdf, "amsr": a})


@app.route("/amsr_changed", methods=["GET", "POST"])
def amsr_changed():
    a = request.form.get("amsr")
    threeD = request.form.get("threeD") == "true"
    flipX = request.form.get("flipX") == "true"
    flipY = request.form.get("flipY") == "true"
    rotationValue = int(request.form.get("rotationValue"))  # degrees
    smiles, svg, sdf = "", "", ""
    dihedral = dict()
    mol = amsr.ToMol(a, dihedral=dihedral)
    if mol_isOK(mol):
        smiles = Chem.MolToSmiles(mol)
        svg = get_svg(mol, flipX, flipY, rotationValue)
        if threeD:
            mol = amsr.GetConformer(mol, dihedral)
            sdf = Chem.MolToMolBlock(mol)
    return json.dumps({"svg": svg, "sdf": sdf, "smiles": smiles})


if __name__ == "__main__":
    app.run(debug=True)
