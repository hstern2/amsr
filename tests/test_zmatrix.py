import os

from rdkit import Chem
from rdkit.Chem import rdMolAlign

import amsr

_data_dir = os.path.join(os.path.dirname(__file__), "data")
_out_dir = os.path.join(os.path.dirname(__file__), "data", "output")


def _load_sdf(name):
    return Chem.MolFromMolFile(os.path.join(_data_dir, name), removeHs=True)


def _roundtrip(mol, name):
    """Encode 3D mol to AMSR, decode, generate z-matrix conformer.
    Align to original, write both SDFs, return RMSD."""
    dihedral = {}
    s = amsr.FromMol(mol)
    print(f"  AMSR: {s}")
    mol2 = amsr.ToMol(s, dihedral=dihedral)
    mol3 = amsr.GetConformer(mol2, dihedral=dihedral)

    # Compute best RMSD (handles molecular symmetry)
    rmsd = rdMolAlign.GetBestRMS(mol3, mol)
    # Align for SDF output
    match = mol.GetSubstructMatch(mol3)
    if match:
        atom_map = [(i, match[i]) for i in range(mol3.GetNumAtoms())]
        rdMolAlign.AlignMol(mol3, mol, atomMap=atom_map)

    os.makedirs(_out_dir, exist_ok=True)
    Chem.MolToMolFile(mol, os.path.join(_out_dir, f"{name}_original.sdf"))
    Chem.MolToMolFile(mol3, os.path.join(_out_dir, f"{name}_zmatrix.sdf"))

    print(f"  RMSD: {rmsd:.3f} Å")
    assert mol3.GetConformer().Is3D()
    return rmsd


# --- Simple chains ---


def test_ethane():
    _roundtrip(_load_sdf("ethane_3d.sdf"), "ethane")


def test_propane():
    _roundtrip(_load_sdf("propane_3d.sdf"), "propane")


def test_butane():
    _roundtrip(_load_sdf("butane_3d.sdf"), "butane")


def test_neopentane():
    _roundtrip(_load_sdf("neopentane_3d.sdf"), "neopentane")


# --- Rings ---


def test_benzene():
    _roundtrip(_load_sdf("benzene_3d.sdf"), "benzene")


def test_cyclohexane():
    _roundtrip(_load_sdf("cyclohexane_3d.sdf"), "cyclohexane")


# --- Mixed ---


def test_aspirin():
    _roundtrip(_load_sdf("aspirin_3d.sdf"), "aspirin")


def test_ibuprofen():
    _roundtrip(_load_sdf("ibuprofen_3d.sdf"), "ibuprofen")


def test_caffeine():
    _roundtrip(_load_sdf("caffeine_3d.sdf"), "caffeine")


# --- Amino acids ---

_AMINO_ACIDS = [
    "glycine",
    "alanine",
    "valine",
    "leucine",
    "isoleucine",
    "proline",
    "phenylalanine",
    "tryptophan",
    "serine",
    "threonine",
    "cysteine",
    "methionine",
    "aspartate",
    "glutamate",
    "asparagine",
    "glutamine",
    "lysine",
    "arginine",
    "histidine",
    "tyrosine",
]


def test_amino_acids():
    for name in _AMINO_ACIDS:
        print(f"{name}:")
        _roundtrip(_load_sdf(f"{name}_3d.sdf"), f"aa_{name}")
