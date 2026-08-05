#!/usr/bin/env python3
"""Pytest suite for morph app."""

import builtins
import sys
from pathlib import Path
from unittest import mock

import pytest

# Add morph-app dir so "import morph_app" works
sys.path.insert(0, str(Path(__file__).parent.resolve()))


def test_default_constants():
    """Default SMILES are non-empty and valid-looking."""
    import morph_app

    assert len(morph_app.DEFAULT_SMILES_1) > 0
    assert len(morph_app.DEFAULT_SMILES_2) > 0
    assert "C" in morph_app.DEFAULT_SMILES_1
    assert "C" in morph_app.DEFAULT_SMILES_2


def test_lookup_known_molecule_supports_names_and_aliases():
    """Curated molecule names and aliases resolve locally."""
    import morph_app

    aspirin = morph_app.lookup_known_molecule("aspirin")
    paracetamol = morph_app.lookup_known_molecule("Paracetamol")
    levodopa = morph_app.lookup_known_molecule("L-DOPA")

    assert aspirin.name == "Aspirin"
    assert aspirin.smiles == "CC(=O)OC1=CC=CC=C1C(=O)O"
    assert paracetamol.name == "Acetaminophen"
    assert levodopa.name == "Levodopa"


def test_resolve_molecule_input_accepts_catalog_selection():
    """A catalog key resolves directly from the single molecule field."""
    import morph_app

    resolved = morph_app.resolve_molecule_input("caffeine")

    assert resolved.name == "Caffeine"
    assert resolved.source == "Natural product, PubChem CID 2519"


def test_resolve_molecule_input_accepts_display_label():
    """Visible selectbox labels resolve to their catalog entries."""
    import morph_app

    resolved = morph_app.resolve_molecule_input("Epibatidine (Natural product)")

    assert resolved.name == "Epibatidine"
    assert resolved.source == "Natural product, PubChem CID 854023"


def test_catalog_has_requested_curated_sets():
    """Catalog includes traceable natural products and FDA-approved drugs."""
    import morph_app

    collections = {molecule.collection for molecule in morph_app.KNOWN_MOLECULES.values()}

    assert len(morph_app.KNOWN_MOLECULES) >= 100
    assert collections == {"Natural product", "FDA-approved drug"}
    assert "benzene" not in morph_app.KNOWN_MOLECULES
    assert "ethanol" not in morph_app.KNOWN_MOLECULES


def test_epibatidine_is_pubchem_verified():
    """Epibatidine stays pinned to the PubChem-resolved SMILES and CID."""
    import morph_app

    epibatidine = morph_app.KNOWN_MOLECULES["epibatidine"]
    resolved = morph_app.lookup_known_molecule("Epibatidine")

    assert epibatidine.pubchem_cid == 854023
    assert epibatidine.smiles == "C1C[C@@H]2[C@H](C[C@H]1N2)C3=CN=C(C=C3)Cl"
    assert morph_app.DEFAULT_SMILES_2 == epibatidine.smiles
    assert resolved.name == "Epibatidine"
    assert resolved.source == "Natural product, PubChem CID 854023"


def test_resolve_molecule_input_accepts_raw_smiles(monkeypatch):
    """Obvious SMILES strings do not go through the network lookup."""
    import morph_app

    def fail_lookup(name):
        raise AssertionError(f"Unexpected lookup for {name}")

    monkeypatch.setattr(morph_app, "lookup_pubchem_molecule", fail_lookup)
    resolved = morph_app.resolve_molecule_input("CCO")

    assert resolved.name == "Custom SMILES"
    assert resolved.smiles == "CCO"
    assert resolved.source == "Manual input"


def test_smiles_classifier_leaves_hyphenated_names_for_lookup():
    """Hyphenated molecule names should not be mistaken for raw SMILES."""
    import morph_app

    assert morph_app.looks_like_smiles("CCO")
    assert morph_app.looks_like_smiles("C-C")
    assert not morph_app.looks_like_smiles("L-DOPA")


def test_resolve_molecule_input_uses_pubchem_fallback(monkeypatch):
    """Unknown non-SMILES names can be resolved by the PubChem fallback."""
    import morph_app

    monkeypatch.setattr(
        morph_app,
        "lookup_pubchem_molecule",
        lambda name: morph_app.MoleculeResolution(
            "Adenosine", "C1=NC=NC2=C1N=CN2", "PubChem CID 60961"
        ),
    )

    resolved = morph_app.resolve_molecule_input("adenosine")

    assert resolved.name == "Adenosine"
    assert resolved.source == "PubChem CID 60961"


def test_resolve_molecule_input_raises_for_unknown_name(monkeypatch):
    """A failed name lookup produces a clear validation error before morphing."""
    import morph_app

    monkeypatch.setattr(morph_app, "lookup_pubchem_molecule", lambda name: None)

    with pytest.raises(ValueError, match="Could not look up 'not a molecule'"):
        morph_app.resolve_molecule_input("not a molecule")


def test_run_morph_requires_amsr():
    """run_morph raises ImportError when amsr is not available."""
    import morph_app

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "amsr":
            raise ImportError("No module named 'amsr'")
        return real_import(name, *args, **kwargs)

    with mock.patch.object(builtins, "__import__", side_effect=fake_import):
        with pytest.raises(ImportError):
            morph_app.run_morph("C", "CC")


def test_run_morph_integration():
    """With amsr installed, run_morph returns (morph, smiles_text); morph.mol are RDKit Mols."""
    pytest.importorskip("amsr")

    import morph_app

    smiles_1 = morph_app.DEFAULT_SMILES_1
    smiles_2 = morph_app.DEFAULT_SMILES_2

    morph, smiles_text = morph_app.run_morph(smiles_1, smiles_2)

    assert isinstance(smiles_text, str)
    assert hasattr(morph, "mol")
    assert len(morph.mol) >= 2


def test_mols_to_svgs_renders_rdkit_molecules():
    """Molecule rendering uses RDKit's width/height MolToSVG signature."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mol = Chem.MolFromSmiles("CCO")
    svgs = morph_app.mols_to_svgs([mol], mol_size=120)

    assert len(svgs) == 1
    assert "<svg" in svgs[0]
