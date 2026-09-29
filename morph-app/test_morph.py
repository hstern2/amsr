#!/usr/bin/env python3
"""Pytest suite for morph app."""

import builtins
import subprocess
import sys
from pathlib import Path
from unittest import mock
from xml.etree import ElementTree

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
    """Legacy category labels still resolve to their catalog entries."""
    import morph_app

    resolved = morph_app.resolve_molecule_input("Epibatidine (Natural product)")

    assert resolved.name == "Epibatidine"
    assert resolved.source == "Natural product, PubChem CID 854023"


def test_autocomplete_options_are_name_only_and_broad():
    """Autocomplete is a large name-only list, not category-labeled options."""
    import morph_app

    options = morph_app.MOLECULE_AUTOCOMPLETE_OPTIONS

    assert len(options) >= 2000
    assert len({option.casefold() for option in options}) == len(options)
    assert all(option == option.lower() for option in options)
    assert all(not option[0].isdigit() for option in options)
    assert "epibatidine" in options
    assert "Ibogaine" not in options
    assert "ibogaine" in options
    assert "metformin" in options
    assert "abemaciclib" in options
    assert "1-octacosanol" not in options
    assert "Epibatidine (Natural product)" not in options


def test_autocomplete_data_file_is_plain_lowercase_text():
    """Vendored autocomplete data is fixed text, not executable Python."""
    import morph_app

    lines = morph_app.AUTOCOMPLETE_NAMES_PATH.read_text(encoding="ascii").splitlines()
    names = [line.strip() for line in lines if line.strip() and not line.startswith("#")]

    assert len(names) >= 2000
    assert len(set(names)) == len(names)
    assert all(name == name.lower() for name in names)
    assert all(not name[0].isdigit() for name in names)
    assert all("\t" not in name for name in names)


def test_autocomplete_matches_prioritizes_prefixes_without_requiring_selection():
    """Autocomplete suggestions do not replace free-form molecule lookup."""
    import morph_app

    matches = morph_app.autocomplete_matches("met")

    assert "metformin" in matches
    assert len(matches) <= 6
    assert morph_app.autocomplete_matches("metformin") == ()
    assert morph_app.autocomplete_matches("C1=CC=CC=C1") == ()


def test_morph_input_snapshot_normalizes_committed_text_values():
    """Committed input changes are compared independent of surrounding whitespace."""
    import morph_app

    assert morph_app.morph_input_snapshot(" metformin ", " CCO\n", True) == (
        "metformin",
        "CCO",
        True,
    )


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


def test_resolve_molecule_input_removes_counterions_from_raw_smiles(monkeypatch):
    """Only the largest connected endpoint component is retained."""
    pytest.importorskip("rdkit.Chem")

    import morph_app

    def fail_lookup(name):
        raise AssertionError(f"Unexpected lookup for {name}")

    monkeypatch.setattr(morph_app, "lookup_pubchem_molecule", fail_lookup)
    resolved = morph_app.resolve_molecule_input("[Na+].CC(=O)[O-]")

    assert resolved.smiles == "CC(=O)[O-]"
    assert "." not in resolved.smiles


def test_resolve_molecule_input_removes_counterions_from_name_lookup(monkeypatch):
    """PubChem-resolved salts are reduced before they become endpoints."""
    pytest.importorskip("rdkit.Chem")

    import morph_app

    monkeypatch.setattr(
        morph_app,
        "lookup_pubchem_molecule",
        lambda name: morph_app.MoleculeResolution(
            "Trimethylammonium chloride", "C[NH+](C)C.[Cl-]", "PubChem"
        ),
    )

    resolved = morph_app.resolve_molecule_input("trimethylammonium chloride")

    assert resolved.smiles == "C[NH+](C)C"


def test_smiles_classifier_leaves_hyphenated_names_for_lookup():
    """Hyphenated molecule names should not be mistaken for raw SMILES."""
    import morph_app

    assert morph_app.looks_like_smiles("CCO")
    assert morph_app.looks_like_smiles("C-C")
    assert not morph_app.looks_like_smiles("L-DOPA")


def test_lilly_filter_default_is_enabled():
    """The Streamlit option defaults to applying Lilly Medchem Rules."""
    import morph_app

    assert morph_app.DEFAULT_APPLY_LILLY_FILTER is True


def test_filter_mols_with_lilly_uses_relaxed_and_filters(monkeypatch):
    """The Lilly wrapper passes -relaxed and keeps only accepted morph outputs."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mols = [Chem.MolFromSmiles(smiles) for smiles in ("CCCCCCCC", "CCO", "c1ccccc1O")]
    captured = {}

    monkeypatch.setattr(morph_app.shutil, "which", lambda name: f"/usr/local/bin/{name}")

    def fake_run(cmd, cwd, text, capture_output, check, timeout):
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        captured["input"] = Path(cmd[-1]).read_text(encoding="utf-8")
        captured["text"] = text
        captured["capture_output"] = capture_output
        captured["check"] = check
        captured["timeout"] = timeout
        return subprocess.CompletedProcess(
            cmd,
            0,
            stdout="CCCCCCCC morph_0000\nOc1ccccc1 morph_0002\n",
            stderr="",
        )

    monkeypatch.setattr(morph_app.subprocess, "run", fake_run)

    result = morph_app.filter_mols_with_lilly(mols)

    assert captured["cmd"][:2] == ["/usr/local/bin/Lilly_Medchem_Rules.rb", "-relaxed"]
    assert captured["cwd"] in captured["cmd"][-1]
    assert captured["text"] is True
    assert captured["capture_output"] is True
    assert captured["check"] is False
    assert captured["timeout"] == morph_app.LILLY_FILTER_TIMEOUT_SECONDS
    assert "morph_0000" in captured["input"]
    assert "morph_0001" in captured["input"]
    assert "morph_0002" in captured["input"]
    assert len(result.mols) == 2
    assert result.rejected_count == 1
    assert result.smiles_text.splitlines() == ["CCCCCCCC", "Oc1ccccc1"]


def test_pathway_to_smi_text_names_endpoints_and_intermediates():
    """The downloadable SMI includes endpoint and intermediate record names."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mols = [Chem.MolFromSmiles(smiles) for smiles in ("CCO", "CCCC", "c1ccccc1O")]
    molecule_1 = morph_app.MoleculeResolution("Ethanol", "CCO", "Manual input")
    molecule_2 = morph_app.MoleculeResolution("Phenol", "c1ccccc1O", "PubChem CID 996")

    assert morph_app.pathway_to_smi_text(mols, molecule_1, molecule_2).splitlines() == [
        "CCO\tendpoint_1: Ethanol (Manual input)",
        "CCCC\tintermediate_001",
        "Oc1ccccc1\tendpoint_2: Phenol (PubChem CID 996)",
    ]


def test_pathway_to_cdxml_is_labeled_editable_and_round_trips():
    """The ChemDraw download contains separately editable, labeled molecules."""
    Chem = pytest.importorskip("rdkit.Chem")
    if not hasattr(Chem, "MolToCDXMLBlock") or not Chem.HasChemDrawCDXSupport():
        pytest.skip("RDKit was built without ChemDraw CDXML support")

    import morph_app

    smiles = ["CCO", "CCCC", "c1ccccc1O"]
    mols = [Chem.MolFromSmiles(value) for value in smiles]
    molecule_1 = morph_app.MoleculeResolution("Ethanol", smiles[0], "Manual input")
    molecule_2 = morph_app.MoleculeResolution("Phenol", smiles[-1], "PubChem CID 996")

    cdxml = morph_app.pathway_to_cdxml(mols, molecule_1, molecule_2, columns=2)
    root = ElementTree.fromstring(cdxml)

    assert root.tag == "CDXML"
    assert root.get("CreationProgram") == "AMSR Morph"
    assert len(root.findall("./page/fragment")) == 3
    assert [text.find("s").text for text in root.findall("./page/t")] == [
        "1. From: Ethanol",
        "2. Intermediate 001",
        "3. To: Phenol",
    ]
    ids = [element.get("id") for element in root.iter() if element.get("id")]
    assert len(ids) == len(set(ids))

    round_tripped = Chem.MolsFromCDXML(cdxml)
    assert len(round_tripped) == len(mols)
    assert [Chem.MolToSmiles(mol) for mol in round_tripped] == [
        Chem.MolToSmiles(mol) for mol in mols
    ]


def test_pathway_to_cdxml_places_structures_in_requested_grid():
    """The ChemDraw pathway is laid out in rows instead of overlapping structures."""
    Chem = pytest.importorskip("rdkit.Chem")
    if not hasattr(Chem, "MolToCDXMLBlock") or not Chem.HasChemDrawCDXSupport():
        pytest.skip("RDKit was built without ChemDraw CDXML support")

    import morph_app

    mols = [Chem.MolFromSmiles(value) for value in ("CC", "CCC", "CCCC")]
    root = ElementTree.fromstring(morph_app.pathway_to_cdxml(mols, columns=2))
    fragment_centers = []
    for fragment in root.findall("./page/fragment"):
        min_x, min_y, max_x, max_y = morph_app._cdxml_fragment_bounds(fragment)
        fragment_centers.append(((min_x + max_x) / 2, (min_y + max_y) / 2))

    assert fragment_centers[0][1] == pytest.approx(fragment_centers[1][1])
    assert fragment_centers[0][0] < fragment_centers[1][0]
    assert fragment_centers[2][1] > fragment_centers[0][1]


def test_pathway_to_cdxml_preserves_detailed_chemical_information():
    """CDXML export retains stereochemistry, charges, isotopes, and radicals."""
    Chem = pytest.importorskip("rdkit.Chem")
    if not hasattr(Chem, "MolToCDXMLBlock") or not Chem.HasChemDrawCDXSupport():
        pytest.skip("RDKit was built without ChemDraw CDXML support")

    import morph_app

    smiles = [
        "F[C@](Cl)(Br)I",
        "C/C=C\\C",
        "[13CH3][NH3+]",
        "[O-][N+](=O)c1ccccc1",
        "[2H]O[2H]",
        "[CH3]",
    ]
    mols = [Chem.MolFromSmiles(value) for value in smiles]
    round_tripped = Chem.MolsFromCDXML(morph_app.pathway_to_cdxml(mols))

    assert [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in round_tripped] == [
        Chem.MolToSmiles(mol, isomericSmiles=True) for mol in mols
    ]


def test_streamlit_app_offers_cdxml_download():
    """A completed morph exposes both pathway download formats in the UI."""
    pytest.importorskip("rdkit.Chem")
    pytest.importorskip("amsr")
    pytest.importorskip("streamlit.testing.v1")
    from streamlit.testing.v1 import AppTest

    app_path = Path(__file__).with_name("morph_app.py")
    app = AppTest.from_file(str(app_path), default_timeout=20)
    app.run()
    app.text_input(key="molecule_1").set_value("CCO")
    app.text_input(key="molecule_2").set_value("CCN")
    app.checkbox[0].uncheck()
    app.button[0].click()
    app.run()

    assert not app.exception
    assert [button.label for button in app.get("download_button")] == [
        "Download .smi file",
        "Download CDXML for ChemDraw",
    ]


def test_filter_morph_output_with_lilly_preserves_input_endpoints(monkeypatch):
    """The Lilly filter applies only to generated intermediates, not endpoints."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mols = [Chem.MolFromSmiles(smiles) for smiles in ("CCO", "CCCCCCCC", "c1ccccc1O")]
    seen = {}

    def fake_filter(intermediates):
        seen["smiles"] = [Chem.MolToSmiles(mol) for mol in intermediates]
        return morph_app.LillyFilterResult([], "", 1)

    monkeypatch.setattr(morph_app, "filter_mols_with_lilly", fake_filter)

    result = morph_app.filter_morph_output_with_lilly(mols)

    assert seen["smiles"] == ["CCCCCCCC"]
    assert result.rejected_count == 1
    assert result.smiles_text.splitlines() == ["CCO", "Oc1ccccc1"]


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


def test_mols_to_svgs_renders_kekulized_aromatic_molecules():
    """Aromatic systems are drawn as alternating bonds, not dashed aromatic bonds."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mol = Chem.MolFromSmiles("c1ccccc1")
    svgs = morph_app.mols_to_svgs([mol], mol_size=120)

    assert len(svgs) == 1
    assert "<svg" in svgs[0]
    assert "stroke-dasharray" not in svgs[0]
    assert any(atom.GetIsAromatic() for atom in mol.GetAtoms())


def test_mols_to_svgs_omits_molecules_that_cannot_be_kekulized(monkeypatch):
    """Molecule rendering drops structures RDKit cannot kekulize."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    def fail_kekulize(_mol, clearAromaticFlags):
        raise ValueError("cannot kekulize")

    monkeypatch.setattr(Chem, "Kekulize", fail_kekulize)

    mol = Chem.MolFromSmiles("c1ccccc1")

    assert morph_app.kekulized_mol_for_drawing(mol) is None
    assert morph_app.mols_to_svgs([mol], mol_size=120) == []


def test_molecule_properties_match_amsr_2d_labels():
    """Morph hover properties use the same descriptor labels as the AMSR 2D view."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mol = Chem.MolFromSmiles("CCO")
    properties = morph_app.molecule_properties(mol)

    assert [label for label, _value in properties] == [
        "QED score",
        "TPSA",
        "SA score",
        "Heavy atom count",
        "Molecular weight",
        "LogP",
        "H-bond donors",
        "H-bond acceptors",
        "Rotatable bonds",
        "Passes Rule of 5",
    ]
    assert dict(properties)["TPSA"].endswith("&#8491;<sup>2</sup>")
    assert dict(properties)["Passes Rule of 5"] == "Yes"


def test_molecule_properties_use_traditional_lipinski_hbond_counts():
    """HBA/HBD match RDKit's traditional Lipinski definitions, not newer perception."""
    Chem = pytest.importorskip("rdkit.Chem")
    import morph_app
    from rdkit.Chem.Lipinski import NHOHCount, NOCount, NumHAcceptors, NumHDonors

    for smiles in ("CC(=O)N", "CC(=O)Nc1ccc(O)cc1"):
        mol = Chem.MolFromSmiles(smiles)
        properties = dict(morph_app.molecule_properties(mol))

        assert properties["H-bond donors"] == str(NHOHCount(mol))
        assert properties["H-bond acceptors"] == str(NOCount(mol))

    regression_mol = Chem.MolFromSmiles("CC(=O)N")
    assert NHOHCount(regression_mol) != NumHDonors(regression_mol)
    assert NOCount(regression_mol) != NumHAcceptors(regression_mol)


def test_molecule_hover_grid_html_embeds_property_tooltips():
    """Rendered morph molecules expose descriptor tooltips on hover/focus."""
    Chem = pytest.importorskip("rdkit.Chem")

    import morph_app

    mol = Chem.MolFromSmiles("CCO")
    html = morph_app.molecule_hover_grid_html([mol], mol_size=120)

    assert '<div class="mol-card"' in html
    assert 'class="mol-tooltip"' in html
    assert 'role="tooltip"' in html
    assert "<strong>QED score:</strong>" in html
    assert "<strong>Passes Rule of 5:</strong> Yes" in html
    assert "width: 120px" in html
    assert "width: 240px" in html
    assert "position: fixed" in html
    assert "placeTooltip" in html
