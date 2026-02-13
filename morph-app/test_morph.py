#!/usr/bin/env python3
"""Pytest suite for morph app."""

import sys
from pathlib import Path

import pytest

# Add project root so "import app" works
sys.path.insert(0, str(Path(__file__).parent.resolve()))


def test_default_constants():
    """Default SMILES are non-empty and valid-looking."""
    import app

    assert len(app.DEFAULT_SMILES_1) > 0
    assert len(app.DEFAULT_SMILES_2) > 0
    assert "C" in app.DEFAULT_SMILES_1
    assert "C" in app.DEFAULT_SMILES_2


def test_run_morph_requires_amsr():
    """run_morph raises ImportError when amsr is not installed."""
    import app

    try:
        import amsr
    except ImportError:
        amsr = None

    if amsr is None:
        with pytest.raises(ImportError):
            app.run_morph("C", "CC")
    else:
        pytest.skip("amsr is installed; run_morph is tested in test_run_morph_integration")


def test_run_morph_integration():
    """With amsr installed, run_morph returns (morph, smiles_text); morph.mol are RDKit Mols."""
    pytest.importorskip("amsr")

    import app

    smiles_1 = app.DEFAULT_SMILES_1
    smiles_2 = app.DEFAULT_SMILES_2

    morph, smiles_text = app.run_morph(smiles_1, smiles_2)

    assert isinstance(smiles_text, str)
    assert hasattr(morph, "mol")
    assert len(morph.mol) >= 2
