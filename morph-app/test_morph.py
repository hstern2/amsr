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
