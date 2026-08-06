#!/usr/bin/env python3
"""Morph: Streamlit app for molecular morph (two SMILES -> pathway). Based on morph.ipynb."""

import json
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from molecule_catalog import (
    KNOWN_MOLECULE_KEYS,
    KNOWN_MOLECULES,
)


@dataclass(frozen=True)
class MoleculeResolution:
    name: str
    smiles: str
    source: str


@dataclass(frozen=True)
class LillyFilterResult:
    mols: list
    smiles_text: str
    rejected_count: int


DEFAULT_MOLECULE_1_KEY = "ibogaine"
DEFAULT_MOLECULE_2_KEY = "epibatidine"
DEFAULT_SMILES_1 = KNOWN_MOLECULES[DEFAULT_MOLECULE_1_KEY].smiles
DEFAULT_SMILES_2 = KNOWN_MOLECULES[DEFAULT_MOLECULE_2_KEY].smiles
DEFAULT_APPLY_LILLY_FILTER = True
LILLY_MEDCHEM_RULES = "Lilly_Medchem_Rules.rb"
LILLY_FILTER_TIMEOUT_SECONDS = 120


def _molecule_lookup_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


MOLECULE_NAME_INDEX = {
    _molecule_lookup_key(name): key
    for key, molecule in KNOWN_MOLECULES.items()
    for name in (molecule.name, key, *molecule.aliases)
}
SMILES_CHARS = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)\\/#=%.:]+$")
SMILES_ATOM_TOKEN = re.compile(
    r"Br|Cl|Si|Se|Na|Li|Mg|Ca|Zn|Fe|Al|Ag|As|Au|Ba|Bi|Cd|Co|Cu|Hg|Mn|Ni|Pb|Pt|Sn|Ti|"
    r"B|C|N|O|P|S|F|I|H|K|V|Y|W|U|b|c|n|o|p|s"
)


def format_molecule_option(key: str) -> str:
    molecule = KNOWN_MOLECULES[key]
    return f"{molecule.name} ({molecule.collection})"


CATALOG_OPTION_INDEX = {
    _molecule_lookup_key(format_molecule_option(key)): key for key in KNOWN_MOLECULE_KEYS
}


def molecule_source_label(molecule) -> str:
    """Return a concise source label for a curated molecule."""
    if molecule.pubchem_cid:
        return f"{molecule.collection}, PubChem CID {molecule.pubchem_cid}"
    return molecule.collection


def lookup_known_molecule(query: str) -> Optional[MoleculeResolution]:
    """Resolve a curated molecule name or alias to SMILES."""
    key = MOLECULE_NAME_INDEX.get(_molecule_lookup_key(query.strip()))
    if not key:
        return None

    molecule = KNOWN_MOLECULES[key]
    return MoleculeResolution(molecule.name, molecule.smiles, molecule_source_label(molecule))


def looks_like_smiles(value: str) -> bool:
    """Best-effort guard so names go to lookup while obvious SMILES stay local."""
    value = value.strip()
    if not value or " " in value:
        return False
    if not SMILES_CHARS.fullmatch(value):
        return False

    unbracketed = re.sub(r"\[[^\]]+\]", "C", value)
    if not SMILES_ATOM_TOKEN.search(unbracketed):
        return False

    without_atoms = SMILES_ATOM_TOKEN.sub("", unbracketed)
    return not re.search(r"[A-Za-z]", without_atoms)


@lru_cache(maxsize=128)
def lookup_pubchem_molecule(name: str) -> Optional[MoleculeResolution]:
    """Resolve a molecule name through PubChem PUG REST."""
    name = name.strip()
    if not name:
        return None

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{quote(name, safe='')}/property/Title,CanonicalSMILES,IsomericSMILES/JSON"
    )
    request = Request(url, headers={"User-Agent": "morph-app/1.0"})

    try:
        with urlopen(request, timeout=8) as response:
            payload = json.load(response)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError, ValueError):
        return None

    properties = payload.get("PropertyTable", {}).get("Properties", [])
    if not properties:
        return None

    first = properties[0]
    smiles = (
        first.get("SMILES")
        or first.get("IsomericSMILES")
        or first.get("CanonicalSMILES")
        or first.get("ConnectivitySMILES")
    )
    if not smiles:
        return None

    cid = first.get("CID")
    source = f"PubChem CID {cid}" if cid else "PubChem"
    return MoleculeResolution(first.get("Title") or name, smiles, source)


def resolve_molecule_input(value: Optional[str]) -> MoleculeResolution:
    """Resolve a catalog selection, typed molecule name, or SMILES string."""
    value = (value or "").strip()
    if not value:
        raise ValueError("Choose a molecule or enter a valid SMILES string.")

    selected_key = KNOWN_MOLECULES.get(value) and value
    selected_key = selected_key or CATALOG_OPTION_INDEX.get(_molecule_lookup_key(value))
    if selected_key:
        molecule = KNOWN_MOLECULES[selected_key]
        return MoleculeResolution(molecule.name, molecule.smiles, molecule_source_label(molecule))

    known = lookup_known_molecule(value)
    if known:
        return known

    if looks_like_smiles(value):
        return MoleculeResolution("Custom SMILES", value, "Manual input")

    pubchem = lookup_pubchem_molecule(value)
    if pubchem:
        return pubchem

    raise ValueError(
        f"Could not look up '{value}'. Choose a catalog molecule or enter a valid SMILES string."
    )


def default_molecule(default_key: str) -> MoleculeResolution:
    molecule = KNOWN_MOLECULES[default_key]
    return MoleculeResolution(molecule.name, molecule.smiles, molecule_source_label(molecule))


def mols_to_smiles_text(mols) -> str:
    """Return one canonical SMILES string per molecule."""
    from rdkit import Chem

    return "\n".join(Chem.MolToSmiles(mol) for mol in mols)


def _accepted_lilly_ids(stdout: str):
    accepted_ids = set()
    for line in stdout.splitlines():
        fields = line.strip().split()
        if len(fields) >= 2:
            accepted_ids.add(fields[1])

    return accepted_ids


def filter_mols_with_lilly(mols) -> LillyFilterResult:
    """Filter generated molecules through Lilly_Medchem_Rules.rb -relaxed."""
    from rdkit import Chem

    lilly_rules = shutil.which(LILLY_MEDCHEM_RULES)
    if not lilly_rules:
        raise RuntimeError(f"{LILLY_MEDCHEM_RULES} was not found on PATH.")

    entries = []
    for index, mol in enumerate(mols):
        smiles = Chem.MolToSmiles(mol)
        entries.append((f"morph_{index:04d}", mol, smiles))

    if not entries:
        return LillyFilterResult([], "", 0)

    with tempfile.TemporaryDirectory(prefix="morph_lilly_") as tmpdir:
        input_path = Path(tmpdir) / "morph_output.smi"
        input_path.write_text(
            "".join(f"{smiles} {molecule_id}\n" for molecule_id, _mol, smiles in entries),
            encoding="utf-8",
        )
        completed = subprocess.run(
            [lilly_rules, "-relaxed", str(input_path)],
            cwd=tmpdir,
            text=True,
            capture_output=True,
            check=False,
            timeout=LILLY_FILTER_TIMEOUT_SECONDS,
        )

    if completed.returncode != 0:
        details = (completed.stderr or completed.stdout or "").strip()
        message = f"Lilly Medchem filter failed with exit code {completed.returncode}."
        if details:
            message = f"{message} {details}"
        raise RuntimeError(message)

    accepted_ids = _accepted_lilly_ids(completed.stdout)
    filtered_mols = [mol for molecule_id, mol, _smiles in entries if molecule_id in accepted_ids]

    return LillyFilterResult(
        filtered_mols,
        mols_to_smiles_text(filtered_mols),
        len(entries) - len(filtered_mols),
    )


def filter_morph_output_with_lilly(mols) -> LillyFilterResult:
    """Filter generated morph intermediates while preserving input endpoints."""
    mols = list(mols)
    if len(mols) <= 2:
        return LillyFilterResult(mols, mols_to_smiles_text(mols), 0)

    first_mol = mols[0]
    last_mol = mols[-1]
    filter_result = filter_mols_with_lilly(mols[1:-1])
    filtered_mols = [first_mol, *filter_result.mols, last_mol]
    return LillyFilterResult(
        filtered_mols,
        mols_to_smiles_text(filtered_mols),
        filter_result.rejected_count,
    )


def run_morph(smiles_1: str, smiles_2: str, apply_lilly_filter: bool = False):
    """Run morph between two SMILES. Returns (morph, smiles_text).

    morph.mol is the list of RDKit Mol objects for the pathway (use for SVG).
    Uses randomize=True so each run can produce a different pathway.
    Lilly_Medchem_Rules.rb filtering is optional and applies only after morphing.
    Raises ImportError if amsr is not installed.
    """
    import amsr

    s_tok = amsr.FromSmilesToTokens(smiles_1.strip(), randomize=True)
    t_tok = amsr.FromSmilesToTokens(smiles_2.strip(), randomize=True)
    morph = amsr.Morph(s_tok, t_tok)

    if apply_lilly_filter:
        filter_result = filter_morph_output_with_lilly(morph.mol)
        morph.mol = filter_result.mols
        setattr(morph, "lilly_rejected_count", filter_result.rejected_count)
        smiles_text = filter_result.smiles_text
    else:
        setattr(morph, "lilly_rejected_count", 0)
        smiles_text = mols_to_smiles_text(morph.mol)

    return morph, smiles_text


def mols_to_svgs(mols, mol_size: int = 180):
    """Render RDKit molecules to SVG strings."""
    from rdkit.Chem.Draw import rdMolDraw2D

    return [rdMolDraw2D.MolToSVG(mol, mol_size, mol_size) for mol in mols]


def molecule_input(st, label: str, key: str, default_key: str):
    """Render one searchable molecule field."""
    return st.text_input(
        label,
        value=format_molecule_option(default_key),
        key=key,
        placeholder="Name, SMILES, or catalog molecule",
        autocomplete="off",
    )


def main():
    import streamlit as st

    st.set_page_config(page_title="Morph molecules", layout="wide")

    st.markdown(
        """
<style>
    * { font-family: Arial, Helvetica, sans-serif !important; }
    .stMainBlockContainer { font-size: 14px; padding-top: 1rem !important; }
    .block-container { padding-top: 1rem !important; }
    h1 { font-size: 1.5rem !important; font-weight: 600 !important;
         margin-bottom: 0 !important; padding-bottom: 0 !important; }
    h2, h3, [data-testid="stSubheader"] {
        font-size: 1rem !important; font-weight: 700 !important;
        margin-top: 0.75rem !important; margin-bottom: 0.25rem !important; }
</style>
""",
        unsafe_allow_html=True,
    )

    with st.form("morph_form"):
        col1, col2 = st.columns(2)
        with col1:
            molecule_1_value = molecule_input(
                st, "From (name or SMILES)", "molecule_1", DEFAULT_MOLECULE_1_KEY
            )
        with col2:
            molecule_2_value = molecule_input(
                st, "To (name or SMILES)", "molecule_2", DEFAULT_MOLECULE_2_KEY
            )
        apply_lilly_filter = st.checkbox(
            "Filter morph output with Lilly Medchem Rules (-relaxed)",
            value=DEFAULT_APPLY_LILLY_FILTER,
        )
        submitted = st.form_submit_button("morph")

    if submitted:
        try:
            molecule_1 = resolve_molecule_input(molecule_1_value)
            molecule_2 = resolve_molecule_input(molecule_2_value)
        except ValueError as e:
            st.error(str(e))
            st.stop()

        st.caption(
            f"Morphing {molecule_1.name} ({molecule_1.source}) -> "
            f"{molecule_2.name} ({molecule_2.source})"
        )

        try:
            with st.spinner("Computing morph pathway..."):
                morph, smiles_text = run_morph(
                    molecule_1.smiles,
                    molecule_2.smiles,
                    apply_lilly_filter=apply_lilly_filter,
                )
        except Exception as e:
            st.error(f"Morph failed: {e}")
            import traceback

            st.code(traceback.format_exc())
            st.stop()

        if apply_lilly_filter:
            st.caption(
                f"Lilly Medchem Rules rejected "
                f"{getattr(morph, 'lilly_rejected_count', 0)} generated intermediates; "
                f"input endpoints were preserved."
            )

        # Molecules first (SVG in iframe), then SMILES below
        try:
            import streamlit.components.v1 as components

            COLS_PER_ROW = 4
            MOL_SIZE = 180
            svgs = mols_to_svgs(morph.mol, MOL_SIZE)
            if svgs:
                cells = "".join(f'<div style="flex: 0 0 auto;">{s}</div>' for s in svgs)
                html = f"""<!DOCTYPE html><html><body style="margin:0;padding:8px;">
                <div style="display:flex;flex-wrap:wrap;gap:12px;align-items:flex-start;">
                {cells}</div>
                </body></html>"""
                rows = (len(svgs) + COLS_PER_ROW - 1) // COLS_PER_ROW
                iframe_height = 24 + rows * (MOL_SIZE + 12)
                components.html(html, height=iframe_height, scrolling=False)
        except Exception as e:
            st.warning(f"Could not render molecules: {e}")

        # SMILES output: key per content so it updates on each morph (no stale state)
        st.text_area(
            "",
            value=smiles_text,
            height=200,
            key=f"morph_smiles_{hash(smiles_text) & 0xFFFFFFFF:X}",
            label_visibility="collapsed",
        )


if __name__ == "__main__":
    main()
