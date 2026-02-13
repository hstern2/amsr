#!/usr/bin/env python3
"""Morph: Streamlit app for molecular morph (two SMILES → pathway). Based on morph.ipynb."""

import streamlit as st

# Default SMILES (epothilone A, rocaglamide)
DEFAULT_SMILES_1 = "CC1=C(C(=O)CC(C)C1CC=CC2CC(CC(=O)O2)C)C"
DEFAULT_SMILES_2 = "COC1=CC(=C2C(=C1)C(=O)C3CC(CC(O3)C=C(C)C)O2)C4=CC=CC=N4"


def run_morph(smiles_1: str, smiles_2: str):
    """Run morph between two SMILES. Returns (morph, smiles_text).

    morph.mol is the list of RDKit Mol objects for the pathway (use for SVG).
    Uses randomize=True so each run can produce a different pathway.
    Raises ImportError if amsr is not installed.
    """
    from rdkit import Chem

    import amsr

    s_tok = amsr.FromSmilesToTokens(smiles_1.strip(), randomize=True)
    t_tok = amsr.FromSmilesToTokens(smiles_2.strip(), randomize=True)
    morph = amsr.Morph(s_tok, t_tok)

    smiles_text = "\n".join(Chem.MolToSmiles(m) for m in morph.mol)
    return morph, smiles_text


st.set_page_config(page_title="Morph", layout="wide")

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
    smiles_1 = st.text_input(
        "Molecule 1 SMILES", value="", placeholder=DEFAULT_SMILES_1, key="smiles1"
    )
    smiles_2 = st.text_input(
        "Molecule 2 SMILES", value="", placeholder=DEFAULT_SMILES_2, key="smiles2"
    )
    submitted = st.form_submit_button("morph")

if submitted:
    # Form widget values are available on the run after submit
    s1 = (smiles_1 or "").strip() or DEFAULT_SMILES_1
    s2 = (smiles_2 or "").strip() or DEFAULT_SMILES_2

    try:
        with st.spinner("Computing morph pathway..."):
            morph, smiles_text = run_morph(s1, s2)
    except Exception as e:
        st.error(f"Morph failed: {e}")
        import traceback

        st.code(traceback.format_exc())
        st.stop()

    # Molecules first (SVG in iframe), then SMILES below
    try:
        import streamlit.components.v1 as components
        from rdkit.Chem import Draw

        COLS_PER_ROW = 4
        MOL_SIZE = 180
        svgs = [Draw.MolToSVG(mol, MOL_SIZE, MOL_SIZE) for mol in morph.mol]
        if svgs:
            cells = "".join(f'<div style="flex: 0 0 auto;">{s}</div>' for s in svgs)
            html = f"""<!DOCTYPE html><html><body style="margin:0;padding:8px;">
            <div style="display:flex;flex-wrap:wrap;gap:12px;align-items:flex-start;">{cells}</div>
            </body></html>"""
            rows = (len(svgs) + COLS_PER_ROW - 1) // COLS_PER_ROW
            iframe_height = 24 + rows * (MOL_SIZE + 12)
            components.html(html, height=iframe_height, scrolling=False)
    except Exception as e:
        st.warning(f"Could not render molecules: {e}")

    # SMILES output — key per content so it updates on each morph (no stale state)
    st.text_area(
        "",
        value=smiles_text,
        height=200,
        key=f"morph_smiles_{hash(smiles_text) & 0xFFFFFFFF:X}",
        label_visibility="collapsed",
    )
