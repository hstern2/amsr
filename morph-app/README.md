# Morph

Streamlit app for **molecular morphing**: compute the minimum-edit pathway between two molecules. Based on the logic in `morph.ipynb` using the `amsr` library.

## Features

- Two searchable molecule fields with local name autocomplete; typed names and raw SMILES are still accepted
- Local quick lookup catalog with interesting natural products and FDA-approved drugs
- Autocomplete includes a fixed, lowercase vendored DrugCentral FDA-approved drug name list plus curated app molecules; numeric-leading names are excluded
- Catalog entries include PubChem CIDs for SMILES traceability; typed names still fall back to PubChem PUG REST
- Endpoint salts and mixtures are reduced to their largest connected component, excluding counterions
- **Morph** button tries 10 randomized pathways by default (configurable from 1 to 100) and displays the one with the most molecules remaining after filtering
- Optional filtering of generated morph intermediates through `Lilly_Medchem_Rules.rb -relaxed`, enabled by default
- Optional filtering of generated intermediates by RDKit cLogP and heteroatom count, matching the respective Muegge criteria in `~/mtrl`; both enabled by default with editable limits of -2 to 5 and at least two heteroatoms
- Output: rendered molecule grid plus downloadable pathway `.smi` and labeled, editable
  ChemDraw `.cdxml` files

## Requirements

- Python 3.x
- [Streamlit](https://streamlit.io/) (`pip install streamlit`)
- **amsr** — install from your environment (e.g. `pip install -e /path/to/amsr`); not on PyPI
- `Lilly_Medchem_Rules.rb` on `$PATH` to use the default intermediate filter

## Install

```bash
pip install -r requirements.txt
# Install amsr from your project/venv as needed
```

## Run the app

```bash
streamlit run morph_app.py
```

1. Enter **From** and **To** molecule names such as `epibatidine`, catalog labels, or raw SMILES.
2. Set the number of morphs to try and choose the Lilly, cLogP, and heteroatom filters. All three are enabled by default. The selected input endpoints are always preserved.
3. View the molecule grid. Use **Download .smi file** for a text representation, or
   **Download CDXML for ChemDraw** for an editable, labeled pathway laid out four structures
   per row.

## Test

```bash
pip install -r requirements.txt   # includes pytest
pytest
```

Runs the test suite. Tests that require `amsr` are skipped if the package is not installed. Using a virtual environment is recommended: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pytest`
