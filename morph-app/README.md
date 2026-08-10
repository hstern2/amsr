# Morph

Streamlit app for **molecular morphing**: compute the minimum-edit pathway between two molecules. Based on the logic in `morph.ipynb` using the `amsr` library.

## Features

- Two searchable molecule fields with local name autocomplete; typed names and raw SMILES are still accepted
- Local quick lookup catalog with interesting natural products and FDA-approved drugs
- Autocomplete includes a fixed, lowercase vendored DrugCentral FDA-approved drug name list plus curated app molecules; numeric-leading names are excluded
- Catalog entries include PubChem CIDs for SMILES traceability; typed names still fall back to PubChem PUG REST
- **Morph** button runs `amsr.Morph(...)` on the resolved SMILES to get the pathway
- Optional filtering of generated morph intermediates through `Lilly_Medchem_Rules.rb -relaxed`, enabled by default
- Output: downloadable pathway `.smi` file and rendered molecule grid

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
2. Leave the Lilly Medchem Rules intermediate filter enabled, or uncheck it to see the raw morph pathway. The selected input endpoints are not filtered.
3. View the molecule grid and use **Download .smi file** to save endpoint SMILES with displayed name/source labels and numbered intermediate records.

## Test

```bash
pip install -r requirements.txt   # includes pytest
pytest
```

Runs the test suite. Tests that require `amsr` are skipped if the package is not installed. Using a virtual environment is recommended: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pytest`
