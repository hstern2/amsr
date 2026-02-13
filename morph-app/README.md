# Morph

Streamlit app for **molecular morphing**: compute the minimum-edit pathway between two molecules (SMILES). Based on the logic in `morph.ipynb` using the `amsr` library.

## Features

- Two SMILES inputs with optional names (defaults: epithilone, rocaglamide)
- **Morph** button runs `amsr.Morph.fromSmiles(s, t)` to get the pathway
- Output: pathway SMILES in a text box (one per line) and rendered molecule grid

## Requirements

- Python 3.x
- [Streamlit](https://streamlit.io/) (`pip install streamlit`)
- **amsr** — install from your environment (e.g. `pip install -e /path/to/amsr`); not on PyPI

## Install

```bash
pip install -r requirements.txt
# Install amsr from your project/venv as needed
```

## Run the app

```bash
streamlit run morph_app.py
```

1. Optionally edit names and SMILES (defaults: Epithilone, Rocaglamide).
2. Click **Morph** to compute the pathway.
3. View SMILES in the text area and the molecule grid below.

## Test

```bash
pip install -r requirements.txt   # includes pytest
pytest
```

Runs the test suite. Tests that require `amsr` are skipped if the package is not installed. Using a virtual environment is recommended: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pytest`
