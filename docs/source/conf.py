# amsr/docs/source/conf.py

import os, sys
sys.path.insert(0, os.path.abspath("../../"))

project = "amsr"
author = "Harry Stern"
release = "0.1.4"

extensions = [
    "sphinx.ext.autodoc",       # <-- must precede sphinx_autodoc_typehints
    "sphinx.ext.napoleon",      # (optional) Google/NumPy docstrings
    "myst_parser",
    "sphinx_autodoc_typehints",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"
html_static_path = []  # keep empty unless you create _static

# (optional) nice defaults
autodoc_typehints = "description"
autodoc_default_options = {"members": True, "undoc-members": True}
