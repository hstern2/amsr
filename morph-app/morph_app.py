#!/usr/bin/env python3
"""Morph: Streamlit app for molecular morph (two SMILES -> pathway). Based on morph.ipynb."""

import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from html import escape
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
class FilterResult:
    mols: list
    smiles_text: str
    rejected_count: int


LillyFilterResult = FilterResult


DEFAULT_MOLECULE_1_KEY = "ibogaine"
DEFAULT_MOLECULE_2_KEY = "epibatidine"
DEFAULT_SMILES_1 = KNOWN_MOLECULES[DEFAULT_MOLECULE_1_KEY].smiles
DEFAULT_SMILES_2 = KNOWN_MOLECULES[DEFAULT_MOLECULE_2_KEY].smiles
DEFAULT_APPLY_LILLY_FILTER = True
DEFAULT_APPLY_CLOGP_FILTER = True
DEFAULT_APPLY_HETEROATOM_FILTER = True
DEFAULT_MORPH_COUNT = 10
MIN_CLOGP = -2.0
MAX_CLOGP = 5.0
MIN_HETEROATOMS = 2
LILLY_MEDCHEM_RULES = "Lilly_Medchem_Rules.rb"
LILLY_FILTER_TIMEOUT_SECONDS = 120
GITHUB_REPO_URL = "https://github.com/hstern2/amsr"
AUTOCOMPLETE_NAMES_PATH = Path(__file__).with_name("molecule_autocomplete_names.txt")


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


def read_autocomplete_names(path: Path) -> tuple[str, ...]:
    """Read vendored autocomplete names from a non-executable text file."""
    names = []
    for line in path.read_text(encoding="ascii").splitlines():
        name = line.strip()
        if name and not name.startswith("#"):
            names.append(name)
    return tuple(names)


def molecule_autocomplete_options() -> tuple[str, ...]:
    """Return name-only autocomplete suggestions for the molecule inputs."""
    names_by_key = {
        name.casefold(): name.lower() for name in read_autocomplete_names(AUTOCOMPLETE_NAMES_PATH)
    }
    for molecule in KNOWN_MOLECULES.values():
        names_by_key[molecule.name.casefold()] = molecule.name.lower()
        for alias in molecule.aliases:
            names_by_key.setdefault(alias.casefold(), alias.lower())

    return tuple(sorted(names_by_key.values(), key=lambda name: name.casefold()))


MOLECULE_AUTOCOMPLETE_OPTIONS = molecule_autocomplete_options()


def autocomplete_matches(query: Optional[str], limit: int = 6) -> tuple[str, ...]:
    """Return compact name suggestions without constraining free-form input."""
    query = (query or "").strip().lower()
    if len(query) < 2 or looks_like_smiles(query):
        return ()

    prefix_matches = [
        option
        for option in MOLECULE_AUTOCOMPLETE_OPTIONS
        if option.startswith(query) and option != query
    ]
    if len(prefix_matches) >= limit:
        return tuple(prefix_matches[:limit])

    prefix_keys = set(prefix_matches)
    contains_matches = [
        option
        for option in MOLECULE_AUTOCOMPLETE_OPTIONS
        if query in option and option != query and option not in prefix_keys
    ]
    return tuple([*prefix_matches, *contains_matches][:limit])


def set_molecule_input_value(key: str, value: str):
    import streamlit as st

    st.session_state[key] = value


def request_morph_from_input():
    import streamlit as st

    st.session_state["morph_requested"] = True


def morph_input_snapshot(
    molecule_1_value: Optional[str],
    molecule_2_value: Optional[str],
    apply_lilly_filter: bool,
    apply_clogp_filter: bool = DEFAULT_APPLY_CLOGP_FILTER,
    apply_heteroatom_filter: bool = DEFAULT_APPLY_HETEROATOM_FILTER,
    morph_count: int = DEFAULT_MORPH_COUNT,
    min_clogp: float = MIN_CLOGP,
    max_clogp: float = MAX_CLOGP,
    min_heteroatoms: int = MIN_HETEROATOMS,
) -> tuple[str, str, bool, bool, bool, int, float, float, int]:
    """Return the committed input state that should define one morph request."""
    return (
        (molecule_1_value or "").strip(),
        (molecule_2_value or "").strip(),
        apply_lilly_filter,
        apply_clogp_filter,
        apply_heteroatom_filter,
        morph_count,
        min_clogp,
        max_clogp,
        min_heteroatoms,
    )


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


def largest_connected_component_smiles(smiles: str) -> str:
    """Return the largest connected component of a SMILES endpoint."""
    from rdkit import Chem

    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Enter a valid SMILES string.")

    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not fragments:
        raise ValueError("Enter a valid SMILES string.")
    if len(fragments) == 1:
        return smiles

    largest_fragment = max(fragments, key=lambda fragment: fragment.GetNumAtoms())
    return Chem.MolToSmiles(largest_fragment, isomericSmiles=True)


def _without_counterions(molecule: MoleculeResolution) -> MoleculeResolution:
    """Keep only the largest connected component of a resolved endpoint."""
    return MoleculeResolution(
        molecule.name,
        largest_connected_component_smiles(molecule.smiles),
        molecule.source,
    )


def resolve_molecule_input(value: Optional[str]) -> MoleculeResolution:
    """Resolve an endpoint and discard all but its largest connected component."""
    value = (value or "").strip()
    if not value:
        raise ValueError("Choose a molecule or enter a valid SMILES string.")

    selected_key = KNOWN_MOLECULES.get(value) and value
    selected_key = selected_key or CATALOG_OPTION_INDEX.get(_molecule_lookup_key(value))
    if selected_key:
        molecule = KNOWN_MOLECULES[selected_key]
        return _without_counterions(
            MoleculeResolution(molecule.name, molecule.smiles, molecule_source_label(molecule))
        )

    known = lookup_known_molecule(value)
    if known:
        return _without_counterions(known)

    if looks_like_smiles(value):
        return _without_counterions(MoleculeResolution("Custom SMILES", value, "Manual input"))

    pubchem = lookup_pubchem_molecule(value)
    if pubchem:
        return _without_counterions(pubchem)

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


def _clean_smi_title(title: str) -> str:
    """Keep SMI record titles on one tab-separated line."""
    return re.sub(r"[\r\n\t]+", " ", title).strip()


def endpoint_smi_title(endpoint_number: int, molecule: Optional[MoleculeResolution]) -> str:
    if molecule is None:
        return f"endpoint_{endpoint_number}"
    return _clean_smi_title(f"endpoint_{endpoint_number}: {molecule.name} ({molecule.source})")


def pathway_to_smi_text(
    mols,
    molecule_1: Optional[MoleculeResolution] = None,
    molecule_2: Optional[MoleculeResolution] = None,
) -> str:
    """Return SMI file content with stable names for endpoints and intermediates."""
    from rdkit import Chem

    mols = list(mols)
    rows = []
    for index, mol in enumerate(mols):
        if index == 0:
            name = endpoint_smi_title(1, molecule_1)
        elif index == len(mols) - 1:
            name = endpoint_smi_title(2, molecule_2)
        else:
            name = f"intermediate_{index:03d}"
        rows.append(f"{Chem.MolToSmiles(mol)}\t{name}")

    return "\n".join(rows) + ("\n" if rows else "")


def _pathway_molecule_label(
    index: int,
    molecule_count: int,
    molecule_1: Optional[MoleculeResolution],
    molecule_2: Optional[MoleculeResolution],
) -> str:
    """Return a short human-readable caption for a pathway structure."""
    if index == 0:
        name = molecule_1.name if molecule_1 else "Endpoint 1"
        return _clean_smi_title(f"{index + 1}. From: {name}")
    if index == molecule_count - 1:
        name = molecule_2.name if molecule_2 else "Endpoint 2"
        return _clean_smi_title(f"{index + 1}. To: {name}")
    return f"{index + 1}. Intermediate {index:03d}"


def _format_cdxml_number(value: float) -> str:
    """Format a CDXML drawing coordinate compactly and consistently."""
    return f"{value:.2f}".rstrip("0").rstrip(".") or "0"


def _cdxml_fragment_bounds(fragment) -> tuple[float, float, float, float]:
    """Return the atom-coordinate bounds for an ElementTree CDXML fragment."""
    positions = []
    for node in fragment.iter("n"):
        position = node.get("p")
        if position:
            x, y = (float(value) for value in position.split())
            positions.append((x, y))

    if not positions:
        raise ValueError("RDKit produced a CDXML fragment without atom coordinates.")

    xs, ys = zip(*positions)
    return min(xs), min(ys), max(xs), max(ys)


def _translate_cdxml_fragment(fragment, dx: float, dy: float) -> None:
    """Translate all atom positions in a CDXML fragment in place."""
    for node in fragment.iter("n"):
        position = node.get("p")
        if not position:
            continue
        x, y = (float(value) for value in position.split())
        node.set("p", f"{_format_cdxml_number(x + dx)} {_format_cdxml_number(y + dy)}")


def _renumber_cdxml_fragment(fragment, next_id: int) -> int:
    """Give a copied RDKit fragment document-unique object IDs."""
    id_map = {}
    for element in fragment.iter():
        old_id = element.get("id")
        if old_id is not None:
            id_map[old_id] = str(next_id)
            element.set("id", str(next_id))
            next_id += 1

    # These are the object-ID references emitted by RDKit's molecule writer.
    for element in fragment.iter():
        for attribute in ("B", "E", "BondOrdering"):
            references = element.get(attribute)
            if references:
                element.set(
                    attribute,
                    " ".join(id_map.get(reference, reference) for reference in references.split()),
                )

    return next_id


def pathway_to_cdxml(
    mols,
    molecule_1: Optional[MoleculeResolution] = None,
    molecule_2: Optional[MoleculeResolution] = None,
    columns: int = 4,
) -> str:
    """Return the pathway as a labeled, editable ChemDraw CDXML document."""
    from xml.etree import ElementTree

    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    mols = list(mols)
    if not mols:
        raise ValueError("Cannot make a CDXML document from an empty pathway.")
    if columns < 1:
        raise ValueError("CDXML layout columns must be at least one.")
    if not hasattr(Chem, "MolToCDXMLBlock") or not Chem.HasChemDrawCDXSupport():
        raise RuntimeError(
            "This RDKit installation does not include ChemDraw CDXML writing support."
        )

    fragments = []
    bounds = []
    for mol in mols:
        drawing_mol = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(drawing_mol, canonOrient=True, clearConfs=True)
        molecule_document = ElementTree.fromstring(Chem.MolToCDXMLBlock(drawing_mol))
        fragment = molecule_document.find("./page/fragment")
        if fragment is None:
            raise ValueError("RDKit did not produce a molecule fragment in its CDXML output.")
        fragments.append(fragment)
        bounds.append(_cdxml_fragment_bounds(fragment))

    structure_width = max(max_x - min_x for min_x, _min_y, max_x, _max_y in bounds)
    structure_height = max(max_y - min_y for _min_x, min_y, _max_x, max_y in bounds)
    cell_width = max(180.0, structure_width + 72.0)
    cell_height = max(160.0, structure_height + 96.0)
    page_padding = 36.0
    used_columns = min(columns, len(fragments))
    row_count = (len(fragments) + columns - 1) // columns
    page_width = 2 * page_padding + used_columns * cell_width
    page_height = 2 * page_padding + row_count * cell_height

    root = ElementTree.Element(
        "CDXML",
        {
            "id": "1",
            "CreationProgram": "AMSR Morph",
            "Name": "morph_pathway.cdxml",
            "BondLength": "28.8",
            "LabelFont": "3",
            "LabelSize": "10",
            "LabelFace": "96",
            "CaptionFont": "3",
            "CaptionSize": "10",
            "CaptionFace": "0",
        },
    )
    font_table = ElementTree.SubElement(root, "fonttable")
    ElementTree.SubElement(
        font_table,
        "font",
        {"id": "3", "charset": "iso-8859-1", "name": "Arial"},
    )
    page = ElementTree.SubElement(
        root,
        "page",
        {
            "id": "2",
            "BoundingBox": (
                f"0 0 {_format_cdxml_number(page_width)} {_format_cdxml_number(page_height)}"
            ),
        },
    )

    next_id = 10
    for index, (fragment, (min_x, min_y, max_x, max_y)) in enumerate(zip(fragments, bounds)):
        row, column = divmod(index, columns)
        cell_left = page_padding + column * cell_width
        cell_top = page_padding + row * cell_height
        target_x = cell_left + cell_width / 2
        target_y = cell_top + structure_height / 2
        _translate_cdxml_fragment(
            fragment,
            target_x - (min_x + max_x) / 2,
            target_y - (min_y + max_y) / 2,
        )
        next_id = _renumber_cdxml_fragment(fragment, next_id)
        page.append(fragment)

        label_y = cell_top + structure_height + 28.0
        label = ElementTree.SubElement(
            page,
            "t",
            {
                "id": str(next_id),
                "p": f"{_format_cdxml_number(target_x)} {_format_cdxml_number(label_y)}",
                "BoundingBox": (
                    f"{_format_cdxml_number(cell_left + 8)} "
                    f"{_format_cdxml_number(label_y - 12)} "
                    f"{_format_cdxml_number(cell_left + cell_width - 8)} "
                    f"{_format_cdxml_number(label_y + 4)}"
                ),
                "LabelJustification": "Center",
                "Justification": "Center",
                "InterpretChemically": "no",
            },
        )
        next_id += 1
        style = ElementTree.SubElement(
            label,
            "s",
            {"font": "3", "size": "10", "face": "0"},
        )
        style.text = _pathway_molecule_label(index, len(fragments), molecule_1, molecule_2)

    ElementTree.indent(root, space="  ")
    body = ElementTree.tostring(root, encoding="unicode", short_empty_elements=True)
    return (
        '<?xml version="1.0" encoding="UTF-8" ?>\n'
        '<!DOCTYPE CDXML SYSTEM "https://static.chemistry.revvitycloud.com/cdxml/CDXML.dtd" >\n'
        f"{body}\n"
    )


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


def filter_morph_output_by_clogp(
    mols, min_clogp: float = MIN_CLOGP, max_clogp: float = MAX_CLOGP
) -> FilterResult:
    """Keep intermediates in mtrl's Muegge cLogP range, preserving endpoints."""
    from rdkit import Chem
    from rdkit.Chem import Crippen

    mols = list(mols)
    if len(mols) <= 2:
        return FilterResult(mols, mols_to_smiles_text(mols), 0)

    kept = [mols[0]]
    for mol in mols[1:-1]:
        clogp = Crippen.MolLogP(Chem.RemoveHs(mol))
        if math.isfinite(clogp) and min_clogp <= clogp <= max_clogp:
            kept.append(mol)
    kept.append(mols[-1])
    return FilterResult(kept, mols_to_smiles_text(kept), len(mols) - len(kept))


def filter_morph_output_by_heteroatoms(
    mols, min_heteroatoms: int = MIN_HETEROATOMS
) -> FilterResult:
    """Keep intermediates with at least two heteroatoms, preserving endpoints."""
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    mols = list(mols)
    if len(mols) <= 2:
        return FilterResult(mols, mols_to_smiles_text(mols), 0)

    kept = [mols[0]]
    for mol in mols[1:-1]:
        if rdMolDescriptors.CalcNumHeteroatoms(Chem.RemoveHs(mol)) >= min_heteroatoms:
            kept.append(mol)
    kept.append(mols[-1])
    return FilterResult(kept, mols_to_smiles_text(kept), len(mols) - len(kept))


def run_morph(
    smiles_1: str,
    smiles_2: str,
    apply_lilly_filter: bool = False,
    apply_clogp_filter: bool = False,
    apply_heteroatom_filter: bool = False,
    min_clogp: float = MIN_CLOGP,
    max_clogp: float = MAX_CLOGP,
    min_heteroatoms: int = MIN_HETEROATOMS,
):
    """Run morph between two SMILES. Returns (morph, smiles_text).

    morph.mol is the list of RDKit Mol objects for the pathway (use for SVG).
    Uses randomize=True so each run can produce a different pathway.
    Optional filters apply only to generated intermediates after morphing.
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
    else:
        setattr(morph, "lilly_rejected_count", 0)

    if apply_clogp_filter:
        filter_result = filter_morph_output_by_clogp(morph.mol, min_clogp, max_clogp)
        morph.mol = filter_result.mols
        setattr(morph, "clogp_rejected_count", filter_result.rejected_count)
    else:
        setattr(morph, "clogp_rejected_count", 0)

    if apply_heteroatom_filter:
        filter_result = filter_morph_output_by_heteroatoms(morph.mol, min_heteroatoms)
        morph.mol = filter_result.mols
        setattr(morph, "heteroatom_rejected_count", filter_result.rejected_count)
    else:
        setattr(morph, "heteroatom_rejected_count", 0)

    smiles_text = mols_to_smiles_text(morph.mol)

    return morph, smiles_text


def run_best_morph(
    smiles_1: str,
    smiles_2: str,
    morph_count: int = DEFAULT_MORPH_COUNT,
    apply_lilly_filter: bool = DEFAULT_APPLY_LILLY_FILTER,
    apply_clogp_filter: bool = DEFAULT_APPLY_CLOGP_FILTER,
    apply_heteroatom_filter: bool = DEFAULT_APPLY_HETEROATOM_FILTER,
    min_clogp: float = MIN_CLOGP,
    max_clogp: float = MAX_CLOGP,
    min_heteroatoms: int = MIN_HETEROATOMS,
):
    """Run independent randomized morphs and return the most retained pathway."""
    if morph_count < 1:
        raise ValueError("Number of morphs must be at least one.")
    if apply_clogp_filter and min_clogp > max_clogp:
        raise ValueError("Minimum cLogP must be no greater than maximum cLogP.")
    if apply_heteroatom_filter and min_heteroatoms < 0:
        raise ValueError("Minimum heteroatom count cannot be negative.")

    best = None
    best_index = 0
    for index in range(1, morph_count + 1):
        candidate = run_morph(
            smiles_1,
            smiles_2,
            apply_lilly_filter=apply_lilly_filter,
            apply_clogp_filter=apply_clogp_filter,
            apply_heteroatom_filter=apply_heteroatom_filter,
            min_clogp=min_clogp,
            max_clogp=max_clogp,
            min_heteroatoms=min_heteroatoms,
        )
        if best is None or len(candidate[0].mol) > len(best[0].mol):
            best = candidate
            best_index = index

    assert best is not None
    morph, smiles_text = best
    return morph, smiles_text, best_index


def kekulized_mol_for_drawing(mol):
    """Return a drawing-only molecule copy with aromatic flags cleared."""
    from rdkit import Chem

    drawing_mol = Chem.Mol(mol)
    try:
        Chem.Kekulize(drawing_mol, clearAromaticFlags=True)
    except Exception:
        return None
    return drawing_mol


def mol_svg_pairs(mols, mol_size: int = 180) -> list[tuple[object, str]]:
    """Render molecules that can be kekulized to SVG strings."""
    from rdkit.Chem.Draw import rdMolDraw2D

    pairs = []
    for mol in mols:
        drawing_mol = kekulized_mol_for_drawing(mol)
        if drawing_mol is None:
            continue
        pairs.append((mol, rdMolDraw2D.MolToSVG(drawing_mol, mol_size, mol_size)))
    return pairs


def mols_to_svgs(mols, mol_size: int = 180):
    """Render molecules that can be kekulized to SVG strings."""
    return [svg for _mol, svg in mol_svg_pairs(mols, mol_size)]


@lru_cache(maxsize=1)
def _sascorer():
    from rdkit import RDConfig

    sa_score_dir = str(Path(RDConfig.RDContribDir) / "SA_Score")
    if sa_score_dir not in sys.path:
        sys.path.append(sa_score_dir)

    import sascorer

    return sascorer


def molecule_properties(mol) -> list[tuple[str, str]]:
    """Return the same 2D descriptor set shown in the AMSR Flask app."""
    from rdkit.Chem.Crippen import MolLogP
    from rdkit.Chem.Descriptors import TPSA, MolWt
    from rdkit.Chem.Lipinski import (
        HeavyAtomCount,
        NHOHCount,
        NOCount,
        NumRotatableBonds,
    )
    from rdkit.Chem.QED import qed

    hac = HeavyAtomCount(mol)
    mw = MolWt(mol)
    clogp = MolLogP(mol)
    # Use the original/traditional Lipinski Rule-of-Five HBA/HBD definitions.
    hbd = NHOHCount(mol)
    hba = NOCount(mol)
    n_rot_bonds = NumRotatableBonds(mol)
    passes_ro5 = mw <= 500 and clogp <= 5 and hbd <= 5 and hba <= 10

    return [
        ("QED score", f"{qed(mol):.3f}"),
        ("TPSA", f"{TPSA(mol):.3f} &#8491;<sup>2</sup>"),
        ("SA score", f"{_sascorer().calculateScore(mol):.3f}"),
        ("Heavy atom count", str(hac)),
        ("Molecular weight", f"{mw:.2f} Da"),
        ("LogP", f"{clogp:.3f}"),
        ("H-bond donors", str(hbd)),
        ("H-bond acceptors", str(hba)),
        ("Rotatable bonds", str(n_rot_bonds)),
        ("Passes Rule of 5", "Yes" if passes_ro5 else "No"),
    ]


def molecule_hover_grid_html(mols, mol_size: int = 180) -> str:
    """Render molecules as an HTML grid with descriptor tooltips on hover."""
    cells = []
    for index, (mol, svg) in enumerate(mol_svg_pairs(mols, mol_size)):
        properties = "".join(
            f"<div><strong>{escape(label)}:</strong> {value}</div>"
            for label, value in molecule_properties(mol)
        )
        cells.append(
            f"""
            <div class="mol-card" tabindex="0" aria-describedby="mol-props-{index}">
                {svg}
                <div id="mol-props-{index}" class="mol-tooltip" role="tooltip">
                    {properties}
                </div>
            </div>
            """
        )

    return f"""<!DOCTYPE html>
<html>
<head>
<style>
    body {{
        margin: 0;
        padding: 8px;
        font-family: Arial, Helvetica, sans-serif;
        font-size: 12px;
    }}
    .mol-grid {{
        display: flex;
        flex-wrap: wrap;
        gap: 12px;
        align-items: flex-start;
    }}
    .mol-card {{
        position: relative;
        flex: 0 0 auto;
        width: {mol_size}px;
        height: {mol_size}px;
        outline: none;
    }}
    .mol-card svg {{
        display: block;
        width: {mol_size}px;
        height: {mol_size}px;
    }}
    .mol-tooltip {{
        position: fixed;
        top: 8px;
        left: -9999px;
        box-sizing: border-box;
        width: 240px;
        max-width: calc(100vw - 16px);
        max-height: calc(100vh - 16px);
        overflow: auto;
        border: 1px solid #777;
        background: rgba(255, 255, 255, 0.96);
        color: #111;
        line-height: 1.35;
        padding: 10px;
        opacity: 0;
        visibility: hidden;
        pointer-events: none;
        z-index: 10;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.16);
        transition: opacity 100ms ease;
    }}
    .mol-card:hover .mol-tooltip,
    .mol-card:focus .mol-tooltip {{
        opacity: 1;
        visibility: visible;
    }}
</style>
</head>
<body>
    <div class="mol-grid">
        {"".join(cells)}
    </div>
    <script>
        function placeTooltip(card) {{
            const tooltip = card.querySelector('.mol-tooltip');
            if (!tooltip) return;

            const gap = 10;
            const padding = 8;
            const rect = card.getBoundingClientRect();
            const viewWidth = document.documentElement.clientWidth;
            const viewHeight = document.documentElement.clientHeight;
            const tooltipWidth = tooltip.offsetWidth;
            const tooltipHeight = tooltip.offsetHeight;
            const rightX = rect.right + gap;
            const leftX = rect.left - tooltipWidth - gap;
            const belowY = rect.bottom + gap;
            const aboveY = rect.top - tooltipHeight - gap;
            let x;
            let y = rect.top;

            if (rightX + tooltipWidth <= viewWidth - padding) {{
                x = rightX;
            }} else if (leftX >= padding) {{
                x = leftX;
            }} else {{
                x = Math.min(
                    Math.max(padding, rect.left + (rect.width - tooltipWidth) / 2),
                    viewWidth - tooltipWidth - padding
                );
                if (belowY + tooltipHeight <= viewHeight - padding) {{
                    y = belowY;
                }} else if (aboveY >= padding) {{
                    y = aboveY;
                }}
            }}

            if (y + tooltipHeight > viewHeight - padding) {{
                y = viewHeight - tooltipHeight - padding;
            }}
            y = Math.max(padding, y);

            tooltip.style.left = `${{Math.max(padding, x)}}px`;
            tooltip.style.top = `${{y}}px`;
        }}

        document.querySelectorAll('.mol-card').forEach((card) => {{
            card.addEventListener('mouseenter', () => placeTooltip(card));
            card.addEventListener('focus', () => placeTooltip(card));
        }});
    </script>
</body>
</html>"""


def molecule_input(st, label: str, key: str):
    """Render one free-form molecule field with local suggestions."""
    value = st.text_input(
        label,
        key=key,
        placeholder="Name or SMILES",
        autocomplete="off",
        on_change=request_morph_from_input,
    )

    suggestions = autocomplete_matches(value)
    if suggestions:
        cols = st.columns(min(len(suggestions), 3))
        for index, suggestion in enumerate(suggestions):
            with cols[index % len(cols)]:
                st.button(
                    suggestion,
                    key=f"{key}_suggestion_{index}_{_molecule_lookup_key(suggestion)}",
                    type="tertiary",
                    on_click=set_molecule_input_value,
                    args=(key, suggestion),
                )

    return value


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
    .morph-appbar {
        display: flex;
        align-items: center;
        gap: 0.75rem;
        margin-bottom: 0.75rem;
    }
    .morph-appbar-title {
        font-size: 1.5rem;
        font-weight: 600;
        line-height: 1.2;
    }
    .morph-github-link {
        color: inherit !important;
        display: inline-flex;
        align-items: center;
        line-height: 1;
        text-decoration: none;
    }
    .morph-github-link svg {
        width: 18px;
        height: 18px;
        fill: currentColor;
    }
</style>
""",
        unsafe_allow_html=True,
    )

    st.markdown(
        f"""
<div class="morph-appbar">
    <span class="morph-appbar-title">Morph molecules</span>
    <a class="morph-github-link" href="{GITHUB_REPO_URL}" aria-label="GitHub"
       title="source code on GitHub" target="_blank" rel="noopener noreferrer">
        <svg viewBox="0 0 16 16" aria-hidden="true">
            <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38
            0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52
            -.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2
            -3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82
            .64-.18 1.32-.27 2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12
            .51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48
            0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z"/>
        </svg>
    </a>
</div>
""",
        unsafe_allow_html=True,
    )

    col1, col2 = st.columns(2)
    with col1:
        molecule_1_value = molecule_input(st, "From (name or SMILES)", "molecule_1")
    with col2:
        molecule_2_value = molecule_input(st, "To (name or SMILES)", "molecule_2")
    apply_lilly_filter = st.checkbox(
        "Filter morph output with Lilly Medchem Rules (-relaxed)",
        value=DEFAULT_APPLY_LILLY_FILTER,
    )
    apply_clogp_filter = st.checkbox(
        "Filter intermediates by cLogP",
        value=DEFAULT_APPLY_CLOGP_FILTER,
    )
    clogp_columns = st.columns(2)
    with clogp_columns[0]:
        min_clogp = st.number_input("Minimum cLogP", value=MIN_CLOGP, step=0.5, key="min_clogp")
    with clogp_columns[1]:
        max_clogp = st.number_input("Maximum cLogP", value=MAX_CLOGP, step=0.5, key="max_clogp")
    apply_heteroatom_filter = st.checkbox(
        "Filter intermediates by minimum heteroatoms",
        value=DEFAULT_APPLY_HETEROATOM_FILTER,
    )
    min_heteroatoms = st.number_input(
        "Minimum heteroatoms",
        min_value=0,
        value=MIN_HETEROATOMS,
        step=1,
        key="min_heteroatoms",
    )
    morph_count = st.number_input(
        "Number of morphs to try",
        min_value=1,
        max_value=100,
        value=DEFAULT_MORPH_COUNT,
        step=1,
        key="morph_count",
    )
    current_morph_inputs = morph_input_snapshot(
        molecule_1_value,
        molecule_2_value,
        apply_lilly_filter,
        apply_clogp_filter,
        apply_heteroatom_filter,
        morph_count,
        min_clogp,
        max_clogp,
        min_heteroatoms,
    )
    previous_morph_inputs = st.session_state.setdefault(
        "last_morph_input_snapshot", current_morph_inputs
    )
    inputs_changed = current_morph_inputs != previous_morph_inputs
    morph_requested = st.session_state.pop("morph_requested", False)
    submitted = st.button("morph") or morph_requested or inputs_changed

    if submitted:
        st.session_state["last_morph_input_snapshot"] = current_morph_inputs
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
            with st.spinner(f"Computing {morph_count} morph pathway(s)..."):
                morph, _smiles_text, selected_index = run_best_morph(
                    molecule_1.smiles,
                    molecule_2.smiles,
                    morph_count=morph_count,
                    apply_lilly_filter=apply_lilly_filter,
                    apply_clogp_filter=apply_clogp_filter,
                    apply_heteroatom_filter=apply_heteroatom_filter,
                    min_clogp=min_clogp,
                    max_clogp=max_clogp,
                    min_heteroatoms=min_heteroatoms,
                )
        except Exception as e:
            st.error(f"Morph failed: {e}")
            import traceback

            st.code(traceback.format_exc())
            st.stop()

        st.caption(
            f"Showing morph {selected_index} of {morph_count}: "
            f"{len(morph.mol)} molecules retained after filtering."
        )
        if apply_lilly_filter:
            st.caption(
                f"Lilly Medchem Rules rejected "
                f"{getattr(morph, 'lilly_rejected_count', 0)} generated intermediates; "
                f"input endpoints were preserved."
            )
        if apply_clogp_filter:
            st.caption(
                f"cLogP outside {min_clogp:g} to {max_clogp:g} rejected "
                f"{getattr(morph, 'clogp_rejected_count', 0)} generated intermediates; "
                f"input endpoints were preserved."
            )
        if apply_heteroatom_filter:
            st.caption(
                f"Fewer than {min_heteroatoms} heteroatoms rejected "
                f"{getattr(morph, 'heteroatom_rejected_count', 0)} generated intermediates; "
                f"input endpoints were preserved."
            )

        # Molecules first (SVG in iframe), then downloadable SMI below
        try:
            import streamlit.components.v1 as components

            COLS_PER_ROW = 4
            MOL_SIZE = 180
            molecule_count = len(morph.mol)
            if molecule_count:
                html = molecule_hover_grid_html(morph.mol, MOL_SIZE)
                rows = (molecule_count + COLS_PER_ROW - 1) // COLS_PER_ROW
                iframe_height = 24 + rows * (MOL_SIZE + 12)
                components.html(html, height=iframe_height, scrolling=False)
        except Exception as e:
            st.warning(f"Could not render molecules: {e}")

        smi_text = pathway_to_smi_text(morph.mol, molecule_1, molecule_2)
        download_columns = st.columns(2)
        with download_columns[0]:
            st.download_button(
                "Download .smi file",
                data=smi_text,
                file_name="morph_pathway.smi",
                mime="chemical/x-daylight-smiles",
                key=f"morph_smi_{hash(smi_text) & 0xFFFFFFFF:X}",
            )

        try:
            cdxml_text = pathway_to_cdxml(morph.mol, molecule_1, molecule_2)
        except (RuntimeError, ValueError) as e:
            st.warning(f"Could not create ChemDraw file: {e}")
        else:
            with download_columns[1]:
                st.download_button(
                    "Download CDXML for ChemDraw",
                    data=cdxml_text,
                    file_name="morph_pathway.cdxml",
                    mime="chemical/x-cdxml",
                    key=f"morph_cdxml_{hash(cdxml_text) & 0xFFFFFFFF:X}",
                )


if __name__ == "__main__":
    main()
