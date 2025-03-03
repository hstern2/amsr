from typing import Dict, List, Optional
from re import compile, escape, Pattern
from .mreplace import MultipleReplace
from .tokens import ToTokens
import json
from importlib.resources import files


def from_json(fname: str):
    with files("amsr.data").joinpath(fname).open() as f:
        return json.load(f)


_groups: Dict[str, List[str]] = from_json("groups.json") | from_json(
    "aromatic_rings.json"
)


def Groups() -> Dict[str, List[str]]:
    """Keys are functional group abbreviations, values are lists of one or more AMSR strings
    consisting only of atom/bond tokens.  May be modified, but :func:`InitializeGroups`
    must be called after modification.

    :return: Groups dictionary
    """
    return _groups


_mr: Optional[MultipleReplace] = None
_pattern: Optional[Pattern] = None


def InitializeGroups() -> None:
    """Initialize tree and compile regular expression for converting between
    group abbreviations and tokens.  Must be called after modification of
    :func:`Groups` dictionary.
    """
    global _pattern, _mr
    _mr = MultipleReplace([(ToTokens(g), k) for k, v in _groups.items() for g in v])
    _pattern = compile("(" + "|".join([escape(k) for k in _groups.keys()]) + ")")


InitializeGroups()


def DecodeGroups(s: str) -> str:
    return s if _pattern is None else _pattern.sub(lambda m: _groups[m.group(1)][0], s)


def EncodeGroups(s: List[str]) -> List[str]:
    return s if _mr is None else _mr.replace(s)
