from re import compile, escape, Pattern
from .mreplace import MultipleReplace
from .tokens import ToTokens
from typing import Dict, List, Optional

_groups: Dict[str, List[str]] = {
    "[@]": ["cccccc6"],
    "[Cy]": ["CCCCCC6......"],
    "[Pip]": ["NCCCCC6....."],
    "[Ts]": ["S!!:oocccccc6..C...", "S!!:occcccc6..C...o", "S!!:cccccc6..C...oo"],
    "[Tf]": ["S!!:ooCFFF", "S!!:oCFFFo", "S!!:CFFFoo"],
    "[Ms]": ["S!!:ooC.", "S!!:oC.o", "S!!:C.oo"],
    "[Piv]": ["coCC.C.C.", "cCC.C.C.o"],
    "[Boc]": ["coOCC.C.C.", "cOCC.C.C.o"],
    "[Tol]": ["cccccc6..C..."],
    "[Cbz]": ["coOCcccccc6......", "cOCcccccc6......o"],
    "[Bn]": ["Ccccccc6......"],
    "[Bz]": ["cocccccc6....."],
    "[Ph]": ["cccccc6....."],
    "[OEt]": ["OCC.."],
    "[OMe]": ["OC."],
    "[SMe]": ["SC."],
    "[NHAc]": ["NcoC..", "NcC.o."],
    "[NHMe]": ["NC.."],
    "[NMe2]": ["NC.C."],
    "[OAc]": ["OcoC.", "OcC.o"],
    "[COOEt]": ["coOCC..", "cOoCC.."],
    "[COOMe]": ["coOC.", "cOoC."],
    "[COO-]": ["coO-", "cO-o"],
    "[NO2]": ["n+oO-", "n+O-o"],
    "[Ac]": ["coC.", "cC.o"],
    "[COOH]": ["coO.", "cO.o"],
    "[CHO]": ["co."],
    "[Et]": ["CC.."],
    "[nPr]": ["CCC..."],
    "[iPr]": ["CC.C.."],
    "[tBu]": ["CC.C.C."],
    "[nBu]": ["CCCC...."],
    "[sBu]": ["CC.CC...", "CCC..C.."],
    "[iBu]": ["CCC.C..."],
    "[OiBu]": ["OCCC.C..."],
    "[nPent]": ["CCCCC....."],
    "[nHex]": ["CCCCCC......"],
    "[nHept]": ["CCCCCCC......."],
    "[nOct]": ["CCCCCCCC........"],
    "[nNon]": ["CCCCCCCCC........."],
    "[nDec]": ["CCCCCCCCCC.........."],
    "[CN]": ["C:N:"],
    "[NC]": ["N+:C:-"],
    "[CF3]": ["CFFF"],
    "[CCl3]": ["C[Cl][Cl][Cl]"],
    "[PO3]": ["p!oO.O.", "p!O.oO.", "p!O.O.o"],
    "[SO3]": ["S!!:O.oo", "S!!:oO.o", "S!!:ooO."],
}


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
