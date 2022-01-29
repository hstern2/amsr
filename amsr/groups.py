from re import compile, escape
from .mreplace import MultipleReplace
from .tokens import ToTokens
from typing import Dict, List

_groups: Dict[str, List[str]] = {
    "[Ts]": ["S!!:oocccccc6 ..C...", "S!!:occcccc6 ..C...o", "S!!:cccccc6 ..C...oo"],
    "[Tf]": ["S!!:ooCFFF", "S!!:oCFFFo", "S!!:CFFFoo"],
    "[Ms]": ["S!!:ooC.", "S!!:oC.o", "S!!:C.oo"],
    "[Piv]": ["coCC.C.C.", "cCC.C.C.o"],
    "[Boc]": ["coOCC.C.C.", "cOCC.C.C.o"],
    "[Tol]": ["cccccc6 ..C..."],
    "[Cbz]": ["coOCcccccc6 ......", "cOCcccccc6 ......o"],
    "[Bn]": ["Ccccccc6 ......"],
    "[Bz]": ["cocccccc6 ....."],
    "[Ph]": ["cccccc6 ....."],
    "[OEt]": ["OCC.."],
    "[OMe]": ["OC."],
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
    "[tBu]": ["CC.C.C."],
    "[nBu]": ["CCCC...."],
    "[sBu]": ["CC.CC...", "CCC..C.."],
    "[iBu]": ["CCC.C..."],
    "[nPr]": ["CCC..."],
    "[iPr]": ["CC.C."],
    "[Et]": ["CC.."],
    "[CN]": ["C:N:"],
    "[CF3]": ["CFFF"],
    "[CCl3]": ["C[Cl][Cl][Cl]"],
}


def Groups() -> Dict[str, List[str]]:
    """Keys are functional group abbreviations, values are lists of one or more AMSR strings.
    May be modified, but :func InitializeGroups: must be called after modification.

    :return: Groups dictionary
    """


def InitializeGroups() -> None:
    """Initialize tree and compile regular expression for converting between
    group abbreviations and tokens.  Must be called after modification of
    :func Groups: dictionary.
    """
    global _mr, _pattern
    _mr = MultipleReplace([(ToTokens(g), k) for k, v in _groups.items() for g in v])
    _pattern = compile("(" + "|".join([escape(k) for k in _groups.keys()]) + ")")


InitializeGroups()


def DecodeGroups(s):
    return _pattern.sub(lambda m: _groups[m.group(1)][0], s)


def EncodeGroups(s):
    return _mr.replace(s)
