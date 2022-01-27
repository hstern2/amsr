from re import compile, escape
from .mreplace import MultipleReplace
from .tokens import ToTokens

# [Boc]
# [Cbz]
# [Ms]
# [Piv]
# [PMB]
# [PMP]
# [Tf]
# [TMS]
# [Ts]

Groups = {
    "[Bn]": ["Ccccccc6 ......"],
    "[Bz]": ["cocccccc6 ....."],
    "[Ph]": ["cccccc6 ....."],
    "[OEt]": ["OCC.."],
    "[OMe]": ["OC."],
    "[NHAc]": ["NcoC.."],
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

# must be called (again) if Groups is changed
def InitializeGroups():
    global _mr, _pattern
    _mr = MultipleReplace([(ToTokens(g), k) for k, v in Groups.items() for g in v])
    _pattern = compile("(" + "|".join([escape(k) for k in Groups.keys()]) + ")")


InitializeGroups()


def DecodeGroups(s):
    return _pattern.sub(lambda m: Groups[m.group(1)][0], s)


def EncodeGroups(s):
    return _mr.replace(s)
