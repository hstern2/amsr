from re import compile, escape

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
    "[Bn]": [list("Ccccccc6 ......")],
    "[Bz]": [list("cocccccc6 .....")],
    "[Ph]": [list("cccccc6 .....")],
    "[OEt]": [list("OCC..")],
    "[OMe]": [list("OC.")],
    "[NHAc]": [list("NcoC..")],
    "[NHMe]": [list("NC..")],
    "[NMe2]": [list("NC.C.")],
    "[OAc]": [list("OcoC."), list("OcC.o")],
    "[COOEt]": [list("coOCC.."), list("cOoCC..")],
    "[COOMe]": [list("coOC."), list("cOoC.")],
    "[COO-]": [["c", "o", "O-"], ["c", "O-", "o"]],
    "[NO2]": [["n+", "o", "O-"], ["n+", "O-", "o"]],
    "[Ac]": [list("coC."), list("cC.o")],
    "[COOH]": [list("coO."), list("cO.o")],
    "[CHO]": [list("co.")],
    "[tBu]": [list("CC.C.C.")],
    "[nBu]": [list("CCCC....")],
    "[sBu]": [list("CC.CC..."), list("CCC..C..")],
    "[iBu]": [list("CCC.C...")],
    "[nPr]": [list("CCC...")],
    "[iPr]": [list("CC.C.")],
    "[Et]": [list("CC..")],
    "[CN]": [["C:", "N:"]],
    "[CF3]": [list("CFFF")],
    "[CCl3]": [["C"] + ["[Cl]"] * 3],
}

# must be called (again) if Groups is changed
def InitializeGroups():
    global _groupSym, _groupStr, _pattern
    _groupSym = sorted(Groups.keys(), key=lambda k: len(Groups[k][0]), reverse=True)
    _groupStr = {k: "".join(v[0]) for k, v in Groups.items()}
    _pattern = compile("(" + "|".join([escape(k) for k in _groupStr.keys()]) + ")")


InitializeGroups()


def DecodeGroups(s):
    return _pattern.sub(lambda m: _groupStr[m.group(1)], s)


def EncodeGroups(tok):
    ntok = len(tok)
    i = 0
    t = []
    while i < ntok:
        for k in _groupSym:
            for g in Groups[k]:
                n = len(g)
                j = min(n, ntok - i)
                if tok[i : i + j] == g[:j]:
                    t.append(k)
                    i += n
                    break
            else:
                continue
            break
        else:
            t.append(tok[i])
            i += 1
    return t
