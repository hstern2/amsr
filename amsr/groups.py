from re import sub

GROUPS = {
    "[Ph]": list("cccccc6 ....."),
    "[Ac]": list("coC."),
    "[COOH]": list("coO."),
    "[COO-]": list("co") + ["[O-]"],
    "[tBu]": list("CC.C.C."),
    "[nBu]": list("CCCC...."),
    "[iBu]": list("CCC.C.."),
    "[nPr]": list("CCC..."),
    "[iPr]": list("CC.C."),
    "[Et]": list("CC.."),
}

NSYM = {}
for k, v in GROUPS.items():
    j = len(v) - 1
    while j > 0 and v[j] in (".", " "):
        j -= 1
    NSYM[k] = j + 1

GROUPSTR = {k: "".join(v) for k, v in GROUPS.items()}


def DecodeGroups(s):
    for k, v in GROUPSTR.items():
        s = s.replace(k, v)
    return s


def EncodeGroups(tok):
    ntok = len(tok)
    i = 0
    while i < ntok:
        for k, v in GROUPS.items():
            n = len(v)
            j = max(NSYM[k], min(n, ntok - i))
            if tok[i : i + j] == v[:j]:
                yield k
                i += n
                break
        else:
            yield tok[i]
            i += 1
