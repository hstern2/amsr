from re import compile, escape

GROUPS = {
    "[Ph]": [list("cccccc6 .....")],
    "[Ac]": [list("coC."), list("cC.o")],
    "[COOH]": [list("coO."), list("cO.o")],
    "[COO-]": [["c", "o", "[O-]"], ["c", "[O-]", "o"]],
    "[nDec]": [list("CCCCCCCCCC..........")],
    "[nNon]": [list("CCCCCCCCC.........")],
    "[nOct]": [list("CCCCCCCC........")],
    "[nHept]": [list("CCCCCCC.......")],
    "[nHex]": [list("CCCCCC......")],
    "[nPent]": [list("CCCCC.....")],
    "[tBu]": [list("CC.C.C.")],
    "[nBu]": [list("CCCC....")],
    "[iBu]": [list("CCC.C..")],
    "[nPr]": [list("CCC...")],
    "[iPr]": [list("CC.C.")],
    "[Et]": [list("CC..")],
}

groupStr = {k: "".join(v[0]) for k, v in GROUPS.items()}
pattern = compile("(" + "|".join([escape(k) for k in groupStr.keys()]) + ")")


def DecodeGroups(s):
    return pattern.sub(lambda m: groupStr[m.group(1)], s)


def EncodeGroups(tok):
    ntok = len(tok)
    i = 0
    t = []
    while i < ntok:
        for k, v in GROUPS.items():
            for g in v:
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
