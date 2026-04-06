from re import compile, escape

# 12 dihedral symbols: 30 degree discretization
DIHEDRALS = ["^", "^\\", "<\\", "<", "</", "_/", "_", "\\_", "\\>", ">", "/>", "/^"]

# 24 dihedral symbols: 15 degree discretization
# DIHEDRALS = [f"/{chr(x)}" for x in range(ord("a"), ord("y"))]

N_DIHEDRALS = len(DIHEDRALS)
Z = DIHEDRALS[0]
E = DIHEDRALS[N_DIHEDRALS // 2]
CW = "'"
CCW = "`"
PLUS = "+"
MINUS = "-"
EXTRA_PI = ":"
BANG = "!"
RADICAL = "*"
L_BRACKET = "["
R_BRACKET = "]"
L_PAREN = "("
R_PAREN = ")"
SKIP = "@"
MOLSEP = ";"
AMPERSAND = "&"
DOT = "."

DIHEDRAL_FOR_BOND_SYMBOL = {
    s: (360 // N_DIHEDRALS) * (i - N_DIHEDRALS if i > N_DIHEDRALS // 2 else i)
    for i, s in enumerate(DIHEDRALS)
}
BOND_SYMBOL_FOR_DIHEDRAL = {v: k for k, v in DIHEDRAL_FOR_BOND_SYMBOL.items()}
BOND_SYMBOL_FOR_DIHEDRAL[-180] = E

_pbond = f"(?P<bond>{'|'.join(map(escape, sorted(DIHEDRALS, key=len, reverse=True)))})"
_c = f"[{''.join(map(escape, [PLUS,MINUS,RADICAL,EXTRA_PI,BANG,CW,CCW]))}]*"
_patom = (
    f"(?P<atom>{escape(L_BRACKET)}[0-9]*[A-Za-z]+{_c}{escape(R_BRACKET)}|"
    f"{escape(L_PAREN)}[0-9A-Za-z]+{escape(R_PAREN)}|[A-Za-z]{_c})"
)
_pring = f"(?P<ring>({escape(L_BRACKET)}[0-9]+{escape(R_BRACKET)}|[3-9]){escape(SKIP)}*)"
_psaturate = f"(?P<saturate>{escape(DOT)})"
_pmolsep = f"(?P<molsep>{escape(MOLSEP)})"
_pampersand = f"(?P<ampersand>{escape(AMPERSAND)})"

RegExp = compile(f"({_pbond}?({_patom}|({_pring})))|{_psaturate}|{_pmolsep}|{_pampersand}")


def ToTokens(s: str) -> list[str]:
    """Convert AMSR string to a list of tokens

    :param s: AMSR
    :return: list of tokens
    """
    t = []
    for m in RegExp.finditer(s):
        g = m.groupdict()
        for k in ["bond", "atom", "saturate", "molsep"]:
            if g[k] is not None:
                t.append(g[k])
        if g["ring"] is not None:
            t.append(g["ring"].replace(SKIP, ""))
            t.extend(SKIP * g["ring"].count(SKIP))
    return t
