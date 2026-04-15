from re import compile, escape

# 12 dihedral symbols: 30 degree discretization
DIHEDRALS = ["^^", "^\\", "<\\", ">>", "</", "_/", "__", "\\_", "\\>", "<<", "/>", "/^"]

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
    f"(?P<atom>{escape(L_BRACKET)}[0-9]*[A-Za-z][A-Za-z0-9]*{_c}{escape(R_BRACKET)}|"
    f"{escape(L_PAREN)}[0-9A-Za-z]+{escape(R_PAREN)}|[A-Za-z]{_c})"
)
_pring = f"(?P<ring>({escape(L_BRACKET)}[0-9]+{escape(R_BRACKET)}|[3-9]){escape(SKIP)}*)"
_psaturate = f"(?P<saturate>{escape(DOT)})"
_pmolsep = f"(?P<molsep>{escape(MOLSEP)})"
_pampersand = f"(?P<ampersand>{escape(AMPERSAND)})"

RegExp = compile(f"({_pbond}?({_patom}|({_pring})))|{_psaturate}|{_pmolsep}|{_pampersand}")

_dihedral_pat = "|".join(map(escape, sorted(DIHEDRALS, key=len, reverse=True)))
_implicit_c_re = compile(f"({_dihedral_pat})(?={_dihedral_pat})")
_dihedral_set = set(DIHEDRALS)


def _insert_implicit_carbon(s: str) -> str:
    return _implicit_c_re.sub(r"\1C", s)


def _remove_implicit_carbon(t: list[str]) -> list[str]:
    result: list[str] = []
    for i, tok in enumerate(t):
        if (
            tok == "C"
            and len(result) > 0
            and result[-1] in _dihedral_set
            and i + 1 < len(t)
            and t[i + 1] in _dihedral_set
        ):
            continue
        result.append(tok)
    return result


def ToTokens(s: str) -> list[str]:
    """Convert AMSR string to a list of tokens

    :param s: AMSR
    :return: list of tokens
    """
    t = []
    for m in RegExp.finditer(_insert_implicit_carbon(s)):
        g = m.groupdict()
        for k in ["bond", "atom", "saturate", "molsep"]:
            if g[k] is not None:
                t.append(g[k])
        if g["ring"] is not None:
            t.append(g["ring"].replace(SKIP, ""))
            t.extend(SKIP * g["ring"].count(SKIP))
    return t
