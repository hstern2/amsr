from re import compile, escape, Match
from typing import Iterator, List

# tmp

TMPTMP = 'aa'

DOT = "."
E = "_"
Z = "^"
GAUCHE_CW = ">"
GAUCHE_CCW = "<"
SKEW_CW = "\\"
SKEW_CCW = "/"
CW = ")"
CCW = "("
PLUS = "+"
MINUS = "-"
EXTRA_PI = ":"
BANG = "!"
RADICAL = "*"
L_BRACKET = "["
R_BRACKET = "]"
SKIP = "'"
MOLSEP = ";"

DIHEDRAL_FOR_BOND_SYMBOL = {
    Z: 0,
    GAUCHE_CW: -60,
    SKEW_CW: -120,
    E: 180,
    SKEW_CCW: 120,
    GAUCHE_CCW: 60,
}

BOND_SYMBOL_FOR_DIHEDRAL = {
    0: Z,
    -60: GAUCHE_CW,
    -120: SKEW_CW,
    180: E,
    -180: E,
    120: SKEW_CCW,
    60: GAUCHE_CCW,
}

_pbond = f"(?P<bond>[{escape(E)}{escape(Z)}{escape(GAUCHE_CW)}{escape(GAUCHE_CCW)}{escape(SKEW_CW)}{escape(SKEW_CCW)}])"
_c = f"[{escape(PLUS+MINUS+RADICAL+EXTRA_PI+BANG+CW+CCW)}]*"
_patom = (
    f"(?P<atom>{escape(L_BRACKET)}[0-9]*[A-Za-z]+{_c}{escape(R_BRACKET)}|[A-Za-z]{_c})"
)
_pring = (
    f"(?P<ring>({escape(L_BRACKET)}[0-9]+{escape(R_BRACKET)}|[3-9]){escape(SKIP)}*)"
)
_psaturate = f"(?P<saturate>{escape(DOT)})"
_pmolsep = f"(?P<molsep>{escape(MOLSEP)})"

RegExp = compile(f"({_pbond}?({_patom}|({_pring})))|{_psaturate}|{_pmolsep}")


def ToTokens(s: str) -> List[str]:
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
