from re import compile, escape, Match
from typing import Iterator, List

DOT = "."
E = "_"
Z = "^"
CW = ")"
CCW = "("
PLUS = "+"
MINUS = "-"
EXTRA_PI = ":"
BANG = "!"
RADICAL = "*"
L_BRACKET = "["
R_BRACKET = "]"
PERCENT = "%"
SKIP = ";"

_pbond = f"(?P<bond>[{escape(E)}{escape(Z)}])"
_c = f"[{escape(PLUS+MINUS+RADICAL+EXTRA_PI+BANG+CW+CCW)}]*"
_patom = (
    f"(?P<atom>{escape(L_BRACKET)}[0-9]*[A-Za-z]+{_c}{escape(R_BRACKET)}|[A-Za-z]{_c})"
)
_pring = (
    f"(?P<ring>({escape(L_BRACKET)}[0-9]+{escape(R_BRACKET)}|[3-9]){escape(SKIP)}*)"
)
_psaturate = f"(?P<saturate>{escape(DOT)})"

RegExp = compile(f"({_pbond}?({_patom}|({_pring})))|{_psaturate}")


def ToTokens(s: str) -> List[str]:
    """Convert AMSR string to a list of tokens

    :param s: AMSR
    :return: list of tokens
    """
    t = []
    for m in RegExp.finditer(s):
        g = m.groupdict()
        for k in ["bond", "atom", "saturate"]:
            if g[k] is not None:
                t.append(g[k])
        if g["ring"] is not None:
            t.append(g["ring"].replace(SKIP, ""))
            t.extend(SKIP * g["ring"].count(SKIP))
    return t
