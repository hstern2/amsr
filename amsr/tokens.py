from re import compile, escape, Match
from typing import Iterator, List

DOT = "."
NOP = " "
E = "_"
Z = "^"
CW = ")"
CCW = "("
PLUS = "+"
MINUS = "-"
EXTRA_PI = ":"
BANG = "!"
RADICAL = "*"
RING_DIGITS = "3456"
L_BRACKET = "["
R_BRACKET = "]"
AMPERSAND = "&"

_pbond = f"(?P<bond>[{escape(E)}{escape(Z)}])"
_c = f"[{escape(PLUS+MINUS+RADICAL+EXTRA_PI+BANG+CW+CCW)}]*"
_patom = f"(?P<atom>{escape(L_BRACKET)}[^{escape(R_BRACKET)}]+{escape(R_BRACKET)}|[A-Za-z]{_c})"
_pring = f"(?P<ring>[{RING_DIGITS}]+" + escape(DOT) + "*)"
_psaturate = f"(?P<saturate>{escape(DOT)})"
_pdangling = f"(?P<dangling>{escape(L_BRACKET)}{escape(R_BRACKET)})"
_pampersand = f"(?P<ampersand>{escape(AMPERSAND)})"
_pnop = f"(?P<nop>{escape(NOP)}+)"

RegExp = compile(
    f"({_pbond}?({_patom}|({_pring})))|{_psaturate}|{_pdangling}|{_pampersand}|{_pnop}"
)


def ToTokens(s: str) -> List[str]:
    """Convert AMSR string to a list of tokens

    :param s: AMSR
    :return: list of tokens
    """
    t = []
    for m in RegExp.finditer(s):
        g = m.groupdict()
        for k in ["bond", "atom"]:
            if g[k] is not None:
                t.append(g[k])
        if g["ring"] is not None:
            t.extend(g["ring"])
        for k in ["saturate", "dangling", "ampersand", "nop"]:
            if g[k] is not None:
                t.append(g[k])
    return t
