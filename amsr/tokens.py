from re import compile, escape, Match
from typing import Iterator, List

DOT = "."
DIHEDRALS = ["^", "^\\", "<\\", "<", "</", "_/", "_", "\\_", "\\>", ">", "/>", "/^"]
Z = DIHEDRALS[0]
E = DIHEDRALS[6]
CW = "'"
CCW = "`"
PLUS = "+"
MINUS = "-"
EXTRA_PI = ":"
BANG = "!"
RADICAL = "*"
L_BRACKET = "["
R_BRACKET = "]"
SKIP = "@"
MOLSEP = ";"

DIHEDRAL_FOR_BOND_SYMBOL = {
    s: 30 * (i - 12 if i > 6 else i) for i, s in enumerate(DIHEDRALS)
}
BOND_SYMBOL_FOR_DIHEDRAL = {v: k for k, v in DIHEDRAL_FOR_BOND_SYMBOL.items()}
BOND_SYMBOL_FOR_DIHEDRAL[-180] = E

_pbond = f"(?P<bond>{'|'.join(map(escape, sorted(DIHEDRALS, key=len, reverse=True)))})"
_c = f"[{''.join(map(escape, [PLUS,MINUS,RADICAL,EXTRA_PI,BANG,CW,CCW]))}]*"
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
