from re import compile, escape
from typing import List

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
_re = compile(
    f"({_pbond}?({_patom}|({_pring})))|{_psaturate}|{_pdangling}|{_pampersand}|{_pnop}"
)


def ToTokens(s: str) -> List[str]:
    """Convert an AMSR string to a list of tokens

    :param s: AMSR
    :return: list of tokens
    """
    t = []
    for m in Matches(s):
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


def Matches(s):
    return _re.finditer(s)


def _isDot(t):
    return t == DOT


def _isRingDigit(t):
    return t in RING_DIGITS


def _isDotOrRingDigit(t):
    return _isDot(t) or _isRingDigit(t)


def _needsNOP(prv, nxt):
    return (_isRingDigit(prv) and _isDotOrRingDigit(nxt)) or (
        _isDot(prv) and _isDot(nxt)
    )


def RemoveTrailingDots(t):
    i = len(t) - 1
    while i >= 0:
        if not _isDot(t[i]):
            if not _isRingDigit(t[i]):
                del t[i + 1 :]
            break
        i -= 1
    return t


def OmitUnneededNOPs(t):
    n = len(t)
    while n > 0 and t[n - 1] == NOP:
        n -= 1
    return [
        t[i]
        for i in range(n)
        if t[i] != NOP or (0 < i < n - 1 and _needsNOP(t[i - 1], t[i + 1]))
    ]
