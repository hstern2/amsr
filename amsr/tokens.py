from re import compile, escape

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

_pbond = r"(?P<bond>[" + escape(E) + escape(Z) + "])"
_c = "[" + escape(PLUS + MINUS + RADICAL + EXTRA_PI + BANG + CW + CCW) + "]*"
_patom = r"(?P<atom>\[\d*[A-Za-z][a-z]?" + _c + r"\]|[A-Za-z]" + _c + r")"
_pring = r"(?P<ring>[3-6]+)"
_pskip = r"(?P<skip>" + escape(DOT) + r"*)"
_psaturate = r"(?P<saturate>" + escape(DOT) + r")"
_pdangling = r"(?P<dangling>\[\])"
_pampersand = r"(?P<ampersand>\&)"
_pnop = r"(?P<nop>\s+)"
_re = compile(
    f"({_pbond}?({_patom}|({_pring}{_pskip})))|{_psaturate}|{_pdangling}|{_pampersand}|{_pnop}"
)


def ToTokens(s):
    return [m.group() for m in Matches(s)]


def Matches(s):
    return _re.finditer(s)


def _isDot(t):
    return t == DOT


def _isRingDigit(t):
    return t in ("3", "4", "5", "6")


def _isDotOrRingDigit(t):
    return _isDot(t) or _isRingDigit(t)


def _needsNOP(prv, nxt):
    return (_isRingDigit(prv) and _isDotOrRingDigit(nxt)) or (
        _isDot(prv) and _isDot(nxt)
    )


def RemoveTrailingDots(tok):
    i = len(tok) - 1
    while i >= 0:
        if not _isDot(tok[i]):
            if not _isRingDigit(tok[i]):
                del tok[i + 1 :]
            break
        i -= 1
    return tok


def OmitUnneededNOPs(tok):
    ntok = len(tok)
    return [
        t
        for i, t in enumerate(tok)
        if t != NOP or (0 < i < ntok - 1 and _needsNOP(tok[i - 1], tok[i + 1]))
    ]
