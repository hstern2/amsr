from rdkit import Chem
from re import match, sub
from .valence import VALENCE, BANGS
from .parity import IsEvenParity
from .tokens import CW, CCW, PLUS, MINUS, RADICAL, EXTRA_PI, BANG


def GetSeenIndex(a):
    return a.GetIntProp("_seenIndex")


def SetSeenIndex(a, i):
    return a.SetIntProp("_seenIndex", i)


def IsSeen(a):
    return a.HasProp("_seenIndex")


class Atom:
    def __init__(self, sym):
        self.sym = sym
        self.bangs = sym.count(BANG)
        m = match(r"\[(\d+)", sym)
        self.isotope = None if m is None else int(m.group(1))
        self.chg = sym.count(PLUS) - sym.count(MINUS)
        self.nrad = sym.count(RADICAL)
        if CCW in sym:
            self.ct = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        elif CW in sym:
            self.ct = Chem.ChiralType.CHI_TETRAHEDRAL_CW
        else:
            self.ct = Chem.ChiralType.CHI_UNSPECIFIED
        self.atomSym = sub(r"[^A-Za-z]", "", sym)
        if self.atomSym[0].islower():
            self.atomSym = self.atomSym[0].upper() + self.atomSym[1:]
            self.maxPiBonds = 1
        else:
            self.maxPiBonds = 0
        self.maxPiBonds += 2 * sym.count(EXTRA_PI)
        self.nPiBonds = 0
        self.maxNeighbors = (
            VALENCE[(self.atomSym, self.chg, self.bangs)] - self.nrad - self.maxPiBonds
        )
        self.nNeighbors = 0
        self.isSaturated = False

    def canBond(self):
        return (not self.isSaturated) and self.nNeighbors < self.maxNeighbors

    def nAvailablePiBonds(self):
        return self.maxPiBonds - self.nPiBonds

    def asRDAtom(self):
        a = Chem.Atom(self.atomSym)
        a.SetFormalCharge(self.chg)
        a.SetNumRadicalElectrons(self.nrad)
        a.SetChiralTag(self.ct)
        if self.isotope:
            a.SetIsotope(self.isotope)
        return a

    def isCarbon(self):
        return self.atomSym == "C"

    def symWith(self, s):
        sym = self.sym
        if sym.startswith("["):
            return "[" + sym[1:-1] + s + "]"
        else:
            return sym + s

    def asToken(self, a, mol):
        ct = a.GetChiralTag()
        isEven = IsEvenParity([GetSeenIndex(b) for b in a.GetNeighbors()])
        if ct == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return self.symWith(CCW if isEven else CW)
        elif ct == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return self.symWith(CW if isEven else CCW)
        else:
            return self.sym

    @classmethod
    def fromRDAtom(cls, a):
        atomSym = a.GetSymbol()
        chg = a.GetFormalCharge()
        valence = a.GetTotalValence()
        nrad = a.GetNumRadicalElectrons()
        isotope = a.GetIsotope()
        if atomSym == "He" or nrad > 4:  # RDKit weirdness?
            nrad = 0
        bangs = BANGS.get((atomSym, chg, valence), 0)
        q, r = divmod(VALENCE[(atomSym, chg, bangs)] - nrad - a.GetTotalDegree(), 2)
        c = f"{PLUS*chg if chg > 0 else ''}{MINUS*(-chg) if chg < 0 else ''}{RADICAL*nrad}{BANG*bangs}{EXTRA_PI*q}"
        sym = (f"{isotope}" if isotope else "") + (atomSym.lower() if r else atomSym)
        return cls(f"[{sym}{c}]" if len(atomSym) == 2 or isotope else f"{sym}{c}")
