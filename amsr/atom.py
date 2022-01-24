from rdkit import Chem
from re import match, search, sub
from .valence import VALENCE, BANGS
from .parity import IsEvenParity


def visitedIndex(rda, atom):
    return atom[rda.GetIdx()].visitedIndex


class Atom:
    def __init__(self, sym):
        self.sym = sym
        self.bangs = sym.count("!")
        m = match(r"\[(\d+)", sym)
        self.isotope = None if m is None else int(m.group(1))
        m = search(r"[\+\-]\d?", sym)
        if m is None:
            self.chg = 0
        elif m.group(0) == "+":
            self.chg = 1
        elif m.group(0) == "-":
            self.chg = -1
        else:
            self.chg = int(m.group(0))
        self.nrad = sym.count("*")
        if "(" in sym:
            self.ct = Chem.ChiralType.CHI_TETRAHEDRAL_CCW
        elif ")" in sym:
            self.ct = Chem.ChiralType.CHI_TETRAHEDRAL_CW
        else:
            self.ct = Chem.ChiralType.CHI_UNSPECIFIED
        self.atomSym = sub(r"[^A-Za-z]", "", sym)
        if self.atomSym[0].islower():
            self.atomSym = self.atomSym[0].upper() + self.atomSym[1:]
            self.maxPiBonds = 1
        else:
            self.maxPiBonds = 0
        self.maxPiBonds += 2 * sym.count(":")
        self.nPiBonds = 0
        self.maxNeighbors = (
            VALENCE[(self.atomSym, self.chg, self.bangs)] - self.nrad - self.maxPiBonds
        )
        self.nNeighbors = 0
        self.visitedIndex = None
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

    def symWith(self, s):
        sym = self.sym
        if sym.startswith("["):
            return "[" + sym[1:-1] + s + "]"
        else:
            return sym + s

    def asToken(self, a, atom):
        ct = a.GetChiralTag()
        isEven = IsEvenParity([visitedIndex(b, atom) for b in a.GetNeighbors()])
        if ct == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return self.symWith("(" if isEven else ")")
        elif ct == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return self.symWith(")" if isEven else "(")
        else:
            return self.sym

    @classmethod
    def fromRD(cls, a):
        atomSym = a.GetSymbol()
        chg = a.GetFormalCharge()
        valence = a.GetTotalValence()
        nrad = a.GetNumRadicalElectrons()
        isotope = a.GetIsotope()
        if atomSym == "He" or nrad > 4:
            nrad = 0
        bangs = BANGS.get((atomSym, chg, valence), 0)
        q, r = divmod(VALENCE[(atomSym, chg, bangs)] - nrad - a.GetTotalDegree(), 2)
        if chg == 1:
            c = "+"
        elif chg > 1:
            c = f"+{chg}"
        elif chg == -1:
            c = "-"
        elif chg < -1:
            c = f"-{-chg}"
        else:
            c = ""
        c += "!" * bangs + ":" * q + "*" * nrad
        sym = (f"{isotope}" if isotope else "") + (atomSym.lower() if r else atomSym)
        if len(atomSym) == 2 or chg or isotope:
            sym = "[" + sym + c + "]"
        else:
            sym += c
        return cls(sym)
