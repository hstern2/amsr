from itertools import combinations
from .count import Count


def BridgeheadAtoms(mol, ringInfo, n):
    for i, j in combinations([r for r in ringInfo.BondRings() if len(r) <= n], 2):
        bridge = set(i) & set(j)
        if len(bridge) >= 2:
            a = {}
            for b in bridge:
                b = mol.GetBondWithIdx(b)
                Count(a, b.GetBeginAtomIdx())
                Count(a, b.GetEndAtomIdx())
            yield from (k for k, n in a.items() if n == 1)


def RingAtoms(ringInfo, n):
    yield from (i for r in ringInfo.AtomRings() if len(r) == n for i in r)
