from itertools import combinations, permutations
from typing import Dict, List


def _letter_index(i):
    return chr(ord("a") + i)


def AddAromatics(g: Dict[str, List[str]]) -> None:

    # five-membered rings
    for j in range(5):
        for c in ["O", "N", "S"]:
            ring = list("ccccc5")
            ring[j] = c
            idx = _letter_index(j)
            g[f"(5{idx}{c})"] = ["".join(ring)]

    for j in permutations(range(5), 2):
        for c in ["O", "N", "S"]:
            ring = list("ccccc5")
            ring[j[0]] = "n"
            ring[j[1]] = c
            idx0 = _letter_index(j[0])
            idx1 = _letter_index(j[1])
            g[f"(5{idx0}n{idx1}{c})"] = ["".join(ring)]

    # six-membered rings
    for i in range(4):
        for j in combinations(range(6), i):
            ring = list("cccccc6")
            idx = ""
            for k in j:
                ring[k] = "n"
                idx += _letter_index(k)
            g[f"(6{idx}{'n' if idx else ''})"] = ["".join(ring)]
