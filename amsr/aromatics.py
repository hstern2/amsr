from itertools import combinations, permutations


def _letter_index(i: int) -> str:
    return chr(ord("a") + i)


def AddAromatics(g: dict[str, list[str]]) -> None:
    # five-membered rings
    for j in range(5):
        for c in ["O", "N", "S"]:
            ring = list("ccccc5")
            ring[j] = c
            idx = _letter_index(j)
            g[f"(5{idx}{c})"] = ["".join(ring)]

    for perm in permutations(range(5), 2):
        for c in ["O", "N", "S"]:
            ring = list("ccccc5")
            ring[perm[0]] = "n"
            ring[perm[1]] = c
            idx0 = _letter_index(perm[0])
            idx1 = _letter_index(perm[1])
            g[f"(5{idx0}N{idx1}{c})"] = ["".join(ring)]

    # six-membered rings
    for i in range(4):
        for comb in combinations(range(6), i):
            ring = list("cccccc6")
            idx = ""
            for k in comb:
                ring[k] = "n"
                idx += _letter_index(k)
            g[f"(6{idx}{'N' if idx else ''})"] = ["".join(ring)]
