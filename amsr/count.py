from typing import Any


def Count(d: dict, k: Any) -> None:
    if k in d:
        d[k] += 1
    else:
        d[k] = 1
