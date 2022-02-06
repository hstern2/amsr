import anytree
import rdkit.Chem as rd
from .atom import GetSeenIndex


def _idx(n):
    return n.name.GetIdx()


def BFSPath(n):
    if n is not None:
        yield _idx(n)
        yield from BFSPath(n.parent)


def BFSTree(a, n):
    ai = a.GetIdx()
    if n == 0:
        return [ai]
    root = anytree.Node(a)
    q = [root]
    seen = {ai}
    tree = []
    while q:
        node = q.pop(0)
        a = node.name
        for b in a.GetBonds():
            c = b.GetOtherAtom(a)
            ci = c.GetIdx()
            if ci not in seen:
                seen.add(ci)
                child = anytree.Node(c, parent=node)
                if child.depth < n:
                    q.append(child)
                elif ci < ai:
                    tree.append(child)
    return sorted(tree, key=_idx, reverse=True)


def BFSFind(a, targetIdx, seenBonds):
    ai = a.GetIdx()
    root = anytree.Node(a)
    q = [root]
    seen = {ai}
    targetNode = None
    while q:
        node = q.pop(0)
        a = node.name
        ai = a.GetIdx()
        for b in a.GetBonds():
            c = b.GetOtherAtom(a)
            ci = c.GetIdx()
            if frozenset([ai, ci]) in seenBonds and ci not in seen:
                seen.add(ci)
                child = anytree.Node(c, parent=node)
                if ci == targetIdx:
                    targetNode = child
                if targetNode is None:
                    q.append(child)
    return (
        (n.name.GetIdx(), n.depth)
        for n in sorted(
            anytree.LevelOrderIter(root, filter_=lambda n: n.depth == targetNode.depth),
            key=lambda n: GetSeenIndex(n.name),
            reverse=True,
        )
    )
