class MultipleReplace:
    def __init__(self, a):
        self.root = ({}, {})
        for k, v in a:
            n, val = self.root
            for i in k[:-1]:
                if i not in n:
                    n[i] = ({}, {})
                n, val = n[i]
            val[k[-1]] = v

    def __repr__(self):
        return self.root.__repr__()

    def find(self, k, nk, i):
        n, val = self.root
        target_val = None
        target_depth = None
        depth = 1
        while i < nk:
            if k[i] in val:
                target_val = val[k[i]]
                target_depth = depth
            if k[i] in n:
                n, val = n[k[i]]
                depth += 1
            else:
                break
            i += 1
        return target_val, target_depth

    def replace(self, s):
        t = []
        i = 0
        n = len(s)
        while i < n:
            val, depth = self.find(s, n, i)
            if val is not None:
                t.append(val)
                i += depth
            else:
                t.append(s[i])
                i += 1
        return t
