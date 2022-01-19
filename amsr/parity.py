def IsEvenParity(a):
    s = 0
    for i, ai in enumerate(a):
        for j, aj in enumerate(a):
            if i < j and ai > aj:
                s += 1
    return s % 2 == 0
