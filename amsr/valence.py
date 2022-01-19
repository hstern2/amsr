# key: symbol, charge, number of !'s (hypervalent marker)
VALENCE = {
    ("H", 0, 0): 1,
    ("H", 1, 0): 0,
    ("He", 0, 0): 0,
    ("Li", 1, 0): 0,
    ("Li", 0, 0): 1,
    ("Be", 2, 0): 0,
    ("B", 0, 0): 3,
    ("B", -1, 0): 4,
    ("C", 0, 0): 4,
    ("C", -1, 0): 3,
    ("C", 1, 0): 3,
    ("N", 0, 0): 3,
    ("N", 1, 0): 4,
    ("N", -1, 0): 2,
    ("O", 0, 0): 2,
    ("O", -1, 0): 1,
    ("O", -2, 0): 0,
    ("O", 1, 0): 3,
    ("F", 0, 0): 1,
    ("F", -1, 0): 0,
    ("F", 1, 0): 2,
    ("Na", 0, 0): 1,
    ("Na", 1, 0): 0,
    ("Mg", 0, 0): 2,
    ("Mg", 1, 0): 1,
    ("Mg", 2, 0): 0,
    ("Al", 0, 0): 3,
    ("Al", 3, 0): 0,
    ("Al", -3, 0): 6,
    ("Xe", 0, 0): 0,
    ("Si", 0, 0): 4,
    ("Si", -1, 0): 5,
    ("P", 0, 0): 3,
    ("P", 0, 1): 5,
    ("P", 0, 2): 7,
    ("P", 1, 0): 4,
    ("P", -1, 1): 6,
    ("S", 0, 0): 2,
    ("S", 0, 1): 4,
    ("S", 0, 2): 6,
    ("S", 1, 0): 3,
    ("S", 1, 1): 5,
    ("S", -1, 0): 1,
    ("S", -1, 1): 3,
    ("S", -2, 0): 0,
    ("Cl", 0, 0): 1,
    ("Cl", -1, 0): 0,
    ("Cl", 1, 0): 2,
    ("Cl", 2, 0): 3,
    ("Cl", 3, 0): 4,
    ("K", 0, 0): 1,
    ("K", 1, 0): 0,
    ("Ca", 0, 0): 2,
    ("Ca", 2, 0): 0,
    ("Zn", 0, 0): 2,
    ("Zn", 1, 0): 1,
    ("Zn", 2, 0): 1,
    ("Zn", -2, 0): 4,
    ("As", 0, 0): 3,
    ("As", 0, 1): 5,
    ("As", 0, 2): 7,
    ("As", 1, 0): 4,
    ("As", -1, 1): 6,
    ("Se", 0, 0): 2,
    ("Se", 0, 1): 4,
    ("Se", 0, 2): 6,
    ("Se", 1, 0): 3,
    ("Se", 1, 1): 5,
    ("Se", -1, 0): 1,
    ("Se", -2, 0): 0,
    ("Cs", 1, 0): 0,
    ("Cs", 0, 0): 1,
    ("Ba", 0, 0): 2,
    ("Ba", 2, 0): 0,
    ("Bi", 0, 0): 3,
    ("Bi", 3, 0): 0,
    ("Br", 0, 0): 1,
    ("Br", -1, 0): 0,
    ("Br", 2, 0): 3,
    ("Kr", 0, 0): 0,
    ("Rb", 0, 0): 1,
    ("Rb", 1, 0): 0,
    ("Sr", 0, 0): 2,
    ("Sr", 2, 0): 0,
    ("Ag", 0, 0): 1,
    ("Ag", 1, 0): 0,
    ("Ag", -4, 0): 2,
    ("Te", 0, 0): 2,
    ("Te", 1, 0): 3,
    ("Te", 0, 1): 4,
    ("Te", 0, 2): 6,
    ("Te", -1, 1): 3,
    ("Te", -1, 2): 5,
    ("I", 0, 0): 1,
    ("I", 0, 1): 3,
    ("I", 0, 2): 5,
    ("I", -1, 0): 0,
    ("I", 1, 0): 2,
    ("I", 2, 1): 3,
    ("I", 3, 0): 4,
    ("Ra", 0, 0): 2,
    ("Ra", 2, 0): 0,
}

# key: symbol, charge, valence
BANGS = {
    ("S", 0, 4): 1,
    ("S", 1, 5): 1,
    ("S", -1, 3): 1,
    ("S", 0, 6): 2,
    ("P", 0, 5): 1,
    ("P", 0, 7): 2,
    ("P", -1, 6): 1,
    ("As", 0, 5): 1,
    ("As", 0, 7): 2,
    ("As", -1, 6): 1,
    ("Se", 0, 4): 1,
    ("Se", 1, 5): 1,
    ("Se", -1, 3): 1,
    ("Se", 0, 6): 2,
    ("Te", 0, 4): 1,
    ("Te", 0, 6): 2,
    ("Te", -1, 3): 1,
    ("Te", -1, 5): 2,
    ("I", 0, 3): 1,
    ("I", 0, 5): 2,
    ("I", 2, 3): 1,
}
