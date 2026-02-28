from typing import Optional

import numpy as np
from rdkit import Chem

# Ideal bond lengths in Angstroms: (element1, element2, bond_order) -> length
# Elements are alphabetically sorted
_BOND_LENGTHS = {
    ("Br", "C", 1): 1.94,
    ("C", "C", 1): 1.54,
    ("C", "C", 1.5): 1.40,
    ("C", "C", 2): 1.34,
    ("C", "C", 3): 1.20,
    ("C", "Cl", 1): 1.77,
    ("C", "F", 1): 1.35,
    ("C", "I", 1): 2.14,
    ("C", "N", 1): 1.47,
    ("C", "N", 1.5): 1.34,
    ("C", "N", 2): 1.29,
    ("C", "N", 3): 1.16,
    ("C", "O", 1): 1.43,
    ("C", "O", 1.5): 1.33,
    ("C", "O", 2): 1.23,
    ("C", "S", 1): 1.81,
    ("C", "S", 1.5): 1.73,
    ("C", "S", 2): 1.60,
    ("N", "N", 1): 1.45,
    ("N", "N", 1.5): 1.35,
    ("N", "N", 2): 1.25,
    ("N", "N", 3): 1.10,
    ("N", "O", 1): 1.40,
    ("N", "O", 1.5): 1.28,
    ("N", "O", 2): 1.21,
    ("N", "S", 1): 1.65,
    ("O", "O", 1): 1.48,
    ("O", "P", 1): 1.63,
    ("O", "P", 2): 1.48,
    ("O", "S", 2): 1.43,
    ("N", "P", 1): 1.68,
}

# Covalent radii fallback (Angstroms)
_COVALENT_RADII = {
    "H": 0.31,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Br": 1.20,
    "I": 1.39,
}

# Bond angles by hybridization (degrees)
_BOND_ANGLES = {
    Chem.HybridizationType.SP3: 109.5,
    Chem.HybridizationType.SP2: 120.0,
    Chem.HybridizationType.SP: 180.0,
}


def _get_bond_length(mol, i, j):
    bond = mol.GetBondBetweenAtoms(i, j)
    sym_i = mol.GetAtomWithIdx(i).GetSymbol()
    sym_j = mol.GetAtomWithIdx(j).GetSymbol()
    if bond.GetIsAromatic():
        bo = 1.5
    elif bond.GetBondType() == Chem.BondType.DOUBLE:
        bo = 2
    elif bond.GetBondType() == Chem.BondType.TRIPLE:
        bo = 3
    else:
        bo = 1
    key = (min(sym_i, sym_j), max(sym_i, sym_j), bo)
    if key in _BOND_LENGTHS:
        return _BOND_LENGTHS[key]
    r_i = _COVALENT_RADII.get(sym_i, 1.5)
    r_j = _COVALENT_RADII.get(sym_j, 1.5)
    return r_i + r_j


def _get_bond_angle(mol, g, p, i):
    """Get bond angle g-p-i, using ring polygon angle for SP2 atoms."""
    hyb = mol.GetAtomWithIdx(p).GetHybridization()
    if hyb == Chem.HybridizationType.SP2:
        ri = mol.GetRingInfo()
        bond_gp = mol.GetBondBetweenAtoms(g, p)
        bond_pi = mol.GetBondBetweenAtoms(p, i)
        if bond_gp is not None and bond_pi is not None:
            rings_gp = set(ri.BondRingSizes(bond_gp.GetIdx()))
            rings_pi = set(ri.BondRingSizes(bond_pi.GetIdx()))
            common = rings_gp & rings_pi
            if common:
                n = min(common)
                return (n - 2) * 180.0 / n
    return _BOND_ANGLES.get(hyb, 109.5)


def _place_atom(A, B, C, d, theta_deg, omega_deg):
    """Place atom D given reference points A, B, C.

    d: bond length C-D
    theta_deg: bond angle B-C-D in degrees
    omega_deg: torsion angle A-B-C-D in degrees
    """
    theta = np.radians(theta_deg)
    omega = np.radians(omega_deg)
    BC = C - B
    bc = np.linalg.norm(BC)
    if bc < 1e-10:
        BC = np.array([1.0, 0.0, 0.0])
    else:
        BC = BC / bc
    AB = B - A
    n = np.cross(AB, BC)
    nn = np.linalg.norm(n)
    if nn < 1e-10:
        if abs(BC[0]) < 0.9:
            n = np.cross(BC, np.array([1.0, 0.0, 0.0]))
        else:
            n = np.cross(BC, np.array([0.0, 1.0, 0.0]))
        n = n / np.linalg.norm(n)
    else:
        n = n / nn
    m = np.cross(n, BC)
    return C + d * (
        -np.cos(theta) * BC + np.sin(theta) * np.cos(omega) * m + np.sin(theta) * np.sin(omega) * n
    )


def _synthetic_gg(coords, g, p):
    """Create a synthetic great-grandparent reference point."""
    gp = coords[p] - coords[g]
    gp_norm = np.linalg.norm(gp)
    if gp_norm < 1e-10:
        gp = np.array([1.0, 0.0, 0.0])
    else:
        gp = gp / gp_norm
    if abs(gp[0]) < 0.9:
        perp = np.cross(gp, np.array([1.0, 0.0, 0.0]))
    else:
        perp = np.cross(gp, np.array([0.0, 1.0, 0.0]))
    perp = perp / np.linalg.norm(perp)
    return coords[g] - perp


def _default_torsion(mol, p, child_count, in_ring, gg_in_ring=True):
    """Default torsion for placing a new child of atom p.

    gg_in_ring: whether the great-grandparent is inside the ring.
    When False (external-to-ring transition), sp2 ring torsion is flipped
    so the ring continues away from the external atom.
    """
    hyb = mol.GetAtomWithIdx(p).GetHybridization()
    if hyb == Chem.HybridizationType.SP3:
        return 180.0 + 120.0 * child_count
    elif hyb == Chem.HybridizationType.SP2:
        if in_ring:
            if gg_in_ring:
                return 0.0 if child_count == 0 else 180.0
            else:
                return 180.0 if child_count == 0 else 0.0
        else:
            return 180.0 if child_count == 0 else 0.0
    else:
        return 180.0


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
) -> Chem.Mol:
    """Z-matrix conformer from AMSR dihedrals and ideal geometry.

    :param mol: RDKit Mol (from ToMol)
    :param dihedral: dictionary of dihedral angle constraints {(mi,i,j,mj): angle}
    :return: RDKit Mol with 3D conformer
    """
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 0:
        return mol

    coords = np.zeros((n_atoms, 3))

    # Build dihedral lookup: (g, p) -> (mi, mj, angle)
    # so we can look up AMSR dihedral by the g-p bond axis
    bond_dihedral = {}
    if dihedral:
        for (mi, i, j, mj), angle in dihedral.items():
            bond_dihedral[(i, j)] = (mi, mj, angle)
            bond_dihedral[(j, i)] = (mj, mi, angle)

    # Build parent array: parent[i] = max(neighbor indices < i), or None
    parent = [None] * n_atoms
    for i in range(1, n_atoms):
        neighbors_before = [
            n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() < i
        ]
        if neighbors_before:
            parent[i] = max(neighbors_before)

    child_count = [0] * n_atoms
    first_child_idx = [None] * n_atoms
    # Track torsion used for first child of each atom (for computing offsets)
    first_child_torsion = [None] * n_atoms

    # Place atoms sequentially
    for i in range(n_atoms):
        p = parent[i]
        if p is None:
            continue

        g = parent[p]
        bond_len = _get_bond_length(mol, p, i)
        bond_pi = mol.GetBondBetweenAtoms(p, i)
        in_ring = bond_pi is not None and bond_pi.IsInRing()

        if g is None:
            # Only parent exists — no grandparent
            if child_count[p] == 0:
                coords[i] = coords[p] + np.array([bond_len, 0.0, 0.0])
                first_child_idx[p] = i
            else:
                c1 = first_child_idx[p]
                hyb = mol.GetAtomWithIdx(p).GetHybridization()
                a = _BOND_ANGLES.get(hyb, 109.5)

                # Try AMSR dihedral for the p-c1 bond axis
                placed = False
                if child_count[p] == 1 and (p, c1) in bond_dihedral:
                    mi, mj, amsr_angle = bond_dihedral[(p, c1)]
                    if mi == i:
                        coords[i] = _place_atom(
                            coords[mj],
                            coords[c1],
                            coords[p],
                            bond_len,
                            a,
                            amsr_angle,
                        )
                        placed = True

                if not placed:
                    # Default: use _place_atom with synthetic gg for
                    # consistent reference frame with subtree placement
                    ref_gg = _synthetic_gg(coords, p, c1)
                    if hyb == Chem.HybridizationType.SP2:
                        omega = 180.0
                    else:
                        omega = 120.0 * child_count[p]
                    coords[i] = _place_atom(
                        ref_gg,
                        coords[c1],
                        coords[p],
                        bond_len,
                        a,
                        omega,
                    )
            child_count[p] += 1
            continue

        gg = parent[g]
        angle = _get_bond_angle(mol, g, p, i)
        hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

        # Determine torsion
        torsion = None

        # Check if gg is inside the ring (for default torsion logic)
        if gg is not None:
            bond_gg_g = mol.GetBondBetweenAtoms(gg, g)
            gg_in_ring = bond_gg_g is not None and bond_gg_g.IsInRing()
        else:
            gg_in_ring = True  # synthetic gg: treat as in-ring (works for standalone rings)

        # First child: try AMSR dihedral for the g-p axis
        if child_count[p] == 0:
            if (g, p) in bond_dihedral:
                mi, mj, amsr_angle = bond_dihedral[(g, p)]
                if mj == i and (gg is None or gg == mi):
                    torsion = amsr_angle
            if torsion is None:
                torsion = _default_torsion(mol, p, 0, in_ring, gg_in_ring)
            first_child_torsion[p] = torsion
        else:
            # Subsequent children: offset from first child
            base = first_child_torsion[p]
            if hyb_p == Chem.HybridizationType.SP3:
                torsion = base + 120.0 * child_count[p]
            elif hyb_p == Chem.HybridizationType.SP2:
                torsion = base + 180.0
            else:
                torsion = base + 180.0

        if gg is None:
            ref_gg = _synthetic_gg(coords, g, p)
        else:
            ref_gg = coords[gg]

        coords[i] = _place_atom(ref_gg, coords[g], coords[p], bond_len, angle, torsion)
        child_count[p] += 1

    # Create conformer
    conf = Chem.Conformer(n_atoms)
    conf.Set3D(True)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords[i].tolist())

    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
