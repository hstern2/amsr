"""Z-matrix conformer generation from AMSR dihedrals.

Place every atom one-by-one in DFS order (numerical).  For each atom we need:
  - three reference atoms (A, B, C) already placed
  - bond length C-D (ideal)
  - bond angle B-C-D (ideal)
  - dihedral A-B-C-D (from AMSR, chirality, planarity, or default)
"""

import logging
from typing import Optional

import numpy as np
from rdkit import Chem

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Ideal geometry tables
# ---------------------------------------------------------------------------

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

_HYBRID_ANGLES = {
    Chem.HybridizationType.SP3: 109.5,
    Chem.HybridizationType.SP2: 120.0,
    Chem.HybridizationType.SP: 180.0,
}


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------


def _get_bond_length(mol, i, j):
    bond = mol.GetBondBetweenAtoms(i, j)
    s1, s2 = mol.GetAtomWithIdx(i).GetSymbol(), mol.GetAtomWithIdx(j).GetSymbol()
    if bond.GetIsAromatic():
        bo = 1.5
    elif bond.GetBondType() == Chem.BondType.DOUBLE:
        bo = 2
    elif bond.GetBondType() == Chem.BondType.TRIPLE:
        bo = 3
    else:
        bo = 1
    key = (min(s1, s2), max(s1, s2), bo)
    if key in _BOND_LENGTHS:
        return _BOND_LENGTHS[key]
    return _COVALENT_RADII.get(s1, 1.5) + _COVALENT_RADII.get(s2, 1.5)


def _get_bond_angle(mol, a, b, c):
    """Ideal bond angle a-b-c in degrees."""
    hyb = mol.GetAtomWithIdx(b).GetHybridization()
    if hyb == Chem.HybridizationType.SP2:
        ri = mol.GetRingInfo()
        b_ab = mol.GetBondBetweenAtoms(a, b)
        b_bc = mol.GetBondBetweenAtoms(b, c)
        if b_ab is not None and b_bc is not None:
            common = set(ri.BondRingSizes(b_ab.GetIdx())) & set(ri.BondRingSizes(b_bc.GetIdx()))
            if common:
                n = min(common)
                return (n - 2) * 180.0 / n
    return _HYBRID_ANGLES.get(hyb, 109.5)


def _place_atom(A, B, C, d, theta_deg, omega_deg):
    """Place atom D given three reference points A, B, C.

    d:         bond length C-D
    theta_deg: bond angle  B-C-D
    omega_deg: dihedral    A-B-C-D
    """
    theta = np.radians(theta_deg)
    omega = np.radians(omega_deg)
    BC = C - B
    bc = np.linalg.norm(BC)
    BC = BC / bc if bc > 1e-10 else np.array([1.0, 0.0, 0.0])
    AB = B - A
    n = np.cross(AB, BC)
    nn = np.linalg.norm(n)
    if nn < 1e-10:
        perp = np.array([1.0, 0.0, 0.0]) if abs(BC[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        n = np.cross(BC, perp)
        n /= np.linalg.norm(n)
    else:
        n /= nn
    m = np.cross(n, BC)
    return C + d * (
        -np.cos(theta) * BC + np.sin(theta) * np.cos(omega) * m + np.sin(theta) * np.sin(omega) * n
    )


def _measure_torsion(p0, p1, p2, p3):
    """Torsion angle (degrees) from four 3-D points."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    n1n, n2n = np.linalg.norm(n1), np.linalg.norm(n2)
    if n1n < 1e-10 or n2n < 1e-10:
        return 0.0
    n1, n2 = n1 / n1n, n2 / n2n
    b2u = b2 / np.linalg.norm(b2)
    return np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2u), np.dot(n1, n2)))


def _synthetic_ref(coords, b, c):
    """Synthetic reference point A when no great-grandparent exists."""
    bc = coords[c] - coords[b]
    bc_n = np.linalg.norm(bc)
    bc = bc / bc_n if bc_n > 1e-10 else np.array([1.0, 0.0, 0.0])
    perp = np.array([1.0, 0.0, 0.0]) if abs(bc[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp = np.cross(bc, perp)
    perp /= np.linalg.norm(perp)
    return coords[b] - perp


# ---------------------------------------------------------------------------
# Dihedral selection
# ---------------------------------------------------------------------------


def _choose_dihedral(
    mol, i, p, g, gg, nth_child, first_child_torsion, in_ring, coords, bond_dihedral
):
    """Choose the dihedral angle for placing atom i from parent p.

    Reference frame: gg-g-p-i (or ref_override-g-p-i).
    Returns (dihedral_angle, ref_override_or_None).
    """
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    # --- first child: use AMSR dihedral if available -------------------------
    if nth_child == 0:
        if (g, p) in bond_dihedral:
            mi, mj, angle = bond_dihedral[(g, p)]
            if mj == i:
                return angle, mi if mi != gg else None
            # AMSR dihedral exists but for a different neighbor mj.
            # Offset from mj's torsion based on hybridization and chirality.
            if np.any(coords[mj]) or mj == 0:
                actual_mj = _measure_torsion(
                    coords[gg] if gg is not None else _synthetic_ref(coords, g, p),
                    coords[g],
                    coords[p],
                    coords[mj],
                )
                if hyb_p == Chem.HybridizationType.SP2:
                    return actual_mj + 180.0, None
                chiral = mol.GetAtomWithIdx(p).GetChiralTag()
                if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    return actual_mj - 120.0, None
                return actual_mj + 120.0, None

        # Default torsion: check if gg-g-p-i are all in the same ring
        same_ring = False
        if gg is not None and in_ring:
            for ring in mol.GetRingInfo().AtomRings():
                if gg in ring and g in ring and p in ring and i in ring:
                    same_ring = True
                    break
        return _default_torsion(hyb_p, in_ring, same_ring, 0), None

    # --- subsequent children: offset from first child -----------------------
    base = first_child_torsion[p] if first_child_torsion[p] is not None else 0.0
    if hyb_p == Chem.HybridizationType.SP2:
        return base + 180.0, None
    if hyb_p == Chem.HybridizationType.SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return base + 120.0 * nth_child, None
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return base - 120.0 * nth_child, None
        # Non-chiral: heuristic
        if base > 0:
            return base - 120.0 * nth_child, None
        return base + 120.0 * nth_child, None
    return base + 180.0, None


def _default_torsion(hyb, in_ring, same_ring, nth_child):
    """Default torsion when no AMSR dihedral is available."""
    if in_ring and same_ring:
        return 0.0
    if hyb == Chem.HybridizationType.SP2:
        if in_ring:
            return 180.0  # different ring at junction
        return 180.0
    if hyb == Chem.HybridizationType.SP3:
        return 180.0
    return 180.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
) -> Chem.Mol:
    """Generate 3D conformer by z-matrix atom-by-atom placement.

    Atoms placed in DFS order (numerical).  parent[i] = max(neighbor < i).
    """
    n = mol.GetNumAtoms()
    if n == 0:
        return mol

    coords = np.zeros((n, 3))

    bond_dihedral = {}
    if dihedral:
        for (mi, i, j, mj), angle in dihedral.items():
            bond_dihedral[(i, j)] = (mi, mj, angle)
            bond_dihedral[(j, i)] = (mj, mi, angle)

    # parent[i] = max(neighbor < i), matching AMSR DFS order
    parent = [None] * n
    for i in range(1, n):
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() < i]
        if nbrs:
            parent[i] = max(nbrs)

    placed = set()
    child_count = [0] * n
    first_child_torsion = [None] * n
    first_child_idx = [None] * n

    for i in range(n):
        p = parent[i]
        if p is None:
            placed.add(i)
            continue

        g = parent[p]
        bond_len = _get_bond_length(mol, p, i)
        bp = mol.GetBondBetweenAtoms(p, i)
        in_ring = bp is not None and bp.IsInRing()

        if g is None:
            if child_count[p] == 0:
                coords[i] = coords[p] + np.array([bond_len, 0.0, 0.0])
            else:
                c1 = first_child_idx[p]
                hyb = mol.GetAtomWithIdx(p).GetHybridization()
                ang = _HYBRID_ANGLES.get(hyb, 109.5)
                done = False
                if child_count[p] == 1 and (p, c1) in bond_dihedral:
                    mi, mj, amsr_a = bond_dihedral[(p, c1)]
                    if mi == i:
                        coords[i] = _place_atom(
                            coords[mj], coords[c1], coords[p], bond_len, ang, amsr_a
                        )
                        done = True
                if not done:
                    ref = _synthetic_ref(coords, p, c1)
                    omega = 180.0 if hyb == Chem.HybridizationType.SP2 else 120.0 * child_count[p]
                    coords[i] = _place_atom(ref, coords[c1], coords[p], bond_len, ang, omega)
        else:
            gg = parent[g]
            bond_angle = _get_bond_angle(mol, g, p, i)

            torsion, ref_override = _choose_dihedral(
                mol,
                i,
                p,
                g,
                gg,
                child_count[p],
                first_child_torsion,
                in_ring,
                coords,
                bond_dihedral,
            )
            if ref_override is not None:
                ref = coords[ref_override]
            elif gg is not None:
                ref = coords[gg]
            elif child_count[p] == 0 and (g, p) in bond_dihedral:
                mi, mj, _ = bond_dihedral[(g, p)]
                ref = coords[mi] if mj == i else _synthetic_ref(coords, g, p)
            else:
                ref = _synthetic_ref(coords, g, p)
            coords[i] = _place_atom(ref, coords[g], coords[p], bond_len, bond_angle, torsion)

            if child_count[p] == 0:
                std_ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                first_child_torsion[p] = _measure_torsion(std_ref, coords[g], coords[p], coords[i])

        if first_child_idx[p] is None:
            first_child_idx[p] = i
        child_count[p] += 1
        placed.add(i)

    # Build RDKit conformer
    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
