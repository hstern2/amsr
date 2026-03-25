"""Conformer generation from AMSR dihedrals.

Ring atoms are placed by optimizing a Cartesian-space cost function that
enforces ideal bond lengths, bond angles, planarity at SP2 centers,
chirality at SP3 centers, and AMSR dihedral restraints.

Non-ring atoms are placed sequentially via z-matrix from their parent chain.

Geometry primitives and cost-function components use only numpy arrays
(no RDKit), making them suitable for reimplementation in C/C++.
"""

import math
from typing import Optional

import numpy as np
from rdkit import Chem

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
    ("N", "S", 1.5): 1.63,
    ("N", "S", 2): 1.54,
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

SP2 = Chem.HybridizationType.SP2
SP3 = Chem.HybridizationType.SP3
CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW


def _is_placed(coords, k):
    """True if atom k has been placed (has nonzero coords, or is atom 0)."""
    return k == 0 or np.any(coords[k])


def _sp3_offset(mol, p, g, gg, i, k, coords):
    """Compute torsion and alternatives for atom i offset from placed neighbor k.

    Uses graph-order chirality convention at SP3 parent p.
    Returns (torsion, [alternatives]).
    """
    ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
    actual_k = measure_torsion(ref, coords[g], coords[p], coords[k])
    children = [nb.GetIdx() for nb in mol.GetAtomWithIdx(p).GetNeighbors() if nb.GetIdx() != g]
    n_children = len(children)
    steps = (children.index(i) - children.index(k)) % n_children
    sign = -1.0 if mol.GetAtomWithIdx(p).GetChiralTag() == CCW else 1.0
    torsion = actual_k + sign * 120.0 * steps
    # Alternatives: all other 120° positions around k.
    alts = [actual_k + sign * 120.0 * s for s in range(1, n_children) if s != steps]
    # Also include the opposite-chirality positions.
    alts.extend(actual_k - sign * 120.0 * s for s in range(1, n_children))
    return torsion, alts


# ---------------------------------------------------------------------------
# Geometry primitives (pure numpy — C-portable)
# ---------------------------------------------------------------------------


def _cross3(a, b):
    """Cross product for 3-element arrays (avoids numpy.cross overhead)."""
    return np.array(
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    )


def _norm3(v):
    """Euclidean norm for 3-element array (avoids numpy.linalg.norm overhead)."""
    return np.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def place_atom(A, B, C, d, theta_deg, omega_deg):
    """Place atom D given refs A, B, C, bond length d, angle B-C-D, torsion A-B-C-D."""
    theta = math.radians(theta_deg)
    omega = math.radians(omega_deg)
    BC = C - B
    bc = _norm3(BC)
    if bc > 1e-10:
        BC = BC / bc
    else:
        BC = np.array([1.0, 0.0, 0.0])
    n = _cross3(B - A, BC)
    nn = _norm3(n)
    if nn < 1e-10:
        perp = np.array([1.0, 0.0, 0.0]) if abs(BC[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        n = _cross3(BC, perp)
        n /= _norm3(n)
    else:
        n /= nn
    m = _cross3(n, BC)
    st = math.sin(theta)
    return C + d * (-math.cos(theta) * BC + st * math.cos(omega) * m + st * math.sin(omega) * n)


def measure_torsion(p0, p1, p2, p3):
    """Torsion angle (degrees) for four 3-D points."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = _cross3(b1, b2), _cross3(b2, b3)
    n1n, n2n = _norm3(n1), _norm3(n2)
    if n1n < 1e-10 or n2n < 1e-10:
        return 0.0
    n1, n2 = n1 / n1n, n2 / n2n
    return np.degrees(np.arctan2(np.dot(_cross3(n1, n2), b2 / _norm3(b2)), np.dot(n1, n2)))


def measure_angle(coords, a, b, c):
    """Angle a-b-c (degrees) from coordinates."""
    v1, v2 = coords[a] - coords[b], coords[c] - coords[b]
    cos_a = np.dot(v1, v2) / (_norm3(v1) * _norm3(v2) + 1e-10)
    return math.degrees(math.acos(max(-1.0, min(1.0, float(cos_a)))))


# ---------------------------------------------------------------------------
# Batch geometry primitives (vectorized numpy)
# ---------------------------------------------------------------------------


def _batch_norm3(v):
    """Euclidean norm for Nx3 array, returns shape (N,)."""
    return np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)


def _batch_cross3(a, b):
    """Cross product for Nx3 arrays, returns Nx3."""
    return np.column_stack(
        [
            a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
            a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
            a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
        ]
    )


def _batch_measure_torsion(p0, p1, p2, p3):
    """Torsion angles (degrees) for N sets of four 3-D points (each Nx3)."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1 = _batch_cross3(b1, b2)
    n2 = _batch_cross3(b2, b3)
    n1n = _batch_norm3(n1)
    n2n = _batch_norm3(n2)
    # Avoid division by zero
    safe = (n1n > 1e-10) & (n2n > 1e-10)
    result = np.zeros(len(p0))
    if not np.any(safe):
        return result
    n1s = n1[safe] / n1n[safe, None]
    n2s = n2[safe] / n2n[safe, None]
    b2s = b2[safe]
    b2n = _batch_norm3(b2s)
    b2s = b2s / b2n[:, None]
    cross_n = _batch_cross3(n1s, n2s)
    sin_val = np.sum(cross_n * b2s, axis=1)
    cos_val = np.sum(n1s * n2s, axis=1)
    result[safe] = np.degrees(np.arctan2(sin_val, cos_val))
    return result


def _batch_measure_angle(coords, triples):
    """Angles a-b-c (degrees) from Nx3 index array. Returns shape (N,)."""
    v1 = coords[triples[:, 0]] - coords[triples[:, 1]]
    v2 = coords[triples[:, 2]] - coords[triples[:, 1]]
    cos_a = np.sum(v1 * v2, axis=1) / (_batch_norm3(v1) * _batch_norm3(v2) + 1e-10)
    return np.degrees(np.arccos(np.clip(cos_a, -1, 1)))


# ---------------------------------------------------------------------------
# RDKit helpers (not C-portable)
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


_ELEMENT_ANGLES = {
    # Heteroatoms with bond angles smaller than the regular polygon formula
    "S": {SP2: 102.0, SP3: 96.0},
    "Se": {SP2: 100.0, SP3: 95.0},
}


def _get_bond_angle(mol, a, b, c):
    """Ideal bond angle a-b-c in degrees."""
    atom_b = mol.GetAtomWithIdx(b)
    hyb = atom_b.GetHybridization()
    sym = atom_b.GetSymbol()
    if sym in _ELEMENT_ANGLES and hyb in _ELEMENT_ANGLES[sym]:
        return _ELEMENT_ANGLES[sym][hyb]
    if hyb == SP2:
        ri = mol.GetRingInfo()
        b_ab = mol.GetBondBetweenAtoms(a, b)
        b_bc = mol.GetBondBetweenAtoms(b, c)
        if b_ab is not None and b_bc is not None:
            common = set(ri.BondRingSizes(b_ab.GetIdx())) & set(ri.BondRingSizes(b_bc.GetIdx()))
            if common:
                n = min(common)
                return (n - 2) * 180.0 / n
    return _HYBRID_ANGLES.get(hyb, 109.5)


def _synthetic_ref(coords, b, c):
    """Synthetic reference point when no great-grandparent exists."""
    bc = coords[c] - coords[b]
    bc_n = _norm3(bc)
    if bc_n > 1e-10:
        bc = bc / bc_n
    else:
        bc = np.array([1.0, 0.0, 0.0])
    perp = np.array([1.0, 0.0, 0.0]) if abs(bc[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp = _cross3(bc, perp)
    perp /= _norm3(perp)
    return coords[b] - perp


def _ref_point(coords, parent, g, p, i, bond_dihedral):
    """Compute the reference point (A in A-B-C-D) for placing atom i."""
    gg = parent[g]
    if gg is not None:
        return coords[gg], gg
    if (g, p) in bond_dihedral:
        mi, mj, _ = bond_dihedral[(g, p)]
        if mj == i:
            return coords[mi], mi
    return _synthetic_ref(coords, g, p), None


def _find_ring_systems(mol):
    """Return list of sets of atom indices forming connected ring systems."""
    systems = []
    for ring in mol.GetRingInfo().AtomRings():
        new = set(ring)
        merged = []
        for s in systems:
            if new & s:
                new |= s
            else:
                merged.append(s)
        merged.append(new)
        systems = merged
    return systems


# ============================================================
# Cartesian ring geometry: data collection (RDKit-dependent)
# ============================================================


def _collect_ring_bonds(mol, sys_set, fixed):
    """Collect all bonds within the ring system and to fixed neighbors.

    Returns (pairs, ideal_lengths) — both numpy arrays.
    """
    pairs, ideals = [], []
    for a in sorted(sys_set):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b in sys_set and b > a:
                pairs.append((a, b))
                ideals.append(_get_bond_length(mol, a, b))
            elif b in fixed:
                pairs.append((a, b))
                ideals.append(_get_bond_length(mol, a, b))
    return np.array(pairs, dtype=int), np.array(ideals)


def _collect_ring_angles(mol, sys_set, fixed):
    """Collect all bond angle triples involving ring system atoms.

    Includes angles at ring atoms AND angles at fixed atoms that have
    two ring-system neighbors.
    Returns (triples, ideal_angles) — numpy arrays.
    """
    available = sys_set | set(fixed)
    triples, ideals = [], []
    # Angles at ring atoms
    for b in sorted(sys_set):
        nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(b).GetNeighbors() if nb.GetIdx() in available
        ]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                triples.append((a, b, c))
                ideals.append(_get_bond_angle(mol, a, b, c))
    # Angles at fixed atoms with >=2 ring neighbors
    for fb in sorted(fixed):
        ring_nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(fb).GetNeighbors() if nb.GetIdx() in sys_set
        ]
        if len(ring_nbrs) >= 2:
            for ia in range(len(ring_nbrs)):
                for ic in range(ia + 1, len(ring_nbrs)):
                    a, c = ring_nbrs[ia], ring_nbrs[ic]
                    triples.append((a, fb, c))
                    ideals.append(_get_bond_angle(mol, a, fb, c))
    return np.array(triples, dtype=int) if triples else np.empty((0, 3), dtype=int), np.array(
        ideals
    )


def _collect_planar_atoms(mol, sys_set, fixed):
    """Collect planarity constraints for SP2 atoms in the ring system.

    For each SP2 atom with 3+ neighbors whose coords are available,
    returns (center, a, b, c) tuples — all four should be coplanar.
    """
    available = sys_set | set(fixed)
    groups = []
    for j in sorted(sys_set):
        if mol.GetAtomWithIdx(j).GetHybridization() != SP2:
            continue
        nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() in available
        ]
        if len(nbrs) >= 3:
            groups.append((j, nbrs[0], nbrs[1], nbrs[2]))
    return np.array(groups, dtype=int) if groups else np.empty((0, 4), dtype=int)


def _collect_chiral_atoms(mol, sys_set, fixed):
    """Collect chirality constraints for SP3 chiral atoms in the ring system.

    Returns Nx5 array: (center, a, b, c, sign) where sign is +1 (CW) or -1 (CCW).
    """
    available = sys_set | set(fixed)
    result = []
    for j in sorted(sys_set):
        atom = mol.GetAtomWithIdx(j)
        chiral = atom.GetChiralTag()
        if chiral not in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in available]
        if len(nbrs) < 3:
            continue
        sign = 1 if chiral == CW else -1
        result.append((j, nbrs[0], nbrs[1], nbrs[2], sign))
    return np.array(result, dtype=int) if result else np.empty((0, 5), dtype=int)


def _collect_ring_dihedrals(mol, sys_set, bond_dihedral, fixed):
    """Collect AMSR dihedral restraints for ring bonds and boundary bonds.

    Includes dihedrals for bonds within the ring system and bonds
    connecting ring atoms to fixed (placed) atoms.
    Returns (quads, targets) where quads is Nx4 int array and targets is N float.
    """
    available = sys_set | set(fixed)
    quads, targets = [], []
    seen = set()
    for a in sorted(sys_set):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b not in sys_set and b not in fixed:
                continue
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            seen.add(key)
            bond_key = (a, b) if (a, b) in bond_dihedral else (b, a)
            if bond_key in bond_dihedral:
                mi, mj, angle = bond_dihedral[bond_key]
                if mi in available and mj in available:
                    quads.append((mi, bond_key[0], bond_key[1], mj))
                    targets.append(float(angle))
    return (
        np.array(quads, dtype=int) if quads else np.empty((0, 4), dtype=int),
        np.array(targets) if targets else np.empty(0),
    )


def _collect_ez_constraints(mol, sys_set, fixed):
    """Collect cis/trans dihedral constraints for E/Z double bonds.

    For each double bond with E/Z stereo where at least one end is in
    the ring system, add a dihedral constraint (0° for Z, 180° for E).
    Returns (quads, targets) — Nx4 int, N float.
    """
    available = sys_set | set(fixed)
    quads, targets = [], []
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i not in sys_set and j not in sys_set:
            continue
        # Get stereo atoms (the reference atoms for E/Z)
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) < 2:
            continue
        si, sj = stereo_atoms[0], stereo_atoms[1]
        if si not in available or sj not in available:
            continue
        target = 0.0 if stereo == Chem.BondStereo.STEREOZ else 180.0
        quads.append((si, i, j, sj))
        targets.append(target)
    return (
        np.array(quads, dtype=int) if quads else np.empty((0, 4), dtype=int),
        np.array(targets) if targets else np.empty(0),
    )


# ============================================================
# Cartesian ring geometry: residual functions (vectorized numpy, C-portable)
# ============================================================


def _resolve_coords(indices, all_coords, idx_map, fixed):
    """Look up Nx3 coordinates for an array of atom indices.

    Atoms in idx_map use all_coords (the optimization variable reshaped);
    atoms in fixed use their stored coordinates.
    """
    n = len(indices)
    out = np.empty((n, 3))
    for k in range(n):
        a = int(indices[k])
        if a in idx_map:
            out[k] = all_coords[idx_map[a]]
        else:
            out[k] = fixed[a]
    return out


def _bond_residuals(x_3d, idx_map, fixed, pairs, ideal_lengths, w):
    """Residuals: w * (|r_i - r_j| - d_ij) for each bond."""
    ri = _resolve_coords(pairs[:, 0], x_3d, idx_map, fixed)
    rj = _resolve_coords(pairs[:, 1], x_3d, idx_map, fixed)
    dists = _batch_norm3(ri - rj)
    return w * (dists - ideal_lengths)


def _angle_residuals(x_3d, idx_map, fixed, triples, ideal_angles, w):
    """Residuals in Angstroms: w * 2 * L * sin(delta_angle / 2).

    Converts angular error to Cartesian displacement, making it
    commensurable with bond-length residuals.
    """
    ra = _resolve_coords(triples[:, 0], x_3d, idx_map, fixed)
    rb = _resolve_coords(triples[:, 1], x_3d, idx_map, fixed)
    rc = _resolve_coords(triples[:, 2], x_3d, idx_map, fixed)
    v1, v2 = ra - rb, rc - rb
    n1, n2 = _batch_norm3(v1), _batch_norm3(v2)
    L = 0.5 * (n1 + n2)
    cos_a = np.sum(v1 * v2, axis=1) / (n1 * n2 + 1e-10)
    actual_rad = np.arccos(np.clip(cos_a, -1, 1))
    ideal_rad = np.radians(ideal_angles)
    delta_half = (actual_rad - ideal_rad) / 2.0
    return w * 2.0 * L * np.sin(delta_half)


def _planarity_residuals(x_3d, idx_map, fixed, groups, w):
    """Residuals: w * normalized_volume for each SP2 center.

    Volume = (a-j) . ((b-j) x (c-j)), normalized by product of bond lengths.
    Zero when all four atoms are coplanar.
    """
    rj = _resolve_coords(groups[:, 0], x_3d, idx_map, fixed)
    ra = _resolve_coords(groups[:, 1], x_3d, idx_map, fixed)
    rb = _resolve_coords(groups[:, 2], x_3d, idx_map, fixed)
    rc = _resolve_coords(groups[:, 3], x_3d, idx_map, fixed)
    v1, v2, v3 = ra - rj, rb - rj, rc - rj
    cross = _batch_cross3(v2, v3)
    vol = np.sum(v1 * cross, axis=1)
    norm = _batch_norm3(v1) * _batch_norm3(v2) * _batch_norm3(v3) + 1e-10
    return w * vol / norm


def _chirality_residuals(x_3d, idx_map, fixed, chiral_info, w):
    """Residuals penalizing wrong-sign volume at chiral centers."""
    rj = _resolve_coords(chiral_info[:, 0], x_3d, idx_map, fixed)
    ra = _resolve_coords(chiral_info[:, 1], x_3d, idx_map, fixed)
    rb = _resolve_coords(chiral_info[:, 2], x_3d, idx_map, fixed)
    rc = _resolve_coords(chiral_info[:, 3], x_3d, idx_map, fixed)
    sign = chiral_info[:, 4].astype(float)
    v1, v2, v3 = ra - rj, rb - rj, rc - rj
    vol = np.sum(v1 * _batch_cross3(v2, v3), axis=1)
    return w * np.maximum(0.0, -sign * vol)


def _dihedral_residuals(x_3d, idx_map, fixed, quads, targets, w):
    """Residuals: w * angular_diff for each dihedral restraint."""
    p0 = _resolve_coords(quads[:, 0], x_3d, idx_map, fixed)
    p1 = _resolve_coords(quads[:, 1], x_3d, idx_map, fixed)
    p2 = _resolve_coords(quads[:, 2], x_3d, idx_map, fixed)
    p3 = _resolve_coords(quads[:, 3], x_3d, idx_map, fixed)
    actual = _batch_measure_torsion(p0, p1, p2, p3)
    diff = (actual - targets + 180.0) % 360.0 - 180.0
    return w * diff


# ============================================================
# Cartesian ring geometry: RDKit embedding for initialization
# ============================================================


def _rdkit_embed(mol, n_confs=1, seed=42):
    """Embed molecule with RDKit distance geometry (adds/removes Hs internally).

    Returns list of Nx3 coordinate arrays (one per conformer), or empty list
    on failure.  Heavy-atom indices match the input mol.
    """
    from rdkit.Chem import AllChem

    mol_h = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(
        mol_h, numConfs=n_confs, randomSeed=seed, enforceChirality=True
    )
    results = []
    for cid in cids:
        conf = mol_h.GetConformer(cid)
        c = np.zeros((mol.GetNumAtoms(), 3))
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            c[i] = [pos.x, pos.y, pos.z]
        results.append(c)
    return results


# ============================================================
# Cartesian ring geometry: optimizer
# ============================================================

# Weights for residual terms
_W_BOND = 5.0
_W_ANGLE = 2.0
_W_PLANAR = 3.0
_W_CHIRAL = 2.0
_W_DIHEDRAL = 0.05


_W_EZ = 0.3


def _ring_system_residuals(
    x,
    idx_map,
    fixed,
    bonds,
    ideal_lengths,
    angle_triples,
    ideal_angles,
    planar_groups,
    chiral_info,
    dih_quads,
    dih_targets,
    ez_quads,
    ez_targets,
):
    """Combined residual vector for ring system Cartesian optimization."""
    x_3d = x.reshape(-1, 3)
    parts = []
    if len(bonds):
        parts.append(_bond_residuals(x_3d, idx_map, fixed, bonds, ideal_lengths, _W_BOND))
    if len(angle_triples):
        parts.append(_angle_residuals(x_3d, idx_map, fixed, angle_triples, ideal_angles, _W_ANGLE))
    if len(planar_groups):
        parts.append(_planarity_residuals(x_3d, idx_map, fixed, planar_groups, _W_PLANAR))
    if len(chiral_info):
        parts.append(_chirality_residuals(x_3d, idx_map, fixed, chiral_info, _W_CHIRAL))
    if len(dih_quads):
        parts.append(_dihedral_residuals(x_3d, idx_map, fixed, dih_quads, dih_targets, _W_DIHEDRAL))
    if len(ez_quads):
        parts.append(_dihedral_residuals(x_3d, idx_map, fixed, ez_quads, ez_targets, _W_EZ))
    if not parts:
        return np.array([0.0])
    return np.concatenate(parts)


def _optimize_ring_system(mol, system_atoms, all_rings, bond_dihedral, coords, placed, parent):
    """Place and optimize ring system atoms in Cartesian space.

    Initializes ring atoms as regular polygons, then minimizes a cost
    function enforcing ideal bonds, angles, planarity, chirality, and
    AMSR dihedral restraints.
    """
    from scipy.optimize import least_squares

    sys_set = set(system_atoms)
    sys_list = sorted(system_atoms)
    n_sys = len(sys_list)
    idx_map = {atom: i for i, atom in enumerate(sys_list)}

    # Fixed coordinates: placed atoms adjacent to the ring system
    fixed: dict[int, np.ndarray] = {}
    for a in sys_set:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b not in sys_set and b in placed:
                fixed[b] = coords[b].copy()

    # Collect constraint data
    bonds, ideal_lengths = _collect_ring_bonds(mol, sys_set, fixed)
    angle_triples, ideal_angles = _collect_ring_angles(mol, sys_set, fixed)
    planar_groups = _collect_planar_atoms(mol, sys_set, fixed)
    chiral_info = _collect_chiral_atoms(mol, sys_set, fixed)
    dih_quads, dih_targets = _collect_ring_dihedrals(mol, sys_set, bond_dihedral, fixed)
    ez_quads, ez_targets = _collect_ez_constraints(mol, sys_set, fixed)

    def residual_fn(x):
        return _ring_system_residuals(
            x,
            idx_map,
            fixed,
            bonds,
            ideal_lengths,
            angle_triples,
            ideal_angles,
            planar_groups,
            chiral_info,
            dih_quads,
            dih_targets,
            ez_quads,
            ez_targets,
        )

    # Initial coordinates come from the embedding already stored in coords.
    x0 = np.zeros(3 * n_sys)
    for a in sys_list:
        k = idx_map[a]
        x0[3 * k : 3 * k + 3] = coords[a]

    result = least_squares(residual_fn, x0, method="trf", ftol=1e-10, xtol=1e-10, gtol=1e-10)
    best_cost = result.cost
    best_result = result

    # If cost is still high, try more RDKit embeddings
    if best_cost > 0.5 and n_sys > 0:
        extra = _rdkit_embed(mol, n_confs=4, seed=123)
        for embed_coords in extra:
            x0 = np.zeros(3 * n_sys)
            for a in sys_list:
                x0[3 * idx_map[a] : 3 * idx_map[a] + 3] = embed_coords[a]
            r = least_squares(residual_fn, x0, method="trf", ftol=1e-10, xtol=1e-10, gtol=1e-10)
            if r.cost < best_cost:
                best_cost = r.cost
                best_result = r

    # Copy optimized coordinates back
    x_opt = best_result.x.reshape(-1, 3)
    for a in sys_list:
        coords[a] = x_opt[idx_map[a]]


# ============================================================
# Chain-atom dihedral selection (AMSR dihedrals + chirality)
# ============================================================


def _choose_chain_dihedral(
    mol, i, p, g, gg, nth_child, first_child_torsion, first_child_idx, coords, bond_dihedral
):
    """Choose torsion for a non-ring atom.

    Returns (torsion, ref_override, alternatives).
    """
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    # AMSR dihedral on backward bond (g, p)
    if (g, p) in bond_dihedral:
        mi, mj, angle = bond_dihedral[(g, p)]
        if mj == i:
            return angle, mi if mi != gg else None, []
        if nth_child == 0 and _is_placed(coords, mj):
            if hyb_p == SP2:
                ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                actual_mj = measure_torsion(ref, coords[g], coords[p], coords[mj])
                return actual_mj + 180.0, None, []
            torsion, alts = _sp3_offset(mol, p, g, gg, i, mj, coords)
            return torsion, None, alts

    # First child: SP3 offset from placed neighbor
    if nth_child == 0:
        if hyb_p == SP3:
            for nb in mol.GetAtomWithIdx(p).GetNeighbors():
                k = nb.GetIdx()
                if k != g and k != i and _is_placed(coords, k):
                    torsion, alts = _sp3_offset(mol, p, g, gg, i, k, coords)
                    return torsion, None, alts
        return 180.0, None, [0.0]

    # Subsequent children: offset from first child
    base = first_child_torsion[p] if first_child_torsion[p] is not None else 0.0
    if hyb_p == SP2:
        return base + 180.0, None, []
    if hyb_p == SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral in (CW, CCW):
            fc = first_child_idx[p]
            children = [
                nb.GetIdx() for nb in mol.GetAtomWithIdx(p).GetNeighbors() if nb.GetIdx() != g
            ]
            swapped = (
                len(children) >= 2
                and fc is not None
                and fc in children
                and i in children
                and children.index(fc) > children.index(i)
            )
            if chiral == CW:
                sign = -1 if swapped else 1
            else:
                sign = 1 if swapped else -1
            return base + sign * 120.0 * nth_child, None, [base - sign * 120.0 * nth_child]
        if g is not None and (g, p) in bond_dihedral:
            mi_bwd, mj_bwd, angle_bwd = bond_dihedral[(g, p)]
            if mj_bwd != i and not _is_placed(coords, mj_bwd) and mi_bwd == gg:
                for s in (1, -1):
                    candidate = base + s * 120.0 * nth_child
                    diff = abs((candidate - angle_bwd + 180.0) % 360.0 - 180.0)
                    if diff > 30.0:
                        return candidate, None, [base - s * 120.0 * nth_child]
        sign = -1 if base > 0 else 1
        return base + sign * 120.0 * nth_child, None, [base - sign * 120.0 * nth_child]
    return base + 180.0, None, []


# ============================================================
# Chain-atom placement (z-matrix)
# ============================================================


def _has_collision(mol, j, coords, placed, threshold=0.5):
    """Check if atom j collides with any placed non-bonded atom."""
    for other in placed:
        if other == j or mol.GetBondBetweenAtoms(j, other) is not None:
            continue
        if _norm3(coords[j] - coords[other]) < threshold:
            return True
    return False


def _place_chain_atom(
    mol,
    i,
    coords,
    parent,
    child_count,
    first_child_torsion,
    first_child_idx,
    bond_dihedral,
    torsion_override=None,
):
    """Place a non-ring atom via z-matrix. Returns alternative torsions."""
    p = parent[i]
    if p is None:
        return []
    g = parent[p]
    bond_len = _get_bond_length(mol, p, i)
    alternatives: list[float] = []

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
                    coords[i] = place_atom(coords[mj], coords[c1], coords[p], bond_len, ang, amsr_a)
                    done = True
            if not done:
                # Prefer placed neighbors of p over synthetic reference —
                # essential when p is in a ring (synthetic ref is perpendicular
                # to the ring plane, giving degenerate chain atom placement).
                other_placed = [
                    nb.GetIdx()
                    for nb in mol.GetAtomWithIdx(p).GetNeighbors()
                    if nb.GetIdx() != c1 and nb.GetIdx() != i and _is_placed(coords, nb.GetIdx())
                ]
                if other_placed:
                    k = other_placed[0]
                    if hyb == SP2:
                        omega = 180.0
                    else:
                        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
                        sign = -1.0 if chiral == CCW else 1.0
                        children_of_p = [
                            nb.GetIdx()
                            for nb in mol.GetAtomWithIdx(p).GetNeighbors()
                            if nb.GetIdx() != c1
                        ]
                        if i in children_of_p and k in children_of_p:
                            steps = (children_of_p.index(i) - children_of_p.index(k)) % len(
                                children_of_p
                            )
                        else:
                            steps = child_count[p]
                        omega = sign * 120.0 * steps
                    coords[i] = place_atom(coords[k], coords[c1], coords[p], bond_len, ang, omega)
                else:
                    ref = _synthetic_ref(coords, p, c1)
                    omega = 180.0 if hyb == SP2 else 120.0 * child_count[p]
                    coords[i] = place_atom(ref, coords[c1], coords[p], bond_len, ang, omega)
    else:
        gg = parent[g]
        bond_angle = _get_bond_angle(mol, g, p, i)
        if torsion_override is not None:
            torsion = torsion_override
            ref_override = None
        else:
            torsion, ref_override, alternatives = _choose_chain_dihedral(
                mol,
                i,
                p,
                g,
                gg,
                child_count[p],
                first_child_torsion,
                first_child_idx,
                coords,
                bond_dihedral,
            )
        if ref_override is not None and _is_placed(coords, ref_override):
            ref = coords[ref_override]
        else:
            ref, _ = _ref_point(coords, parent, g, p, i, bond_dihedral)
        coords[i] = place_atom(ref, coords[g], coords[p], bond_len, bond_angle, torsion)

        if child_count[p] == 0:
            std_ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
            first_child_torsion[p] = measure_torsion(std_ref, coords[g], coords[p], coords[i])

    if first_child_idx[p] is None:
        first_child_idx[p] = i
    child_count[p] += 1
    return alternatives


# ============================================================
# Public API
# ============================================================


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
) -> Chem.Mol:
    """Generate 3D conformer.

    1. Cartesian-optimize the ring core (rings + bridging chains) from
       an RDKit distance-geometry embedding.
    2. Z-matrix place branch atoms off the core using AMSR dihedrals.
    """
    n = mol.GetNumAtoms()
    if n == 0:
        return mol

    coords = np.zeros((n, 3))

    bond_dihedral: dict[tuple[int, int], tuple[int, int, int]] = {}
    if dihedral:
        for (mi, i, j, mj), angle in dihedral.items():
            bond_dihedral[(i, j)] = (mi, mj, angle)
            bond_dihedral[(j, i)] = (mj, mi, angle)

    parent: list[Optional[int]] = [None] * n
    for i in range(1, n):
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() < i]
        if nbrs:
            parent[i] = max(nbrs)

    child_count = [0] * n
    first_child_torsion: list[Optional[float]] = [None] * n
    first_child_idx: list[Optional[int]] = [None] * n

    ring_systems = _find_ring_systems(mol)
    ring_atoms: set[int] = set()
    for sys in ring_systems:
        ring_atoms.update(sys)

    all_rings = [tuple(r) for r in mol.GetRingInfo().AtomRings()]

    # --- Phase 1: Cartesian optimize the ring core ---
    # The "core" is the minimal connected subgraph spanning all ring atoms:
    # ring atoms + chain atoms bridging between ring systems, but NOT
    # terminal branches.  Identified by iteratively pruning non-ring leaves.
    core_atoms: set[int] = set()
    if ring_atoms:
        # Core = ring atoms + bridging chain atoms between ring systems.
        # Identified by pruning non-ring leaves iteratively.
        core_atoms = set(range(n))
        changed = True
        while changed:
            changed = False
            for a in list(core_atoms):
                if a in ring_atoms:
                    continue
                nbrs_in_core = sum(
                    1 for nb in mol.GetAtomWithIdx(a).GetNeighbors() if nb.GetIdx() in core_atoms
                )
                if nbrs_in_core <= 1:
                    core_atoms.discard(a)
                    changed = True

        embeddings = _rdkit_embed(mol, n_confs=1)
        if embeddings:
            coords[:] = embeddings[0]
            # Optimize core atoms; branch atoms from embedding serve as
            # fixed anchors (though few branches touch the core directly).
            fixed_for_opt = {i: coords[i].copy() for i in range(n) if i not in core_atoms}
            _optimize_ring_system(
                mol,
                core_atoms,
                all_rings,
                bond_dihedral,
                coords,
                fixed_for_opt,
                parent,
            )

    # --- Phase 2: z-matrix place branch atoms outward from core ---
    placed: set[int] = set(core_atoms)
    # Only seed atom 0 when there are no ring atoms.  When rings exist,
    # atom 0 will be discovered by BFS from the ring so that branch
    # chains run outward from the ring — this lets the AMSR dihedrals
    # (which reference ring atoms) be consumed correctly.
    if not ring_atoms:
        placed.add(0)

    # Build outward parent tree: BFS from core, each branch atom's parent
    # is its neighbor closest to the core.  The AMSR dihedral for bond
    # (g, p) gives torsion(mi, g, p, mj) — when mj is the atom being
    # placed, this is a direct match for _choose_chain_dihedral.
    outward_parent: list[Optional[int]] = [None] * n
    branch_order: list[int] = []

    # Within the core, build a BFS tree so core atoms have parents.
    core_root = min(ring_atoms) if ring_atoms else (min(placed) if placed else 0)
    core_visited = {core_root}
    bfs = [core_root]
    qi = 0
    while qi < len(bfs):
        curr = bfs[qi]
        qi += 1
        for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
            b = nb.GetIdx()
            if b in placed and b not in core_visited:
                outward_parent[b] = curr
                core_visited.add(b)
                bfs.append(b)

    # BFS outward from core into branches
    for a in sorted(placed):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b not in placed and outward_parent[b] is None:
                outward_parent[b] = a
                branch_order.append(b)
                bfs.append(b)
    while qi < len(bfs):
        curr = bfs[qi]
        qi += 1
        for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
            b = nb.GetIdx()
            if b not in placed and outward_parent[b] is None:
                outward_parent[b] = curr
                branch_order.append(b)
                bfs.append(b)

    # Handle disconnected components not reachable from the core
    branch_set = set(branch_order)
    for i in range(n):
        if i in placed or i in branch_set:
            continue
        # Seed this disconnected component
        placed.add(i)
        for nb in mol.GetAtomWithIdx(i).GetNeighbors():
            b = nb.GetIdx()
            if b not in placed and b not in branch_set:
                outward_parent[b] = i
                branch_order.append(b)
                branch_set.add(b)
                bfs.append(b)
        while qi < len(bfs):
            curr = bfs[qi]
            qi += 1
            for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
                b = nb.GetIdx()
                if b not in placed and b not in branch_set:
                    outward_parent[b] = curr
                    branch_order.append(b)
                    branch_set.add(b)
                    bfs.append(b)

    # Initialize bookkeeping from core atom geometry
    for a in sorted(placed):
        p = outward_parent[a]
        if p is None or p not in placed:
            continue
        if first_child_idx[p] is None:
            first_child_idx[p] = a
            g = outward_parent[p]
            if g is not None and g in placed:
                gg = outward_parent[g]
                std_ref = (
                    coords[gg]
                    if (gg is not None and gg in placed)
                    else _synthetic_ref(coords, g, p)
                )
                first_child_torsion[p] = measure_torsion(std_ref, coords[g], coords[p], coords[a])
        child_count[p] += 1

    # Place branch atoms in BFS order (outward from core)
    for i in branch_order:
        alts = _place_chain_atom(
            mol,
            i,
            coords,
            outward_parent,
            child_count,
            first_child_torsion,
            first_child_idx,
            bond_dihedral,
        )
        placed.add(i)
        p = outward_parent[i]
        if alts and p is not None and _has_collision(mol, i, coords, placed, threshold=1.0):
            for alt in alts:
                child_count[p] -= 1
                _place_chain_atom(
                    mol,
                    i,
                    coords,
                    outward_parent,
                    child_count,
                    first_child_torsion,
                    first_child_idx,
                    bond_dihedral,
                    torsion_override=alt,
                )
                if not _has_collision(mol, i, coords, placed, threshold=1.0):
                    break

    # Build RDKit conformer
    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
