"""Z-matrix conformer generation from AMSR dihedrals.

Place ring atoms first (completing one ring before starting the next),
with DFS backtracking on collisions and bad closure angles.
Optimize non-planar torsions and bond angles jointly across each ring
system to close all rings.  Then place non-ring atoms.

The code is structured so that geometry primitives and the optimization
cost function use only numpy arrays (no RDKit), making them suitable for
reimplementation in C.
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


def _ring_visit_order(mol, system, parent):
    """Order atoms within a ring system: complete one ring before starting next."""
    ri = mol.GetRingInfo()
    ordered = []
    remaining = set(system)
    placed: set[int] = set()

    while remaining:
        ready = [
            a
            for a in remaining
            if parent[a] is None or parent[a] not in remaining or parent[a] in placed
        ]
        if not ready:
            ready = [min(remaining)]
        if ordered:
            last = ordered[-1]
            same = [a for a in ready if ri.AreAtomsInSameRing(a, last)]
            if same:
                ready = same
        pick = min(ready)
        ordered.append(pick)
        placed.add(pick)
        remaining.discard(pick)

    return ordered


# ---------------------------------------------------------------------------
# Dihedral selection
# ---------------------------------------------------------------------------


def _choose_dihedral(
    mol,
    i,
    p,
    g,
    gg,
    nth_child,
    first_child_torsion,
    first_child_idx,
    in_ring,
    coords,
    bond_dihedral,
):
    """Choose the dihedral angle for placing atom i from parent p.

    Returns (torsion, ref_override, alternatives) where alternatives is a
    list of other torsion values to try if the first causes a collision.
    """
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    # Check for AMSR dihedral on backward bond (g, p)
    if (g, p) in bond_dihedral:
        mi, mj, angle = bond_dihedral[(g, p)]
        if mj == i:
            return angle, mi if mi != gg else None, []
        # Offset from mj based on graph-order position of i relative to mj
        if nth_child == 0 and _is_placed(coords, mj):
            if hyb_p == SP2:
                ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                actual_mj = measure_torsion(ref, coords[g], coords[p], coords[mj])
                return actual_mj + 180.0, None, []
            torsion, alts = _sp3_offset(mol, p, g, gg, i, mj, coords)
            return torsion, None, alts

    # No AMSR — first child: search placed neighbors of g for a reference
    # atom that gives a known dihedral (same-ring → 0°).
    if nth_child == 0 and in_ring:
        for nb in mol.GetAtomWithIdx(g).GetNeighbors():
            gg_c = nb.GetIdx()
            if gg_c == p or not _is_placed(coords, gg_c):
                continue
            for ring in mol.GetRingInfo().AtomRings():
                if gg_c in ring and g in ring and p in ring and i in ring:
                    ref_ovr = gg_c if gg_c != gg else None
                    return 0.0, ref_ovr, []

    # Fallback: default torsion (ambiguous)
    if nth_child == 0:
        coplanar = False
        if in_ring:
            if gg is not None:
                for ring in mol.GetRingInfo().AtomRings():
                    if gg in ring and g in ring and p in ring and i in ring:
                        coplanar = True
                        break
            # Fused aromatic ring junctions: gg may not share a ring with
            # i, but the bond is still planar.  Only applies when gg is
            # in a fully SP2 (aromatic) ring.  Use 0° unless gg's ring
            # has <6 members (tighter interior angles flip the torsion).
            if (
                not coplanar
                and gg is not None
                and hyb_p == SP2
                and mol.GetAtomWithIdx(g).GetHybridization() == SP2
                and mol.GetAtomWithIdx(i).GetHybridization() == SP2
                and any(
                    gg in ring
                    and all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring)
                    for ring in mol.GetRingInfo().AtomRings()
                )
            ):
                coplanar = True
                for ring in mol.GetRingInfo().AtomRings():
                    if gg in ring and g in ring and i not in ring:
                        if len(ring) < 6:
                            coplanar = False
                        break
        # SP3 parent with already-placed neighbor: offset using chirality.
        if not coplanar and hyb_p == SP3:
            for nb in mol.GetAtomWithIdx(p).GetNeighbors():
                k = nb.GetIdx()
                if k != g and k != i and _is_placed(coords, k):
                    torsion, alts = _sp3_offset(mol, p, g, gg, i, k, coords)
                    return torsion, None, alts

        torsion = 0.0 if coplanar else 180.0
        return torsion, None, [torsion + 180.0]

    # Subsequent children: offset from first child
    base = first_child_torsion[p] if first_child_torsion[p] is not None else 0.0
    if hyb_p == SP2:
        return base + 180.0, None, []
    if hyb_p == SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral in (CW, CCW):
            # CW/CCW is defined relative to graph neighbor order.  If
            # ring-system placement caused the first child to be placed
            # out of graph order, flip the sign.
            fc = first_child_idx[p]
            children = [
                nb.GetIdx() for nb in mol.GetAtomWithIdx(p).GetNeighbors() if nb.GetIdx() != g
            ]
            swapped = (
                len(children) >= 2 and fc is not None and children.index(fc) > children.index(i)
            )
            if chiral == CW:
                sign = -1 if swapped else 1
            else:
                sign = 1 if swapped else -1
            return base + sign * 120.0 * nth_child, None, [base - sign * 120.0 * nth_child]
        sign = -1 if base > 0 else 1
        return base + sign * 120.0 * nth_child, None, [base - sign * 120.0 * nth_child]
    return base + 180.0, None, []


# ---------------------------------------------------------------------------
# Single-atom placement
# ---------------------------------------------------------------------------


def _place_one(
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
    """Place atom i using z-matrix from its parent chain.

    If torsion_override is given, use it instead of _choose_dihedral.
    Returns list of alternative torsions (empty if placement is unambiguous).
    """
    p = parent[i]
    if p is None:
        return []

    g = parent[p]
    bond_len = _get_bond_length(mol, p, i)
    bp = mol.GetBondBetweenAtoms(p, i)
    in_ring = bp is not None and bp.IsInRing()
    alternatives: list[float] = []

    if g is None:
        # No grandparent — special cases for first/second child
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
            torsion, ref_override, alternatives = _choose_dihedral(
                mol,
                i,
                p,
                g,
                gg,
                child_count[p],
                first_child_torsion,
                first_child_idx,
                in_ring,
                coords,
                bond_dihedral,
            )
        if ref_override is not None:
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


# ---------------------------------------------------------------------------
# Ring-system placement (DFS backtracking)
# ---------------------------------------------------------------------------


def _has_collision(mol, j, coords, placed, threshold=0.5):
    """Check if atom j collides with any placed non-bonded atom."""
    for other in placed:
        if other == j or mol.GetBondBetweenAtoms(j, other) is not None:
            continue
        if _norm3(coords[j] - coords[other]) < threshold:
            return True
    return False


def _has_bad_closure(
    mol, coords, placed, parent, just_placed=None, dist_threshold=1.5, angle_threshold=80.0
):
    """Check if a just-completed ring has a bad closure bond or angle.

    Only examines closure bonds involving `just_placed` to avoid
    false positives from previously placed rings.
    """
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ring_set = set(ring)
        if not ring_set.issubset(placed):
            continue
        if just_placed is not None and just_placed not in ring_set:
            continue
        for idx in range(len(ring)):
            a, b = ring[idx], ring[(idx + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond and parent[a] != b and parent[b] != a:
                if just_placed is not None and just_placed != a and just_placed != b:
                    continue
                dist = _norm3(coords[a] - coords[b])
                ideal = _get_bond_length(mol, a, b)
                if abs(dist - ideal) > dist_threshold:
                    return True
                for endpoint, other in [(a, b), (b, a)]:
                    for nb in mol.GetAtomWithIdx(endpoint).GetNeighbors():
                        c = nb.GetIdx()
                        if c == other or c not in placed:
                            continue
                        ideal_ang = _get_bond_angle(mol, c, endpoint, other)
                        actual_ang = measure_angle(coords, c, endpoint, other)
                        if abs(actual_ang - ideal_ang) > angle_threshold:
                            return True
    return False


def _save_state(coords, child_count, first_child_torsion, first_child_idx, placed):
    """Snapshot mutable placement state for backtracking."""
    return (
        coords.copy(),
        child_count[:],
        first_child_torsion[:],
        first_child_idx[:],
        placed.copy(),
    )


def _restore_state(snap, coords, child_count, first_child_torsion, first_child_idx, placed):
    """Restore mutable placement state from snapshot."""
    coords[:] = snap[0]
    child_count[:] = list(snap[1])
    first_child_torsion[:] = list(snap[2])
    first_child_idx[:] = list(snap[3])
    placed.clear()
    placed.update(snap[4])


def _place_all_default(
    mol,
    visit_order,
    coords,
    parent,
    child_count,
    first_child_torsion,
    first_child_idx,
    bond_dihedral,
    placed,
):
    """Place all unplaced atoms in visit_order with default torsions (no backtracking)."""
    for j in visit_order:
        if j not in placed and parent[j] is not None:
            _place_one(
                mol,
                j,
                coords,
                parent,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
            )
            placed.add(j)


def _place_ring_system_dfs(
    mol,
    visit_order,
    coords,
    parent,
    child_count,
    first_child_torsion,
    first_child_idx,
    bond_dihedral,
    placed,
):
    """Place ring system atoms with backtracking on collisions/bad closures."""
    clean = _save_state(coords, child_count, first_child_torsion, first_child_idx, placed)
    # Stack: (atom_idx, remaining_alternatives, pre-placement_snapshot)
    stack: list[tuple[int, list[float], tuple]] = []
    k = 0

    while k < len(visit_order):
        j = visit_order[k]
        if j in placed or parent[j] is None:
            k += 1
            continue

        snap = _save_state(coords, child_count, first_child_torsion, first_child_idx, placed)
        alts = _place_one(
            mol, j, coords, parent, child_count, first_child_torsion, first_child_idx, bond_dihedral
        )
        placed.add(j)
        stack.append((j, alts, snap))

        if not (
            _has_collision(mol, j, coords, placed)
            or _has_bad_closure(mol, coords, placed, parent, just_placed=j)
        ):
            k += 1
            continue

        # Backtrack to most recent atom with untried alternatives
        resolved = False
        while stack and not resolved:
            bj, balts, bsnap = stack[-1]
            if not balts:
                stack.pop()
                continue
            alt = balts.pop(0)
            _restore_state(bsnap, coords, child_count, first_child_torsion, first_child_idx, placed)
            _place_one(
                mol,
                bj,
                coords,
                parent,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
                torsion_override=alt,
            )
            placed.add(bj)
            # Re-check: if alternative still collides, keep backtracking.
            if _has_collision(mol, bj, coords, placed):
                continue
            stack[-1] = (bj, balts, bsnap)
            k = visit_order.index(bj) + 1
            resolved = True

        if not resolved:
            _restore_state(clean, coords, child_count, first_child_torsion, first_child_idx, placed)
            _place_all_default(
                mol,
                visit_order,
                coords,
                parent,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
                placed,
            )
            return


# ---------------------------------------------------------------------------
# Ring-system closure optimization
#
# After z-matrix placement, closure bonds (ring bonds that aren't in the
# parent chain) may have gaps.  We optimize non-planar torsions and bond
# angles jointly to close every ring simultaneously.  Planar (aromatic)
# ring torsions are kept fixed to avoid distorting flat rings.
#
# The cost function is a sum of squared residuals — natural for future
# gradient-based (analytical Jacobian) or Gauss-Newton/LM solvers.
# Currently uses L-BFGS-B with numerical gradients.  A single optimization
# from the initial z-matrix placement is sufficient.
# ---------------------------------------------------------------------------


def _get_ref(coords, ref_idx, g, p):
    """Get reference point: coords[ref_idx] if valid, else synthetic."""
    if ref_idx >= 0:
        return coords[ref_idx]
    return _synthetic_ref(coords, g, p)


def _replace_atoms(
    torsions, atom_indices, ref_indices, g_indices, p_indices, bond_lens, bond_angles, coords
):
    """Re-place atoms given new torsion angles. Handles synthetic refs (ref_index == -1)."""
    for k in range(len(atom_indices)):
        coords[atom_indices[k]] = place_atom(
            _get_ref(coords, ref_indices[k], g_indices[k], p_indices[k]),
            coords[g_indices[k]],
            coords[p_indices[k]],
            bond_lens[k],
            bond_angles[k],
            torsions[k],
        )


def _prepare_adjustable(mol, atoms, coords, parent, sp2_rings):
    """Precompute z-matrix arrays for all atoms that have a grandparent.

    Returns (atom_indices, ref_indices, g_indices, p_indices, bond_lens,
    bond_angles, torsions, planar_mask) as numpy arrays.

    planar_mask[k] is True if atom k's torsion is locked (both atom and
    parent lie in the same fully-SP2 ring).  The optimizer should adjust
    torsions only for non-planar atoms, but may adjust bond angles for all.
    """
    a_idx, r_idx, g_idx, p_idx, bls, bas, tors, planar = (
        [],
        [],
        [],
        [],
        [],
        [],
        [],
        [],
    )
    for i in sorted(atoms):
        p = parent[i]
        if p is None or parent[p] is None:
            continue
        g = parent[p]
        ref_idx = parent[g] if parent[g] is not None else -1

        is_planar = any({i, p}.issubset(sr) for sr in sp2_rings)

        t = measure_torsion(_get_ref(coords, ref_idx, g, p), coords[g], coords[p], coords[i])
        a_idx.append(i)
        r_idx.append(ref_idx)
        g_idx.append(g)
        p_idx.append(p)
        bls.append(_get_bond_length(mol, p, i))
        bas.append(_get_bond_angle(mol, g, p, i))
        tors.append(t)
        planar.append(is_planar)
    return (
        np.array(a_idx, dtype=np.intp),
        np.array(r_idx, dtype=np.intp),
        np.array(g_idx, dtype=np.intp),
        np.array(p_idx, dtype=np.intp),
        np.array(bls),
        np.array(bas),
        np.array(tors),
        np.array(planar, dtype=bool),
    )


def _find_closure_bonds(mol, system_atoms, all_rings, parent):
    """Find all closure bonds (non-parent ring bonds) in a ring system.

    Returns list of (atom_a, atom_b) pairs and their ideal bond lengths.
    """
    seen = set()
    pairs = []
    ideals = []
    system_set = set(system_atoms)
    for ring in all_rings:
        if not set(ring).issubset(system_set):
            continue
        for idx in range(len(ring)):
            a, b = ring[idx], ring[(idx + 1) % len(ring)]
            key = (min(a, b), max(a, b))
            if (
                key not in seen
                and mol.GetBondBetweenAtoms(a, b)
                and parent[a] != b
                and parent[b] != a
            ):
                seen.add(key)
                pairs.append((a, b))
                ideals.append(_get_bond_length(mol, a, b))
    return pairs, ideals


def _collect_closure_constraints(mol, closure_pairs, system_set, bond_dihedral, placed, sp2_rings):
    """Collect angle and dihedral constraints for closure optimization.

    Returns (angle_triples, angle_ideals, dihedral_quads, dihedral_targets,
             chiral_quads, chiral_targets).
    """
    closure_set = {(min(a, b), max(a, b)) for a, b in closure_pairs}

    # Bond angle constraints at closure points.  Only include neighbors
    # that are already placed (have valid coordinates).
    angle_triples = []
    angle_ideals = []
    for a, b in closure_pairs:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            c = nb.GetIdx()
            if c != b and c in placed:
                angle_triples.append((c, a, b))
                angle_ideals.append(_get_bond_angle(mol, c, a, b))
        for nb in mol.GetAtomWithIdx(b).GetNeighbors():
            c = nb.GetIdx()
            if c != a and c in placed:
                angle_triples.append((a, b, c))
                angle_ideals.append(_get_bond_angle(mol, a, b, c))

    # AMSR dihedral constraints.
    dihedral_quads = []
    dihedral_targets = []
    for (i, j), (mi, mj, angle) in bond_dihedral.items():
        if i >= j:
            continue
        key = (min(i, j), max(i, j))
        if (
            mi in placed
            and mj in placed
            and (key in closure_set or (i in system_set and j in system_set))
        ):
            dihedral_quads.append((mi, i, j, mj))
            dihedral_targets.append(float(angle))

    # SP2 planarity: improper torsion n1-center-n2-n3 = ±180° for SP2
    # ring atoms with 3 placed neighbors.  Only apply when all neighbors
    # share at least one all-SP2 ring with the center — at ring junctions
    # bridging non-planar rings, the atom may be slightly pyramidal.
    for a in system_set:
        atom = mol.GetAtomWithIdx(a)
        if atom.GetHybridization() != SP2:
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
        if len(nbrs) != 3:
            continue
        if not all(n in placed for n in nbrs):
            continue
        if not any({a} | set(nbrs) <= sr for sr in sp2_rings):
            continue
        dihedral_quads.append((nbrs[0], a, nbrs[1], nbrs[2]))
        dihedral_targets.append(180.0)

    # SP3 chirality: improper torsion n0-center-n1-n2 = ±120° for chiral
    # SP3 ring atoms with 3 placed neighbors.  Returned separately so the
    # optimizer can weight them independently.
    chiral_quads = []
    chiral_targets = []
    for a in system_set:
        atom = mol.GetAtomWithIdx(a)
        chiral = atom.GetChiralTag()
        if chiral not in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
        if len(nbrs) != 3:
            continue
        if not all(n in placed for n in nbrs):
            continue
        target = 120.0 if chiral == CW else -120.0
        chiral_quads.append((nbrs[0], a, nbrs[1], nbrs[2]))
        chiral_targets.append(target)

    return (
        angle_triples,
        angle_ideals,
        dihedral_quads,
        dihedral_targets,
        chiral_quads,
        chiral_targets,
    )


_W_REG = math.sqrt(3e-4)
_W_DIHEDRAL = math.sqrt(3e-3)
_W_CHIRAL = math.sqrt(3e-2)


def _closure_residuals(
    coords,
    closure_pairs,
    closure_ideals,
    angle_triples,
    angle_ideals,
    dihedral_quads,
    dihedral_targets,
    chiral_quads,
    chiral_targets,
    torsions,
    init_torsions,
    bond_angles,
    init_angles,
    residuals,
):
    """Compute residual vector for ring closure (writes into pre-allocated residuals array).

    Components (weights):
      - Closure bond gap (1.0)
      - Torsion regularization toward initial values (sqrt(3e-4))
      - Bond angle regularization toward ideal values (sqrt(3e-4))
      - Bond angle deviation at closure points (sqrt(3e-4))
      - AMSR dihedral deviation (sqrt(3e-3))
      - SP3 chirality improper torsion (sqrt(3e-2))
    """
    off = 0

    # Closure bond gaps (weight 1.0)
    n_cp = len(closure_pairs)
    diffs = coords[closure_pairs[:, 0]] - coords[closure_pairs[:, 1]]
    residuals[off : off + n_cp] = _batch_norm3(diffs) - closure_ideals
    off += n_cp

    # Torsion regularization
    n_tor = len(torsions)
    residuals[off : off + n_tor] = _W_REG * (torsions - init_torsions)
    off += n_tor

    # Bond angle regularization toward ideal values
    n_ang = len(bond_angles)
    residuals[off : off + n_ang] = _W_REG * (bond_angles - init_angles)
    off += n_ang

    # Angle constraints at closure points
    n_at = len(angle_triples)
    if n_at > 0:
        residuals[off : off + n_at] = _W_REG * (
            _batch_measure_angle(coords, angle_triples) - angle_ideals
        )
    off += n_at

    # AMSR dihedral + SP2 planarity constraints
    n_dq = len(dihedral_quads)
    if n_dq > 0:
        measured = _batch_measure_torsion(
            coords[dihedral_quads[:, 0]],
            coords[dihedral_quads[:, 1]],
            coords[dihedral_quads[:, 2]],
            coords[dihedral_quads[:, 3]],
        )
        residuals[off : off + n_dq] = _W_DIHEDRAL * (
            (measured - dihedral_targets + 180.0) % 360.0 - 180.0
        )
    off += n_dq

    # SP3 chirality constraints
    n_cq = len(chiral_quads)
    if n_cq > 0:
        measured = _batch_measure_torsion(
            coords[chiral_quads[:, 0]],
            coords[chiral_quads[:, 1]],
            coords[chiral_quads[:, 2]],
            coords[chiral_quads[:, 3]],
        )
        residuals[off : off + n_cq] = _W_CHIRAL * (
            (measured - chiral_targets + 180.0) % 360.0 - 180.0
        )
    off += n_cq

    return residuals


def _close_ring_system(mol, system_atoms, all_rings, coords, parent, bond_dihedral, placed):
    """Optimize non-planar torsions and bond angles to close all rings.

    Planar (aromatic) ring torsions are kept fixed.  Bond angles are
    adjustable with regularization toward ideal values.
    """
    from scipy.optimize import least_squares

    closure_pairs, closure_ideals = _find_closure_bonds(mol, system_atoms, all_rings, parent)
    if len(closure_pairs) < 2:
        return

    system_set = set(system_atoms)
    ri = mol.GetRingInfo()
    sp2_rings = [
        set(ring)
        for ring in ri.AtomRings()
        if all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring)
    ]

    (
        atom_indices,
        ref_indices,
        g_indices,
        p_indices,
        bond_lens,
        bond_angles,
        torsions,
        planar_mask,
    ) = _prepare_adjustable(mol, system_atoms, coords, parent, sp2_rings)
    if len(atom_indices) == 0:
        return

    (
        angle_triples,
        angle_ideals,
        dihedral_quads,
        dihedral_targets,
        chiral_quads,
        chiral_targets,
    ) = _collect_closure_constraints(
        mol, closure_pairs, system_set, bond_dihedral, placed, sp2_rings
    )

    def _to_intp(lst, cols):
        return (
            np.array(lst, dtype=np.intp).reshape(-1, cols)
            if lst
            else np.empty((0, cols), dtype=np.intp)
        )

    cp = np.array(closure_pairs, dtype=np.intp).reshape(-1, 2)
    ci = np.array(closure_ideals)
    at = _to_intp(angle_triples, 3)
    ai = np.array(angle_ideals) if angle_ideals else np.empty(0)
    dq = _to_intp(dihedral_quads, 4)
    dt = np.array(dihedral_targets) if dihedral_targets else np.empty(0)
    cq = _to_intp(chiral_quads, 4)
    ct = np.array(chiral_targets) if chiral_targets else np.empty(0)

    free_tor_idx = np.where(~planar_mask)[0]
    n_free_tor = len(free_tor_idx)

    init_free_tor = torsions[free_tor_idx].copy()
    init_angles = bond_angles.copy()
    init_x = np.concatenate([init_free_tor, init_angles])

    # Pre-allocate residual and torsion buffers
    n_residuals = len(cp) + len(torsions) + len(bond_angles) + len(at) + len(dq) + len(cq)
    residuals_buf = np.empty(n_residuals)
    full_torsions = torsions.copy()

    sys_list = sorted(system_atoms)
    saved = coords[sys_list].copy()

    def residual_fn(x):
        free_t = x[:n_free_tor]
        angles = x[n_free_tor:]
        full_torsions[:] = torsions
        full_torsions[free_tor_idx] = free_t
        _replace_atoms(
            full_torsions,
            atom_indices,
            ref_indices,
            g_indices,
            p_indices,
            bond_lens,
            angles,
            coords,
        )
        _closure_residuals(
            coords,
            cp,
            ci,
            at,
            ai,
            dq,
            dt,
            cq,
            ct,
            full_torsions,
            torsions,
            angles,
            init_angles,
            residuals_buf,
        )
        return residuals_buf

    def residual_fn_copy(x):
        """Return a copy for least_squares (which retains references)."""
        return residual_fn(x).copy()

    init_r = residual_fn(init_x)
    init_cost = np.dot(init_r, init_r)
    best = least_squares(residual_fn_copy, init_x, method="lm", ftol=1e-10, xtol=1e-10, gtol=1e-10)

    if best.cost * 2.0 < init_cost:
        residual_fn(best.x)  # apply the best solution to coords
    else:
        coords[sys_list] = saved


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
    refine_rings: bool = True,
) -> Chem.Mol:
    """Generate 3D conformer by z-matrix atom-by-atom placement.

    1. Place ring atoms first (completing one ring before starting the next).
    2. Optionally optimize torsions jointly to close rings.
    3. Place non-ring atoms.
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

    # Ring systems (sorted by lowest atom) and individual rings
    ring_systems = _find_ring_systems(mol)
    ring_systems.sort(key=lambda s: min(s))
    atom_to_system: dict[int, int] = {}
    for si, sys in enumerate(ring_systems):
        for a in sys:
            atom_to_system[a] = si

    all_rings = [tuple(r) for r in mol.GetRingInfo().AtomRings()]

    placed = {0}
    placed_systems: set[int] = set()

    for i in range(1, n):
        if parent[i] is None:
            continue

        if i in atom_to_system:
            si = atom_to_system[i]
            if si in placed_systems:
                continue

            # Place all atoms with DFS backtracking on collisions
            visit = _ring_visit_order(mol, ring_systems[si], parent)
            _place_ring_system_dfs(
                mol,
                visit,
                coords,
                parent,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
                placed,
            )

            # Temporarily place non-ring children of ring atoms so that
            # closure constraints can use their angles (important for
            # bridged ring systems like dibenzazepines).
            temp_placed: list[int] = []
            for ra in ring_systems[si]:
                for nb in mol.GetAtomWithIdx(ra).GetNeighbors():
                    ci = nb.GetIdx()
                    if ci not in placed and ci not in atom_to_system:
                        _place_one(
                            mol,
                            ci,
                            coords,
                            parent,
                            child_count,
                            first_child_torsion,
                            first_child_idx,
                            bond_dihedral,
                        )
                        placed.add(ci)
                        temp_placed.append(ci)

            # Optimize closure bonds jointly across the ring system
            if refine_rings:
                _close_ring_system(
                    mol, ring_systems[si], all_rings, coords, parent, bond_dihedral, placed
                )

            # Undo temporary placements so non-ring atoms are re-placed
            # with correct coordinates after ring refinement.
            for ci in temp_placed:
                placed.discard(ci)
                pi = parent[ci]
                if pi is not None:
                    child_count[pi] -= 1
                    if first_child_idx[pi] == ci:
                        first_child_idx[pi] = None
                        first_child_torsion[pi] = None

            # Refresh first_child_torsion — closure optimization moved atoms.
            for p in ring_systems[si]:
                fc = first_child_idx[p]
                if fc is None:
                    continue
                g = parent[p]
                if g is None:
                    continue
                gg = parent[g]
                std_ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                first_child_torsion[p] = measure_torsion(std_ref, coords[g], coords[p], coords[fc])

            placed_systems.add(si)
        else:
            alts = _place_one(
                mol,
                i,
                coords,
                parent,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
            )
            placed.add(i)
            p = parent[i]
            if alts and p is not None and _has_collision(mol, i, coords, placed, threshold=1.0):
                for alt in alts:
                    child_count[p] -= 1
                    _place_one(
                        mol,
                        i,
                        coords,
                        parent,
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
