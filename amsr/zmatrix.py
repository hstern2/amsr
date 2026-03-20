"""Z-matrix conformer generation from AMSR dihedrals.

Place ring atoms first (completing one ring before starting the next),
with DFS backtracking on collisions and bad closure angles.
Optimize non-planar torsions and bond angles jointly across each ring
system to close all rings.  Then place non-ring atoms.

The code is structured so that geometry primitives and the optimization
cost function use only numpy arrays (no RDKit), making them suitable for
reimplementation in C.
"""

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


# ---------------------------------------------------------------------------
# Geometry primitives (pure numpy — C-portable)
# ---------------------------------------------------------------------------


def place_atom(A, B, C, d, theta_deg, omega_deg):
    """Place atom D given refs A, B, C, bond length d, angle B-C-D, torsion A-B-C-D."""
    theta = np.radians(theta_deg)
    omega = np.radians(omega_deg)
    BC = C - B
    bc = np.linalg.norm(BC)
    if bc > 1e-10:
        BC = BC / bc
    else:
        BC = np.array([1.0, 0.0, 0.0])
    n = np.cross(B - A, BC)
    nn = np.linalg.norm(n)
    if nn < 1e-10:
        perp = np.array([1.0, 0.0, 0.0]) if abs(BC[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        n = np.cross(BC, perp)
        n /= np.linalg.norm(n)
    else:
        n /= nn
    m = np.cross(n, BC)
    st = np.sin(theta)
    return C + d * (-np.cos(theta) * BC + st * np.cos(omega) * m + st * np.sin(omega) * n)


def measure_torsion(p0, p1, p2, p3):
    """Torsion angle (degrees) for four 3-D points."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    n1n, n2n = np.linalg.norm(n1), np.linalg.norm(n2)
    if n1n < 1e-10 or n2n < 1e-10:
        return 0.0
    n1, n2 = n1 / n1n, n2 / n2n
    return np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)), np.dot(n1, n2)))


def measure_angle(coords, a, b, c):
    """Angle a-b-c (degrees) from coordinates."""
    v1, v2 = coords[a] - coords[b], coords[c] - coords[b]
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
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
    bc_n = np.linalg.norm(bc)
    if bc_n > 1e-10:
        bc = bc / bc_n
    else:
        bc = np.array([1.0, 0.0, 0.0])
    perp = np.array([1.0, 0.0, 0.0]) if abs(bc[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp = np.cross(bc, perp)
    perp /= np.linalg.norm(perp)
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
        # Offset from mj — only for first child placed
        if nth_child == 0 and (np.any(coords[mj]) or mj == 0):
            ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
            actual_mj = measure_torsion(ref, coords[g], coords[p], coords[mj])
            if hyb_p == SP2:
                return actual_mj + 180.0, None, []
            chiral = mol.GetAtomWithIdx(p).GetChiralTag()
            if chiral == CCW:
                return actual_mj - 120.0, None, []
            if chiral == CW:
                return actual_mj + 120.0, None, []
            # Unspecified chirality — ambiguous
            return actual_mj + 120.0, None, [actual_mj - 120.0]

    # No AMSR — first child: search placed neighbors of g for a reference
    # atom that gives a known dihedral (same-ring → 0°).
    if nth_child == 0 and in_ring:
        for nb in mol.GetAtomWithIdx(g).GetNeighbors():
            gg_c = nb.GetIdx()
            if gg_c == p or not (np.any(coords[gg_c]) or gg_c == 0):
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
            return base + sign * 120.0 * nth_child, None, []
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
        if np.linalg.norm(coords[j] - coords[other]) < threshold:
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
                dist = np.linalg.norm(coords[a] - coords[b])
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
# Currently uses L-BFGS-B with numerical gradients.
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


def _prepare_adjustable(mol, atoms, coords, parent):
    """Precompute z-matrix arrays for all atoms that have a grandparent.

    Returns (atom_indices, ref_indices, g_indices, p_indices, bond_lens,
    bond_angles, torsions, planar_mask) as numpy arrays.

    planar_mask[k] is True if atom k's torsion is locked (all four reference
    atoms lie in the same fully-SP2 ring).  The optimizer should adjust
    torsions only for non-planar atoms, but may adjust bond angles for all.
    """
    ri = mol.GetRingInfo()
    sp2_rings = [
        set(ring)
        for ring in ri.AtomRings()
        if all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring)
    ]

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

        check = {i, p, g}
        if ref_idx >= 0:
            check.add(ref_idx)
        is_planar = any(check.issubset(sr) for sr in sp2_rings)

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


def _collect_closure_constraints(mol, closure_pairs, system_set, bond_dihedral):
    """Collect angle and dihedral constraints for closure optimization.

    Returns (angle_triples, angle_ideals, dihedral_quads, dihedral_targets).
    """
    closure_set = {(min(a, b), max(a, b)) for a, b in closure_pairs}

    # Bond angle constraints at closure points
    angle_triples = []
    angle_ideals = []
    for a, b in closure_pairs:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            c = nb.GetIdx()
            if c != b and c in system_set:
                angle_triples.append((c, a, b))
                angle_ideals.append(_get_bond_angle(mol, c, a, b))
        for nb in mol.GetAtomWithIdx(b).GetNeighbors():
            c = nb.GetIdx()
            if c != a and c in system_set:
                angle_triples.append((a, b, c))
                angle_ideals.append(_get_bond_angle(mol, a, b, c))

    # AMSR dihedral constraints for bonds within the ring system
    dihedral_quads = []
    dihedral_targets = []
    for (i, j), (mi, mj, angle) in bond_dihedral.items():
        if i >= j:
            continue
        key = (min(i, j), max(i, j))
        if (key in closure_set or (i in system_set and j in system_set)) and (
            mi in system_set and mj in system_set
        ):
            dihedral_quads.append((mi, i, j, mj))
            dihedral_targets.append(float(angle))

    return angle_triples, angle_ideals, dihedral_quads, dihedral_targets


def _closure_residuals(
    coords,
    closure_pairs,
    closure_ideals,
    angle_triples,
    angle_ideals,
    dihedral_quads,
    dihedral_targets,
    torsions,
    init_torsions,
):
    """Compute residual vector for ring closure.

    Returns array of weighted residuals whose sum-of-squares is the cost.
    Each residual is independent — suitable for Gauss-Newton/LM solvers.

    Components:
      - Closure bond gap (weight 1.0)
      - Torsion regularization toward initial values (weight sqrt(3e-4))
      - Bond angle deviation at closure points (weight sqrt(3e-4))
      - AMSR dihedral deviation (weight sqrt(3e-3))
    """
    residuals = []

    # Closure bond gaps (weight 1.0)
    for k in range(len(closure_pairs)):
        a, b = closure_pairs[k]
        gap = np.linalg.norm(coords[a] - coords[b]) - closure_ideals[k]
        residuals.append(gap)

    # Torsion regularization (weight sqrt(3e-4) ≈ 0.0173)
    w_torsion = np.sqrt(3e-4)
    for k in range(len(torsions)):
        residuals.append(w_torsion * (torsions[k] - init_torsions[k]))

    # Angle constraints (weight sqrt(3e-4))
    w_angle = np.sqrt(3e-4)
    for k in range(len(angle_triples)):
        a, b, c = angle_triples[k]
        residuals.append(w_angle * (measure_angle(coords, a, b, c) - angle_ideals[k]))

    # AMSR dihedral constraints (weight sqrt(3e-3) ≈ 0.0548)
    w_dihedral = np.sqrt(3e-3)
    for k in range(len(dihedral_quads)):
        mi, i, j, mj = dihedral_quads[k]
        measured = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
        diff = (measured - dihedral_targets[k] + 180.0) % 360.0 - 180.0
        residuals.append(w_dihedral * diff)

    return np.array(residuals)


def _close_ring_system(mol, system_atoms, all_rings, coords, parent, bond_dihedral):
    """Optimize non-planar torsions and bond angles to close all rings.

    Planar (aromatic) ring torsions are kept fixed.  Bond angles are
    adjustable with regularization toward ideal values.
    """
    from scipy.optimize import minimize

    closure_pairs, closure_ideals = _find_closure_bonds(mol, system_atoms, all_rings, parent)
    if len(closure_pairs) < 2:
        return

    system_set = set(system_atoms)
    (
        atom_indices,
        ref_indices,
        g_indices,
        p_indices,
        bond_lens,
        bond_angles,
        torsions,
        planar_mask,
    ) = _prepare_adjustable(mol, system_atoms, coords, parent)
    if len(atom_indices) == 0:
        return

    angle_triples, angle_ideals, dihedral_quads, dihedral_targets = _collect_closure_constraints(
        mol, closure_pairs, system_set, bond_dihedral
    )

    cp = np.array(closure_pairs, dtype=np.intp)
    ci = np.array(closure_ideals)
    at = (
        np.array(angle_triples, dtype=np.intp) if angle_triples else np.empty((0, 3), dtype=np.intp)
    )
    ai = np.array(angle_ideals) if angle_ideals else np.empty(0)
    dq = (
        np.array(dihedral_quads, dtype=np.intp)
        if dihedral_quads
        else np.empty((0, 4), dtype=np.intp)
    )
    dt = np.array(dihedral_targets) if dihedral_targets else np.empty(0)

    free_tor_idx = np.where(~planar_mask)[0]
    n_free_tor = len(free_tor_idx)

    init_free_tor = torsions[free_tor_idx].copy()
    init_angles = bond_angles.copy()
    init_x = np.concatenate([init_free_tor, init_angles])

    sys_list = sorted(system_atoms)
    saved = coords[sys_list].copy()

    def cost(x):
        free_t = x[:n_free_tor]
        angles = x[n_free_tor:]
        full_torsions = torsions.copy()
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
        r = _closure_residuals(coords, cp, ci, at, ai, dq, dt, full_torsions, torsions)
        # Angle regularization toward ideal values
        w_ang_reg = np.sqrt(3e-4)
        ang_resid = w_ang_reg * (angles - init_angles)
        return np.dot(r, r) + np.dot(ang_resid, ang_resid)

    # Multi-start optimization: try initial values + perturbations
    init_cost = cost(init_x)
    best = minimize(cost, init_x, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-10})

    if best.fun >= init_cost - 1e-8:
        for delta in [15.0, -15.0, 30.0, -30.0]:
            x0 = init_x.copy()
            x0[:n_free_tor] += delta  # perturb torsions only
            r = minimize(cost, x0, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-10})
            if r.fun < best.fun:
                best = r

    if best.fun < init_cost:
        cost(best.x)  # apply the best solution to coords
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

            # Optimize closure bonds jointly across the ring system
            if refine_rings:
                _close_ring_system(mol, ring_systems[si], all_rings, coords, parent, bond_dihedral)

            placed_systems.add(si)
        else:
            _place_one(
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

    # Build RDKit conformer
    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
