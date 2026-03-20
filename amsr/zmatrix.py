"""Z-matrix conformer generation from AMSR dihedrals.

Place ring atoms first (completing one ring before starting the next).
When a ring closes, fix chirality clashes and optimize torsions to close
the ring.  Then place non-ring atoms.

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


def _get_bond_angle(mol, a, b, c):
    """Ideal bond angle a-b-c in degrees."""
    hyb = mol.GetAtomWithIdx(b).GetHybridization()
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
    mol, i, p, g, gg, nth_child, first_child_torsion, in_ring, coords, bond_dihedral
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
        same_ring = False
        if gg is not None and in_ring:
            for ring in mol.GetRingInfo().AtomRings():
                if gg in ring and g in ring and p in ring and i in ring:
                    same_ring = True
                    break
        torsion = 0.0 if (in_ring and same_ring) else 180.0
        return torsion, None, [torsion + 180.0]

    # Subsequent children: offset from first child
    base = first_child_torsion[p] if first_child_torsion[p] is not None else 0.0
    if hyb_p == SP2:
        return base + 180.0, None, []
    if hyb_p == SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral == CW:
            return base + 120.0 * nth_child, None, []
        if chiral == CCW:
            return base - 120.0 * nth_child, None, []
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
# Ring closure
# ---------------------------------------------------------------------------


def _has_collision(mol, j, coords, placed, threshold=0.5):
    """Check if atom j collides with any placed non-bonded atom."""
    for other in placed:
        if other == j or mol.GetBondBetweenAtoms(j, other) is not None:
            continue
        if np.linalg.norm(coords[j] - coords[other]) < threshold:
            return True
    return False


def _has_bad_closure(mol, coords, placed, parent, threshold=1.5):
    """Check if any ring that just closed has a catastrophically bad closure bond."""
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ring_set = set(ring)
        if not ring_set.issubset(placed):
            continue
        for idx in range(len(ring)):
            a, b = ring[idx], ring[(idx + 1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond and parent[a] != b and parent[b] != a:
                dist = np.linalg.norm(coords[a] - coords[b])
                ideal = _get_bond_length(mol, a, b)
                if abs(dist - ideal) > threshold:
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
            _has_collision(mol, j, coords, placed) or _has_bad_closure(mol, coords, placed, parent)
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


def _prepare_adjustable(mol, new_atoms, coords, parent):
    """Precompute arrays for adjustable ring atoms (those with a grandparent).

    Returns (atom_indices, ref_indices, g_indices, p_indices, bond_lens,
    bond_angles, torsions) as numpy arrays for _replace_atoms().
    """
    atoms, refs, gs, ps, bls, bas, torsions = [], [], [], [], [], [], []
    for i in sorted(new_atoms):
        p = parent[i]
        if p is None or parent[p] is None:
            continue
        g = parent[p]
        ref_idx = parent[g] if parent[g] is not None else -1
        t = measure_torsion(_get_ref(coords, ref_idx, g, p), coords[g], coords[p], coords[i])
        atoms.append(i)
        refs.append(ref_idx)
        gs.append(g)
        ps.append(p)
        bls.append(_get_bond_length(mol, p, i))
        bas.append(_get_bond_angle(mol, g, p, i))
        torsions.append(t)
    return (
        np.array(atoms, dtype=np.intp),
        np.array(refs, dtype=np.intp),
        np.array(gs, dtype=np.intp),
        np.array(ps, dtype=np.intp),
        np.array(bls),
        np.array(bas),
        np.array(torsions),
    )


def _close_ring(mol, ring, new_atoms, coords, parent, placed):
    """Fix chirality clashes, then optimize torsions to close a ring."""
    from scipy.optimize import minimize

    ring_set = set(ring)

    # Precompute adjustable atom data
    atom_indices, ref_indices, g_indices, p_indices, bond_lens, bond_angles, torsions = (
        _prepare_adjustable(mol, new_atoms, coords, parent)
    )
    if len(atom_indices) == 0:
        return

    # Fix chirality clashes: if any new atom clashes with a placed atom,
    # flip the first adjustable SP3 atom's chirality (-240° torsion shift).
    has_clash = any(
        np.linalg.norm(coords[i] - coords[j]) < 0.8
        for i in atom_indices
        for j in placed
        if j not in ring_set and mol.GetBondBetweenAtoms(int(i), int(j)) is None
    )
    if has_clash and mol.GetAtomWithIdx(parent[atom_indices[0]]).GetHybridization() == SP3:
        torsions = torsions.copy()
        torsions[0] -= 240.0
        _replace_atoms(
            torsions,
            atom_indices,
            ref_indices,
            g_indices,
            p_indices,
            bond_lens,
            bond_angles,
            coords,
        )
        for k in range(1, len(atom_indices)):
            ref = _get_ref(coords, ref_indices[k], g_indices[k], p_indices[k])
            torsions[k] = measure_torsion(
                ref, coords[g_indices[k]], coords[p_indices[k]], coords[atom_indices[k]]
            )

    # Find closure bonds
    closures_pairs = []
    closures_ideals = []
    for idx in range(len(ring)):
        a, b = ring[idx], ring[(idx + 1) % len(ring)]
        if mol.GetBondBetweenAtoms(a, b) and parent[a] != b and parent[b] != a:
            closures_pairs.append((a, b))
            closures_ideals.append(_get_bond_length(mol, a, b))

    if not closures_pairs:
        return
    max_gap = max(
        abs(np.linalg.norm(coords[a] - coords[b]) - bl)
        for (a, b), bl in zip(closures_pairs, closures_ideals)
    )
    if max_gap < 0.1:
        return

    # Find new angles at closure bonds
    angle_triples = []
    angle_ideals = []
    for a, b in closures_pairs:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            c = nb.GetIdx()
            if c != b and c in ring_set:
                angle_triples.append((c, a, b))
                angle_ideals.append(_get_bond_angle(mol, c, a, b))
        for nb in mol.GetAtomWithIdx(b).GetNeighbors():
            c = nb.GetIdx()
            if c != a and c in ring_set:
                angle_triples.append((a, b, c))
                angle_ideals.append(_get_bond_angle(mol, a, b, c))

    # Convert to arrays for the cost function
    cp = np.array(closures_pairs, dtype=np.intp)
    ci = np.array(closures_ideals)
    at = (
        np.array(angle_triples, dtype=np.intp) if angle_triples else np.empty((0, 3), dtype=np.intp)
    )
    ai = np.array(angle_ideals) if angle_ideals else np.empty(0)

    init = torsions.copy()
    saved = coords[sorted(new_atoms)].copy()
    new_list = sorted(new_atoms)

    def cost(t):
        _replace_atoms(
            t, atom_indices, ref_indices, g_indices, p_indices, bond_lens, bond_angles, coords
        )
        c = 0.0
        for k in range(len(cp)):
            c += (np.linalg.norm(coords[cp[k, 0]] - coords[cp[k, 1]]) - ci[k]) ** 2
        c += 3e-4 * np.sum((t - init) ** 2)
        for k in range(len(at)):
            c += 3e-4 * (measure_angle(coords, at[k, 0], at[k, 1], at[k, 2]) - ai[k]) ** 2
        return c

    init_cost = cost(init)
    best = minimize(cost, init, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-10})

    # If stuck (no improvement), try perturbed starting points to escape saddle points
    if best.fun >= init_cost - 1e-8:
        for delta in [15.0, -15.0, 30.0, -30.0]:
            perturbed = init + delta
            r = minimize(
                cost, perturbed, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-10}
            )
            if r.fun < best.fun:
                best = r

    if best.fun < init_cost:
        _replace_atoms(
            best.x,
            atom_indices,
            ref_indices,
            g_indices,
            p_indices,
            bond_lens,
            bond_angles,
            coords,
        )
    else:
        for k, i in enumerate(new_list):
            coords[i] = saved[k]


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
    2. Optionally fix chirality clashes and optimize torsions to close rings.
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
    closed_rings: set[int] = set()
    closed_atoms: set[int] = set()

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

            # Close any rings that completed during placement
            if refine_rings:
                for ri_idx, ring in enumerate(all_rings):
                    if ri_idx not in closed_rings and set(ring).issubset(placed):
                        _close_ring(mol, ring, set(ring) - closed_atoms, coords, parent, placed)
                        closed_rings.add(ri_idx)
                        closed_atoms.update(ring)

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
