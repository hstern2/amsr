"""Z-matrix conformer generation from AMSR dihedrals.

Place ring atoms first (completing one ring before starting the next).
When a ring closes, adjust torsions/angles of its new atoms to minimize
the closure bond gap.  Then place non-ring atoms.
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
    """Place atom D given three reference points A, B, C."""
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
    """Choose the dihedral angle for placing atom i from parent p."""
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    # Check for direct AMSR dihedral match (mj == i) regardless of nth_child
    if (g, p) in bond_dihedral:
        mi, mj, angle = bond_dihedral[(g, p)]
        if mj == i:
            return angle, mi if mi != gg else None
        # Offset from mj only for the first child placed
        if nth_child == 0 and (np.any(coords[mj]) or mj == 0):
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
            if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                return actual_mj + 120.0, None
            return actual_mj + 120.0, None

    if nth_child == 0:
        same_ring = False
        if gg is not None and in_ring:
            for ring in mol.GetRingInfo().AtomRings():
                if gg in ring and g in ring and p in ring and i in ring:
                    same_ring = True
                    break
        return _default_torsion(hyb_p, in_ring, same_ring), None

    base = first_child_torsion[p] if first_child_torsion[p] is not None else 0.0
    if hyb_p == Chem.HybridizationType.SP2:
        return base + 180.0, None
    if hyb_p == Chem.HybridizationType.SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return base + 120.0 * nth_child, None
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return base - 120.0 * nth_child, None
        if base > 0:
            return base - 120.0 * nth_child, None
        return base + 120.0 * nth_child, None
    return base + 180.0, None


def _default_torsion(hyb, in_ring, same_ring):
    if in_ring and same_ring:
        return 0.0
    if hyb == Chem.HybridizationType.SP2 and in_ring:
        return 180.0
    return 180.0


# ---------------------------------------------------------------------------
# Single-atom placement
# ---------------------------------------------------------------------------


def _place_one(
    mol, i, coords, parent, child_count, first_child_torsion, first_child_idx, bond_dihedral
):
    """Place atom i using z-matrix from its parent chain."""
    p = parent[i]
    if p is None:
        return

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
            mol, i, p, g, gg, child_count[p], first_child_torsion, in_ring, coords, bond_dihedral
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


# ---------------------------------------------------------------------------
# Ring closure
# ---------------------------------------------------------------------------


def _find_ring_systems(mol):
    """Return list of sets of atom indices forming fused ring systems."""
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        new = set(ring)
        merged = []
        for s in systems:
            if len(new & s) >= 1:
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
    placed = set()

    while remaining:
        ready = [
            a
            for a in remaining
            if parent[a] is None or parent[a] not in remaining or parent[a] in placed
        ]
        if not ready:
            ready = [min(remaining)]

        # Prefer atoms in same ring as the last placed atom
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


def _has_clash(mol, atoms, coords, placed, ring_set, threshold=0.8):
    """Check if any atom in `atoms` clashes with a placed atom outside the ring."""
    for i in atoms:
        for j in placed:
            if j in ring_set or mol.GetBondBetweenAtoms(i, j) is not None:
                continue
            if np.linalg.norm(coords[i] - coords[j]) < threshold:
                return True
    return False


def _close_ring(mol, ring, new_atoms, coords, parent, placed):
    """Adjust torsions of new_atoms to close a single ring.

    Only torsions are varied (bond angles stay at ideal values).
    Only atoms in new_atoms are moved (atoms from previously closed rings are frozen).

    Before optimizing, check for steric clashes with already-placed atoms.
    If a new atom collides, flip its torsion by 240° (chirality flip) and re-place
    downstream atoms in the ring.
    """
    from scipy.optimize import minimize

    ring_set = set(ring)

    # Find closure bonds: ring bonds that are not parent-child
    closures = []
    for idx in range(len(ring)):
        a, b = ring[idx], ring[(idx + 1) % len(ring)]
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond and parent[a] != b and parent[b] != a:
            closures.append((a, b, _get_bond_length(mol, a, b)))

    # Collect adjustable atoms: new_atoms that have a grandparent
    adjustable = []
    for i in sorted(new_atoms):
        p = parent[i]
        if p is None or parent[p] is None:
            continue
        g = parent[p]
        gg = parent[g]
        ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
        t = _measure_torsion(ref, coords[g], coords[p], coords[i])
        adjustable.append((i, t))

    if not adjustable:
        return

    # Before optimizing, check for steric clashes. If any new atom collides
    # with an already-placed atom, and the first adjustable atom's parent is SP3
    # (chirality ambiguity), flip chirality (shift torsion by -240°) and re-place.
    adj_atoms = [i for i, _ in adjustable]
    i0, _t0 = adjustable[0]
    do_flip = (
        _has_clash(mol, adj_atoms, coords, placed, ring_set)
        and mol.GetAtomWithIdx(parent[i0]).GetHybridization() == Chem.HybridizationType.SP3
    )
    if do_flip:
        p = parent[i0]
        g = parent[p]
        gg = parent[g]
        ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
        old_t = _measure_torsion(ref, coords[g], coords[p], coords[i0])
        new_t = old_t - 240.0
        coords[i0] = _place_atom(
            ref,
            coords[g],
            coords[p],
            _get_bond_length(mol, p, i0),
            _get_bond_angle(mol, g, p, i0),
            new_t,
        )
        adjustable[0] = (i0, new_t)
        # Re-place all downstream atoms
        for k2 in range(1, len(adjustable)):
            j = adjustable[k2][0]
            pj = parent[j]
            gj = parent[pj]
            ggj = parent[gj] if gj is not None else None
            ref_j = coords[ggj] if ggj is not None else _synthetic_ref(coords, gj, pj)
            t_j = _measure_torsion(ref_j, coords[gj], coords[pj], coords[j])
            coords[j] = _place_atom(
                ref_j,
                coords[gj],
                coords[pj],
                _get_bond_length(mol, pj, j),
                _get_bond_angle(mol, gj, pj, j),
                t_j,
            )
            adjustable[k2] = (j, t_j)

    # Check if closure optimization is needed
    if not closures:
        return
    max_gap = max(abs(np.linalg.norm(coords[a] - coords[b]) - bl) for a, b, bl in closures)
    if max_gap < 0.1:
        return

    init = np.array([t for _, t in adjustable])
    saved = coords[sorted(new_atoms)].copy()
    new_list = sorted(new_atoms)

    def _replace(torsions):
        for k, (i, _) in enumerate(adjustable):
            p = parent[i]
            g = parent[p]
            gg = parent[g]
            ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
            coords[i] = _place_atom(
                ref,
                coords[g],
                coords[p],
                _get_bond_length(mol, p, i),
                _get_bond_angle(mol, g, p, i),
                torsions[k],
            )

    # Find new angles formed by closure bonds (not set by z-matrix)
    angle_targets = []
    for a, b, _bl in closures:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            c = nb.GetIdx()
            if c != b and c in ring_set:
                angle_targets.append((c, a, b, _get_bond_angle(mol, c, a, b)))
        for nb in mol.GetAtomWithIdx(b).GetNeighbors():
            c = nb.GetIdx()
            if c != a and c in ring_set:
                angle_targets.append((a, b, c, _get_bond_angle(mol, a, b, c)))

    def _measure_angle(a, b, c):
        v1, v2 = coords[a] - coords[b], coords[c] - coords[b]
        cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-10)
        return np.degrees(np.arccos(np.clip(cos_a, -1, 1)))

    def cost(torsions):
        _replace(torsions)
        closure_cost = sum(
            (np.linalg.norm(coords[a] - coords[b]) - bl) ** 2 for a, b, bl in closures
        )
        reg_cost = 3e-4 * np.sum((torsions - init) ** 2)
        angle_cost = 3e-4 * sum(
            (_measure_angle(a, b, c) - ideal) ** 2 for a, b, c, ideal in angle_targets
        )
        return closure_cost + reg_cost + angle_cost

    init_cost = cost(init)
    result = minimize(cost, init, method="L-BFGS-B", options={"maxiter": 200, "ftol": 1e-10})

    if result.fun < init_cost:
        _replace(result.x)
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
    2. Optionally refine ring-system torsions to close rings.
    3. Place non-ring atoms.
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

    parent = [None] * n
    for i in range(1, n):
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() < i]
        if nbrs:
            parent[i] = max(nbrs)

    child_count = [0] * n
    first_child_torsion = [None] * n
    first_child_idx = [None] * n

    # Find ring systems, sorted by lowest atom index
    ring_systems = _find_ring_systems(mol)
    ring_systems.sort(key=lambda s: min(s))
    atom_to_system = {}
    for si, sys in enumerate(ring_systems):
        for a in sys:
            atom_to_system[a] = si

    # All individual rings (for per-ring closure)
    ri = mol.GetRingInfo()
    all_rings = [tuple(r) for r in ri.AtomRings()]

    placed = {0}
    placed_systems = set()
    closed_rings = set()
    closed_atoms: set[int] = set()

    for i in range(1, n):
        if parent[i] is None:
            continue

        if i in atom_to_system:
            si = atom_to_system[i]
            if si in placed_systems:
                continue
            # Place all atoms in this ring system
            order = _ring_visit_order(mol, ring_systems[si], parent)
            for j in order:
                if j in placed:
                    continue
                if parent[j] is None:
                    continue
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

                # Close any ring that just completed
                if refine_rings:
                    for ri_idx, ring in enumerate(all_rings):
                        if ri_idx in closed_rings:
                            continue
                        if set(ring).issubset(placed):
                            new_atoms = set(ring) - closed_atoms
                            _close_ring(mol, ring, new_atoms, coords, parent, placed)
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
