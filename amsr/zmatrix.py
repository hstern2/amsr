from collections import deque
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

_HYBRID_ANGLES = {
    Chem.HybridizationType.SP3: 109.5,
    Chem.HybridizationType.SP2: 120.0,
    Chem.HybridizationType.SP: 180.0,
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


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
    """Bond angle g-p-i.  Uses polygon angle for SP2 atoms in a ring."""
    hyb = mol.GetAtomWithIdx(p).GetHybridization()
    if hyb == Chem.HybridizationType.SP2:
        ri = mol.GetRingInfo()
        bond_gp = mol.GetBondBetweenAtoms(g, p)
        bond_pi = mol.GetBondBetweenAtoms(p, i)
        if bond_gp is not None and bond_pi is not None:
            common = set(ri.BondRingSizes(bond_gp.GetIdx())) & set(
                ri.BondRingSizes(bond_pi.GetIdx())
            )
            if common:
                n = min(common)
                return (n - 2) * 180.0 / n
    return _HYBRID_ANGLES.get(hyb, 109.5)


def _place_atom(A, B, C, d, theta_deg, omega_deg):
    """Place atom D given three reference points.

    d:         bond length C-D
    theta_deg: bond angle  B-C-D
    omega_deg: torsion     A-B-C-D
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


def _synthetic_ref(coords, g, p):
    """Synthetic great-grandparent when none exists."""
    gp = coords[p] - coords[g]
    gp_n = np.linalg.norm(gp)
    gp = gp / gp_n if gp_n > 1e-10 else np.array([1.0, 0.0, 0.0])
    perp = np.array([1.0, 0.0, 0.0]) if abs(gp[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    perp = np.cross(gp, perp)
    perp /= np.linalg.norm(perp)
    return coords[g] - perp


def _kabsch(P, Q, real_center, ideal_center):
    """Kabsch rotation: find R minimizing ||R @ (P - ideal_center) - (Q - real_center)||.

    P: (N, 3) ideal points, Q: (N, 3) real points.
    Returns rotation matrix R.
    """
    P_c = P - ideal_center
    Q_c = Q - real_center
    H = P_c.T @ Q_c
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1.0, 1.0, d])
    return Vt.T @ D @ U.T


def _measure_torsion(p0, p1, p2, p3):
    """Torsion angle (degrees) from four 3-D points."""
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1n = np.linalg.norm(n1)
    n2n = np.linalg.norm(n2)
    if n1n < 1e-10 or n2n < 1e-10:
        return 0.0
    n1 /= n1n
    n2 /= n2n
    b2u = b2 / np.linalg.norm(b2)
    return np.degrees(np.arctan2(np.dot(np.cross(n1, n2), b2u), np.dot(n1, n2)))


# ---------------------------------------------------------------------------
# Ring-system helpers
# ---------------------------------------------------------------------------


def _find_ring_systems(mol):
    """Return list of frozensets, each a set of atom indices forming a ring system.

    Rings sharing >= 2 atoms (an edge) are merged into one system.
    """
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    if not rings:
        return []
    systems = []
    for ring in rings:
        merged = []
        new = set(ring)
        for s in systems:
            if len(new & s) >= 2:
                new |= s
            else:
                merged.append(s)
        merged.append(new)
        systems = merged
    return [frozenset(s) for s in systems]


def _is_planar(mol, system):
    return all(
        mol.GetAtomWithIdx(a).GetHybridization() == Chem.HybridizationType.SP2 for a in system
    )


def _ideal_ring_coords_3d(mol, system_atoms):
    """3D ideal coordinates for a non-planar ring system using RDKit embedding.

    Returns dict {atom_idx: np.array([x, y, z])}.
    """
    from rdkit.Chem import AllChem

    atom_list = sorted(system_atoms)
    idx_map = {}  # original -> fragment
    emol = Chem.RWMol(Chem.Mol())
    for a in atom_list:
        orig_atom = mol.GetAtomWithIdx(a)
        new_atom = Chem.Atom(orig_atom.GetAtomicNum())
        new_atom.SetFormalCharge(orig_atom.GetFormalCharge())
        new_atom.SetNoImplicit(True)
        new_atom.SetNumExplicitHs(0)
        new_atom.SetChiralTag(orig_atom.GetChiralTag())
        idx_map[a] = emol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a1 in system_atoms and a2 in system_atoms:
            bt = bond.GetBondType()
            if bond.GetIsAromatic():
                bt = Chem.BondType.AROMATIC
            emol.AddBond(idx_map[a1], idx_map[a2], bt)

    # Set aromaticity flags on fragment atoms
    for a in atom_list:
        if mol.GetAtomWithIdx(a).GetIsAromatic():
            emol.GetAtomWithIdx(idx_map[a]).SetIsAromatic(True)

    # Clear chiral tags on atoms with too few neighbors in the fragment
    for a in atom_list:
        fi = idx_map[a]
        if emol.GetAtomWithIdx(fi).GetDegree() < 3:
            emol.GetAtomWithIdx(fi).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

    frag = emol.GetMol()
    try:
        Chem.SanitizeMol(frag)
    except Exception:
        # Kekulization may fail for fragments; sanitize without kekulization
        try:
            Chem.SanitizeMol(
                frag,
                Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
            )
        except Exception:
            return _ideal_ring_coords(mol, system_atoms)
    frag_h = Chem.AddHs(frag)
    res = AllChem.EmbedMolecule(frag_h, randomSeed=42)
    if res < 0:
        res = AllChem.EmbedMolecule(frag_h, randomSeed=42, useRandomCoords=True)
    if res < 0:
        return _ideal_ring_coords(mol, system_atoms)
    AllChem.MMFFOptimizeMolecule(frag_h)

    conf = frag_h.GetConformer()
    ideal = {}
    for a in atom_list:
        pos = conf.GetAtomPosition(idx_map[a])
        ideal[a] = np.array([pos.x, pos.y, pos.z])
    return ideal


def _ideal_ring_coords(mol, system_atoms):
    """Ideal planar coordinates for a ring system (single or fused).

    Returns dict {atom_idx: np.array([x, y, 0])}.
    """
    ri = mol.GetRingInfo()
    rings = [list(r) for r in ri.AtomRings() if set(r) <= set(system_atoms)]
    rings.sort(key=len, reverse=True)

    ideal = {}
    placed = set()

    def avg_bl(ring):
        n = len(ring)
        return sum(_get_bond_length(mol, ring[k], ring[(k + 1) % n]) for k in range(n)) / n

    # First ring — regular polygon
    r0 = rings[0]
    n0 = len(r0)
    bl = avg_bl(r0)
    radius = bl / (2 * np.sin(np.pi / n0))
    start = np.pi / 2 + np.pi / n0
    for k, idx in enumerate(r0):
        theta = start + 2 * np.pi * k / n0
        ideal[idx] = np.array([radius * np.cos(theta), radius * np.sin(theta), 0.0])
        placed.add(idx)

    # Remaining rings — extend from shared edge
    remaining = list(range(1, len(rings)))
    for _ in range(len(remaining) * 2):
        if not remaining:
            break
        for ri_idx in list(remaining):
            ring = rings[ri_idx]
            n = len(ring)
            # find two adjacent placed atoms
            edge_k = None
            for k in range(n):
                if ring[k] in placed and ring[(k + 1) % n] in placed:
                    edge_k = k
                    break
            if edge_k is None:
                continue

            ring_set = set(ring)
            ring = ring[edge_k:] + ring[:edge_k]
            a, b = ring[0], ring[1]
            bl = avg_bl(ring)
            r = bl / (2 * np.sin(np.pi / n))

            mid = (ideal[a] + ideal[b]) / 2
            ev = ideal[b] - ideal[a]
            el = np.linalg.norm(ev)
            ed = ev / el if el > 1e-10 else np.array([1, 0, 0])
            perp = np.array([-ed[1], ed[0], 0.0])

            # new ring goes on opposite side of shared edge from existing ring
            for pr in rings:
                if a in pr and b in pr and set(pr) != ring_set:
                    other = [x for x in pr if x in placed and x != a and x != b]
                    if other and np.dot(ideal[other[0]] - mid, perp) > 0:
                        perp = -perp
                    break

            center = mid + r * np.cos(np.pi / n) * perp
            ang_a = np.arctan2(ideal[a][1] - center[1], ideal[a][0] - center[0])
            ang_b = np.arctan2(ideal[b][1] - center[1], ideal[b][0] - center[0])
            step = 2 * np.pi / n
            diff = (ang_b - ang_a + np.pi) % (2 * np.pi) - np.pi
            if diff < 0:
                step = -step

            for k in range(n):
                atom = ring[k]
                if atom not in placed:
                    t = ang_a + step * k
                    ideal[atom] = center + r * np.array([np.cos(t), np.sin(t), 0.0])
                    placed.add(atom)

            remaining.remove(ri_idx)

    return ideal


# ---------------------------------------------------------------------------
# Torsion logic
# ---------------------------------------------------------------------------


def _default_torsion(mol, p, nth_child, in_ring, gg_in_ring):
    """Default torsion angle for placing the nth child of atom p."""
    hyb = mol.GetAtomWithIdx(p).GetHybridization()
    if hyb == Chem.HybridizationType.SP3:
        return 180.0 + 120.0 * nth_child
    if hyb == Chem.HybridizationType.SP2:
        if in_ring:
            return (
                (0.0 if nth_child == 0 else 180.0)
                if gg_in_ring
                else (180.0 if nth_child == 0 else 0.0)
            )
        return 180.0 if nth_child == 0 else 0.0
    return 180.0


def _pick_torsion_by_ring_closure(mol, p, i, g, gg, t_plus, t_minus, coords, bond_dihedral):
    """Try both candidate torsions for atom i, build the ring, pick best closure.

    Returns the better torsion or None if the check doesn't apply.
    """
    # Find the ring containing both p and i
    ri = mol.GetRingInfo()
    target_ring = None
    for ring in ri.AtomRings():
        if p in ring and i in ring:
            target_ring = list(ring)
            break
    if target_ring is None:
        return None

    # Find ring closure: a bond in the ring where neither end is parent of the other
    # Build parent chain within the ring starting from i
    ring_set = set(target_ring)
    # Trace the ring path from i back to p (through ring atoms)
    # The closure atom is the ring neighbor of p that isn't i
    closure_atom = None
    for nb in mol.GetAtomWithIdx(p).GetNeighbors():
        nidx = nb.GetIdx()
        if nidx in ring_set and nidx != i:
            closure_atom = nidx
            break
    if closure_atom is None:
        return None

    # Build the ring path from i to closure_atom (excluding p)
    # BFS through ring bonds, excluding p
    from collections import deque as _deque

    path_parent = {i: None}
    bfs_q = _deque([i])
    found = False
    while bfs_q and not found:
        u = bfs_q.popleft()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            v = nb.GetIdx()
            if v == p or v not in ring_set or v in path_parent:
                continue
            path_parent[v] = u
            if v == closure_atom:
                found = True
                break
            bfs_q.append(v)
    if not found:
        return None

    # Get ordered path from i to closure_atom
    path = []
    v = closure_atom
    while v is not None:
        path.append(v)
        v = path_parent[v]
    path.reverse()  # path[0] = i, path[-1] = closure_atom

    ref_pt = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
    angle_gpi = _get_bond_angle(mol, g, p, i)
    bl_pi = _get_bond_length(mol, p, i)

    best_t = None
    best_dist = float("inf")

    for t_cand in (t_plus, t_minus):
        # Place atom i at candidate torsion
        tmp_coords = (
            dict(coords)
            if isinstance(coords, dict)
            else {k: coords[k].copy() for k in range(len(coords))}
        )
        tmp_coords[i] = _place_atom(ref_pt, coords[g], coords[p], bl_pi, angle_gpi, t_cand)

        # Build out the ring path using AMSR dihedrals where available
        prev_prev = coords[g]
        prev = coords[p]
        curr = tmp_coords[i]
        for k in range(1, len(path)):
            a_prev = path[k - 1]
            a_curr = path[k]
            bl = _get_bond_length(mol, a_prev, a_curr)
            ang = _get_bond_angle(mol, path[k - 2] if k >= 2 else p, a_prev, a_curr)

            # Look for AMSR dihedral on bond (prev_atom, a_prev)
            gp_key = (path[k - 2] if k >= 2 else p, a_prev)
            torsion = None
            if gp_key in bond_dihedral:
                mi_d, mj_d, angle_d = bond_dihedral[gp_key]
                if mj_d == a_curr:
                    torsion = angle_d
            if torsion is None:
                bp = mol.GetBondBetweenAtoms(a_prev, a_curr)
                ir = bp is not None and bp.IsInRing()
                bp2 = (
                    mol.GetBondBetweenAtoms(gp_key[0], gp_key[1]) if gp_key[0] is not None else None
                )
                gg_ir = bp2 is not None and bp2.IsInRing() if bp2 else False
                torsion = _default_torsion(mol, a_prev, 0, ir, gg_ir)

            new_pos = _place_atom(prev_prev, prev, curr, bl, ang, torsion)
            tmp_coords[a_curr] = new_pos
            prev_prev = prev
            prev = curr
            curr = new_pos

        # Check tetrahedral angle deviation at p: closure_atom should form
        # proper tetrahedral angles with p's other neighbors.
        placed_nbrs = []
        for nb in mol.GetAtomWithIdx(p).GetNeighbors():
            nidx = nb.GetIdx()
            if nidx != i and nidx != closure_atom and (np.any(coords[nidx]) or nidx == 0):
                placed_nbrs.append(coords[nidx])
        # Add i and closure_atom from this candidate
        placed_nbrs.append(tmp_coords[i])
        placed_nbrs.append(tmp_coords[closure_atom])
        total_dev = 0.0
        for ia in range(len(placed_nbrs)):
            for ib in range(ia + 1, len(placed_nbrs)):
                va = placed_nbrs[ia] - coords[p]
                vb = placed_nbrs[ib] - coords[p]
                na, nb_ = np.linalg.norm(va), np.linalg.norm(vb)
                if na < 1e-10 or nb_ < 1e-10:
                    continue
                cos_a = np.dot(va, vb) / (na * nb_)
                ang = np.degrees(np.arccos(np.clip(cos_a, -1, 1)))
                total_dev += abs(ang - 109.5)
        if total_dev < best_dist:
            best_dist = total_dev
            best_t = t_cand

    return best_t


def _choose_torsion(
    mol, p, i, g, gg, nth_child, first_child_torsion, in_ring, coords, bond_dihedral
):
    """Pick torsion for placing atom i as the nth child of p.

    Returns (torsion_angle, ref_override) where ref_override is the atom index
    to use as the torsion reference point, or None to use the default (gg).
    """
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    if nth_child == 0:
        if (g, p) in bond_dihedral:
            mi, mj, angle = bond_dihedral[(g, p)]
            if mj == i:
                return angle, mi if mi != gg else None
            # The dihedral on g-p has mj != i.  mj is a different neighbor
            # of p whose torsion (mi-g-p-mj) is known.  Use it as the
            # reference to compute the torsion for atom i.
            if np.any(coords[mj]) or mj == 0:
                # Compute the AMSR torsion of mj relative to ref, then offset
                # to get the torsion for i
                actual_mj = _measure_torsion(
                    coords[gg] if gg is not None else _synthetic_ref(coords, g, p),
                    coords[g],
                    coords[p],
                    coords[mj],
                )
                hyb_p = mol.GetAtomWithIdx(p).GetHybridization()
                if hyb_p == Chem.HybridizationType.SP3:
                    return actual_mj + 120.0, None
                elif hyb_p == Chem.HybridizationType.SP2:
                    return actual_mj + 180.0, None
                else:
                    return actual_mj + 120.0, None

        gg_in_ring = True
        if gg is not None:
            b = mol.GetBondBetweenAtoms(gg, g)
            gg_in_ring = b is not None and b.IsInRing()
        return _default_torsion(mol, p, 0, in_ring, gg_in_ring), None

    base = first_child_torsion[p]
    if base is None:
        base = 0.0
    if hyb_p == Chem.HybridizationType.SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return base + 120.0 * nth_child, None
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return base - 120.0 * nth_child, None
        # Non-chiral: try both ±120° and pick the sign that gives better
        # ring closure when atom i is in a ring containing p.
        t_plus = base + 120.0 * nth_child
        t_minus = base - 120.0 * nth_child
        bp_pi = mol.GetBondBetweenAtoms(p, i)
        if bp_pi is not None and bp_pi.IsInRing() and g is not None and gg is not None:
            best_t = _pick_torsion_by_ring_closure(
                mol, p, i, g, gg, t_plus, t_minus, coords, bond_dihedral
            )
            if best_t is not None:
                return best_t, None
        # Fallback heuristic
        if base > 0:
            return base - 120.0 * nth_child, None
        return base + 120.0 * nth_child, None
    return base + 180.0, None


# ---------------------------------------------------------------------------
# Rigid ring-system placement
# ---------------------------------------------------------------------------


def _place_rigid_ring(
    mol,
    sys_atoms,
    ideal,
    coords,
    bfs_parent,
    child_count,
    first_child_torsion,
    first_child_idx,
    bond_dihedral,
    ext_parent,
    entry,
):
    """Orient and place a planar ring system as a rigid unit.

    ext_parent : already-placed atom bonded to entry
    entry      : first ring atom reached (not yet placed)
    """
    p = ext_parent
    i = entry
    g = bfs_parent[p]

    # --- place entry atom i via z-matrix from its external parent ----------
    bond_len = _get_bond_length(mol, p, i)
    if g is None:
        if child_count[p] == 0:
            coords[i] = coords[p] + np.array([bond_len, 0.0, 0.0])
        else:
            c1 = first_child_idx[p]
            hyb = mol.GetAtomWithIdx(p).GetHybridization()
            a = _HYBRID_ANGLES.get(hyb, 109.5)
            ref = _synthetic_ref(coords, p, c1)
            omega = 180.0 if hyb == Chem.HybridizationType.SP2 else 120.0 * child_count[p]
            coords[i] = _place_atom(ref, coords[c1], coords[p], bond_len, a, omega)
    else:
        gg = bfs_parent.get(g)
        angle = _get_bond_angle(mol, g, p, i)
        bp = mol.GetBondBetweenAtoms(p, i)
        in_ring = bp is not None and bp.IsInRing()
        torsion, ref_override = _choose_torsion(
            mol, p, i, g, gg, child_count[p], first_child_torsion, in_ring, coords, bond_dihedral
        )
        if child_count[p] == 0:
            first_child_torsion[p] = torsion
        if ref_override is not None:
            ref = coords[ref_override]
        elif gg is not None:
            ref = coords[gg]
        else:
            ref = _synthetic_ref(coords, g, p)
        coords[i] = _place_atom(ref, coords[g], coords[p], bond_len, angle, torsion)

    if first_child_idx[p] is None and p not in sys_atoms:
        first_child_idx[p] = i

    # --- pick a ring neighbor j of i and place it to fix orientation -------
    ring_nbrs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in sys_atoms
    ]

    # prefer the neighbor specified in an AMSR dihedral for the p-i axis
    j = ring_nbrs[0]
    orient_torsion = 0.0  # default: ring in plane of g-p-i
    orient_ref = None  # AMSR reference atom for orient_torsion
    if (p, i) in bond_dihedral:
        mi, mj, angle = bond_dihedral[(p, i)]
        if mj in sys_atoms:
            j = mj
            orient_torsion = angle
            orient_ref = mi  # the AMSR's reference atom

    bl_ij = _get_bond_length(mol, i, j)
    ang_pij = _get_bond_angle(mol, p, i, j)
    # Use the AMSR reference atom if available; otherwise fall back to g or synthetic
    if orient_ref is not None and orient_ref in coords and np.any(coords[orient_ref]):
        ref_j = coords[orient_ref]
    elif g is not None:
        ref_j = coords[g]
    else:
        ref_j = _synthetic_ref(coords, p, i)
    coords_j = _place_atom(ref_j, coords[p], coords[i], bl_ij, ang_pij, orient_torsion)

    # --- build rotation: ideal frame → real frame -------------------------
    # Use Kabsch alignment when we have 3+ known points (e.g., ext_parent
    # is in the ring system — spiro case), otherwise use 2-vector frame.
    if p in sys_atoms and p in ideal:
        # 3-point Kabsch: entry (i), orient neighbor (j), and anchored parent (p)
        ideal_pts = np.array([ideal[i], ideal[j], ideal[p]])
        real_pts = np.array([coords[i], coords_j, coords[p]])
        R = _kabsch(ideal_pts, real_pts, coords[i], ideal[i])
    else:
        ideal_v = ideal[j] - ideal[i]
        ideal_d = ideal_v / np.linalg.norm(ideal_v)
        # Check if ideal coords are planar (z ≈ 0)
        z_spread = max(abs(ideal[a][2]) for a in sys_atoms)
        if z_spread < 0.01:
            ideal_n = np.array([0.0, 0.0, 1.0])
        else:
            pts = np.array([ideal[a] for a in sys_atoms])
            centroid = pts.mean(axis=0)
            _, _, Vt = np.linalg.svd(pts - centroid)
            ideal_n = Vt[-1]
            ideal_n = ideal_n - np.dot(ideal_n, ideal_d) * ideal_d
            nn = np.linalg.norm(ideal_n)
            ideal_n = ideal_n / nn if nn > 1e-10 else np.array([0.0, 0.0, 1.0])
        ideal_p_vec = np.cross(ideal_n, ideal_d)
        ideal_p_vec /= np.linalg.norm(ideal_p_vec)

        real_v = coords_j - coords[i]
        real_d = real_v / np.linalg.norm(real_v)
        pi_vec = coords[i] - coords[p]
        real_n = np.cross(pi_vec, real_v)
        rn = np.linalg.norm(real_n)
        if rn < 1e-10:
            real_n = np.array([0.0, 0.0, 1.0])
        else:
            real_n /= rn
        real_p_vec = np.cross(real_n, real_d)
        real_p_vec /= np.linalg.norm(real_p_vec)
        real_n = np.cross(real_d, real_p_vec)

        M_ideal = np.column_stack([ideal_d, ideal_p_vec, ideal_n])
        M_real = np.column_stack([real_d, real_p_vec, real_n])
        R = M_real @ np.linalg.inv(M_ideal)

    for a in sys_atoms:
        coords[a] = R @ (ideal[a] - ideal[i]) + coords[i]

    # --- set BFS parents for ring atoms (mini-BFS from entry) -------------
    bfs_parent[i] = p
    ring_visited = {i}
    rq = deque([i])
    while rq:
        u = rq.popleft()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            v = nb.GetIdx()
            if v in sys_atoms and v not in ring_visited:
                bfs_parent[v] = u
                ring_visited.add(v)
                rq.append(v)

    # --- back-fill child_count / first_child_torsion for ring atoms -------
    #     so that later substituent placement off ring atoms works
    for a in sys_atoms:
        pa = bfs_parent[a]
        if pa not in sys_atoms:
            continue
        # Skip the entry atom — its relationship with ext_parent is handled
        # by the main loop, not the ring BFS.
        if a == i and pa == p:
            continue
        child_count[pa] += 1
        if first_child_idx[pa] is None:
            first_child_idx[pa] = a
            ga = bfs_parent.get(pa)
            if ga is not None:
                gga = bfs_parent.get(ga)
                # Avoid degenerate torsion from cyclic parent pointers
                if gga is not None and gga == pa:
                    first_child_torsion[pa] = 0.0
                else:
                    ref_a = coords[gga] if gga is not None else _synthetic_ref(coords, ga, pa)
                    first_child_torsion[pa] = _measure_torsion(
                        ref_a, coords[ga], coords[pa], coords[a]
                    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
) -> Chem.Mol:
    """Generate 3D conformer with ideal geometry and AMSR dihedrals.

    Algorithm
    ---------
    1. Identify planar fused ring systems; compute ideal polygon coords.
    2. BFS from atom 0, placing each atom with ideal bond length / angle.
       Ring systems are oriented and placed as rigid units.
    3. Torsion angles come from the *dihedral* dict (AMSR) when available,
       otherwise from hybridisation defaults.
    """
    n_atoms = mol.GetNumAtoms()
    if n_atoms == 0:
        return mol

    coords = np.zeros((n_atoms, 3))

    # dihedral lookup: (i, j) -> (mi, mj, angle)
    bond_dihedral = {}
    if dihedral:
        for (mi, i, j, mj), angle in dihedral.items():
            bond_dihedral[(i, j)] = (mi, mj, angle)
            bond_dihedral[(j, i)] = (mj, mi, angle)

    # find ring systems that need rigid placement
    all_systems = _find_ring_systems(mol)
    ri = mol.GetRingInfo()
    atom_to_rigid = {}  # atom -> frozenset
    rigid_ideal = {}  # frozenset -> {atom: coord}
    for sys in all_systems:
        n_rings = sum(1 for r in ri.AtomRings() if set(r) <= sys)
        planar = _is_planar(mol, sys)
        # Rigid placement for multi-ring systems.  For non-planar single
        # rings without AMSR dihedrals, also use rigid 3D placement since
        # the atom-by-atom ring closure correction distorts bond lengths.
        # Planar single rings and non-planar single rings WITH AMSR
        # dihedrals (e.g. cyclohexane in spiro systems) are handled
        # atom-by-atom.
        if n_rings > 1:
            do_rigid = True
        elif n_rings == 1 and not planar:
            has_dih = any(
                (a, b) in bond_dihedral
                for a in sys
                for b in sys
                if mol.GetBondBetweenAtoms(a, b) is not None
            )
            do_rigid = not has_dih
        else:
            do_rigid = False
        if do_rigid:
            # Count non-SP2 atoms; rings with at most 1 are ~planar.
            n_non_sp2 = sum(
                1
                for a in sys
                if mol.GetAtomWithIdx(a).GetHybridization() != Chem.HybridizationType.SP2
            )
            if planar or n_non_sp2 <= 1:
                ideal = _ideal_ring_coords(mol, sys)
            else:
                ideal = _ideal_ring_coords_3d(mol, sys)
            rigid_ideal[sys] = ideal
            for a in sys:
                atom_to_rigid[a] = sys

    # parent[i] = max(neighbor index < i), matching AMSR's DFS output order
    parent = [None] * n_atoms
    for i in range(1, n_atoms):
        nbrs_before = [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() < i]
        if nbrs_before:
            parent[i] = max(nbrs_before)

    placed_rigid: set[int] = set()
    placed = set()
    child_count = [0] * n_atoms
    first_child_torsion = [None] * n_atoms
    first_child_idx = [None] * n_atoms

    for i in range(n_atoms):
        if i in placed:
            continue
        p = parent[i]
        if p is None:
            placed.add(i)
            continue

        # --- rigid ring system? ----------------------------------------
        sys_key = atom_to_rigid.get(i)
        if sys_key is not None and sys_key not in placed_rigid:
            # Use parent dict as bfs_parent for the rigid placer
            parent_dict = dict(enumerate(parent))
            _place_rigid_ring(
                mol,
                sys_key,
                rigid_ideal[sys_key],
                coords,
                parent_dict,
                child_count,
                first_child_torsion,
                first_child_idx,
                bond_dihedral,
                p,
                i,
            )
            # update parent list from rigid ring's mini-BFS
            for a in sys_key:
                parent[a] = parent_dict[a]
            placed_rigid.add(sys_key)
            for a in sys_key:
                placed.add(a)
            child_count[p] += 1
            if first_child_idx[p] is None:
                first_child_idx[p] = i
                # Measure first child torsion if possible
                g_p = parent[p]
                if g_p is not None:
                    gg_p = parent[g_p]
                    ref_p = coords[gg_p] if gg_p is not None else _synthetic_ref(coords, g_p, p)
                    first_child_torsion[p] = _measure_torsion(
                        ref_p, coords[g_p], coords[p], coords[i]
                    )
                else:
                    first_child_torsion[p] = 0.0
            continue

        # --- single atom z-matrix placement ----------------------------
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
            angle = _get_bond_angle(mol, g, p, i)
            torsion, ref_override = _choose_torsion(
                mol,
                p,
                i,
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
                if mj == i:
                    ref = coords[mi]
                else:
                    ref = _synthetic_ref(coords, g, p)
            else:
                ref = _synthetic_ref(coords, g, p)
            coords[i] = _place_atom(ref, coords[g], coords[p], bond_len, angle, torsion)
            if child_count[p] == 0:
                # Measure first-child torsion in the standard reference frame
                # so subsequent children use a consistent reference
                std_ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                first_child_torsion[p] = _measure_torsion(std_ref, coords[g], coords[p], coords[i])

        if first_child_idx[p] is None:
            first_child_idx[p] = i
        placed.add(i)
        child_count[p] += 1

    # --- ring closure: re-place ring atoms with optimized bond angles ------
    for sys in all_systems:
        if sys in placed_rigid:
            continue
        # Find ALL closure bonds in this ring system
        closures = []
        for bond in mol.GetBonds():
            a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if a1 in sys and a2 in sys and bond.IsInRing():
                if parent[a1] != a2 and parent[a2] != a1:
                    closures.append((a1, a2, _get_bond_length(mol, a1, a2)))
        if not closures:
            continue
        # Check if any closure needs fixing
        max_gap = max(abs(np.linalg.norm(coords[a] - coords[b]) - bl) for a, b, bl in closures)
        if max_gap < 0.05:
            continue

        # Build tree path through ring system (BFS from entry)
        entry = None
        for a in sorted(sys):
            if parent[a] is None or parent[a] not in sys:
                entry = a
                break
        if entry is None:
            continue
        # BFS to get all ring atoms in tree order
        ring_path = [entry]
        visited = {entry}
        bfs_q = deque([entry])
        while bfs_q:
            u = bfs_q.popleft()
            for nb in mol.GetAtomWithIdx(u).GetNeighbors():
                v = nb.GetIdx()
                if v in sys and v not in visited and parent[v] == u:
                    ring_path.append(v)
                    visited.add(v)
                    bfs_q.append(v)

        atoms_to_replace = ring_path[1:]  # entry stays fixed
        if not atoms_to_replace:
            continue
        n_repl = len(atoms_to_replace)
        saved_coords = {a: coords[a].copy() for a in atoms_to_replace}

        def _replace_with_offsets(offsets):
            """Re-place ring atoms using original torsions but adjusted angles."""
            for a in atoms_to_replace:
                coords[a] = saved_coords[a]  # reset
            for k, a in enumerate(atoms_to_replace):
                p_a = parent[a]
                g_a = parent[p_a] if p_a is not None else None
                if g_a is None:
                    continue
                gg_a = parent[g_a] if g_a is not None else None
                bl = _get_bond_length(mol, p_a, a)
                angle = _get_bond_angle(mol, g_a, p_a, a) + offsets[k]
                # Re-use the same torsion as original placement
                bp_a = mol.GetBondBetweenAtoms(p_a, a)
                in_ring_a = bp_a is not None and bp_a.IsInRing()
                torsion, ref_override = _choose_torsion(
                    mol,
                    p_a,
                    a,
                    g_a,
                    gg_a,
                    0,  # treat as first child for torsion lookup
                    first_child_torsion,
                    in_ring_a,
                    coords,
                    bond_dihedral,
                )
                if ref_override is not None:
                    ref = coords[ref_override]
                elif gg_a is not None:
                    ref = coords[gg_a]
                else:
                    ref = _synthetic_ref(coords, g_a, p_a)
                coords[a] = _place_atom(ref, coords[g_a], coords[p_a], bl, angle, torsion)

        # Gauss-Newton: minimize sum of closure distance errors
        offsets = np.zeros(n_repl)
        eps = 0.1
        n_closures = len(closures)

        for _ in range(10):
            _replace_with_offsets(offsets)
            residual = np.array(
                [np.linalg.norm(coords[a] - coords[b]) - bl for a, b, bl in closures]
            )
            if np.max(np.abs(residual)) < 0.05:
                break
            # Jacobian
            J = np.zeros((n_closures, n_repl))
            for k in range(n_repl):
                offsets_k = offsets.copy()
                offsets_k[k] += eps
                _replace_with_offsets(offsets_k)
                for ci, (a, b, bl) in enumerate(closures):
                    J[ci, k] = (np.linalg.norm(coords[a] - coords[b]) - bl - residual[ci]) / eps
            _replace_with_offsets(offsets)  # restore
            delta, _, _, _ = np.linalg.lstsq(J, -residual, rcond=None)
            offsets = offsets + delta

        _replace_with_offsets(offsets)

    # build RDKit conformer
    conf = Chem.Conformer(n_atoms)
    conf.Set3D(True)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)

    return mol.GetMol()
