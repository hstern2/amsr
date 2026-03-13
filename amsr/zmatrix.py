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
                if a in pr and b in pr and pr is not ring:
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


def _choose_torsion(
    mol, p, i, g, gg, nth_child, first_child_torsion, in_ring, coords, bond_dihedral
):
    """Pick torsion for placing atom i as the nth child of p."""
    hyb_p = mol.GetAtomWithIdx(p).GetHybridization()

    if nth_child == 0:
        # check AMSR dihedral for g-p bond
        if (g, p) in bond_dihedral:
            mi, mj, angle = bond_dihedral[(g, p)]
            if mj == i and (gg is None or gg == mi):
                return angle
        gg_in_ring = True
        if gg is not None:
            b = mol.GetBondBetweenAtoms(gg, g)
            gg_in_ring = b is not None and b.IsInRing()
        return _default_torsion(mol, p, 0, in_ring, gg_in_ring)

    base = first_child_torsion[p]
    if hyb_p == Chem.HybridizationType.SP3:
        chiral = mol.GetAtomWithIdx(p).GetChiralTag()
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return base + 120.0 * nth_child
        if chiral == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return base - 120.0 * nth_child
        # Non-chiral: place H anti, heavy-atom children in gauche positions
        # Second child goes to the opposite gauche position from the first
        if base > 0:
            return base - 120.0 * nth_child
        return base + 120.0 * nth_child
    return base + 180.0


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
        torsion = _choose_torsion(
            mol, p, i, g, gg, child_count[p], first_child_torsion, in_ring, coords, bond_dihedral
        )
        if child_count[p] == 0:
            first_child_torsion[p] = torsion
        ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
        coords[i] = _place_atom(ref, coords[g], coords[p], bond_len, angle, torsion)

    if first_child_idx[p] is None:
        first_child_idx[p] = i

    # --- pick a ring neighbor j of i and place it to fix orientation -------
    ring_nbrs = [
        n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in sys_atoms
    ]

    # prefer the neighbor specified in an AMSR dihedral for the p-i axis
    j = ring_nbrs[0]
    orient_torsion = 0.0  # default: ring in plane of g-p-i
    if (p, i) in bond_dihedral:
        mi, mj, angle = bond_dihedral[(p, i)]
        if mj in sys_atoms:
            j = mj
            orient_torsion = angle

    bl_ij = _get_bond_length(mol, i, j)
    ang_pij = _get_bond_angle(mol, p, i, j)
    if g is not None:
        ref_j = coords[g]
    else:
        ref_j = _synthetic_ref(coords, p, i)
    coords_j = _place_atom(ref_j, coords[p], coords[i], bl_ij, ang_pij, orient_torsion)

    # --- build rotation: ideal frame → real frame -------------------------
    ideal_v = ideal[j] - ideal[i]
    ideal_d = ideal_v / np.linalg.norm(ideal_v)
    ideal_n = np.array([0.0, 0.0, 1.0])
    ideal_p = np.cross(ideal_n, ideal_d)
    ideal_p /= np.linalg.norm(ideal_p)

    real_v = coords_j - coords[i]
    real_d = real_v / np.linalg.norm(real_v)
    pi_vec = coords[i] - coords[p]
    real_n = np.cross(pi_vec, real_v)
    rn = np.linalg.norm(real_n)
    if rn < 1e-10:
        real_n = np.array([0.0, 0.0, 1.0])
    else:
        real_n /= rn
    real_p = np.cross(real_n, real_d)
    real_p /= np.linalg.norm(real_p)
    real_n = np.cross(real_d, real_p)

    M_ideal = np.column_stack([ideal_d, ideal_p, ideal_n])
    M_real = np.column_stack([real_d, real_p, real_n])
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
        child_count[pa] += 1
        if first_child_idx[pa] is None:
            first_child_idx[pa] = a
            ga = bfs_parent.get(pa)
            if ga is not None:
                gga = bfs_parent.get(ga)
                ref_a = coords[gga] if gga is not None else _synthetic_ref(coords, ga, pa)
                first_child_torsion[pa] = _measure_torsion(ref_a, coords[ga], coords[pa], coords[a])


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

    # find planar fused ring systems that need rigid placement
    all_systems = _find_ring_systems(mol)
    ri = mol.GetRingInfo()
    atom_to_rigid = {}  # atom -> frozenset
    rigid_ideal = {}  # frozenset -> {atom: coord}
    for sys in all_systems:
        n_rings = sum(1 for r in ri.AtomRings() if set(r) <= sys)
        if n_rings > 1 and _is_planar(mol, sys):
            ideal = _ideal_ring_coords(mol, sys)
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
            torsion = _choose_torsion(
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
            if gg is not None:
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

    # build RDKit conformer
    conf = Chem.Conformer(n_atoms)
    conf.Set3D(True)
    for i in range(n_atoms):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
