"""Conformer generation from AMSR dihedrals.

Ring atoms are placed by optimizing a Cartesian-space cost function that
enforces ideal bond lengths, bond angles, planarity at SP2 centers,
chirality at SP3 centers, and AMSR dihedral restraints.

Non-ring atoms are placed sequentially via z-matrix from their parent chain.

Geometry primitives and cost-function components use only numpy arrays
(no RDKit), making them suitable for reimplementation in C/C++.
"""

import math
from collections import deque
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

SP = Chem.HybridizationType.SP
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


def _set_dihedral(coords, mi, i, j, mj, target, atoms_to_rotate):
    """Rotate atoms_to_rotate around bond i-j to set torsion(mi,i,j,mj) = target."""
    current = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
    delta = np.radians(target - current)
    axis = coords[j] - coords[i]
    axis = axis / (_norm3(axis) + 1e-10)
    cos_d, sin_d = np.cos(delta), np.sin(delta)
    origin = coords[i]
    for a in atoms_to_rotate:
        v = coords[a] - origin
        coords[a] = (
            origin + v * cos_d + _cross3(axis, v) * sin_d + axis * np.dot(axis, v) * (1 - cos_d)
        )


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
    # Single C-N bonds between two SP2 atoms have amide resonance character.
    if bo == 1 and {s1, s2} == {"C", "N"}:
        ai, aj = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        if ai.GetHybridization() == SP2 and aj.GetHybridization() == SP2:
            return 0.5 * (_BOND_LENGTHS[("C", "N", 1)] + _BOND_LENGTHS[("C", "N", 1.5)])
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
    ri = mol.GetRingInfo()
    b_ab = mol.GetBondBetweenAtoms(a, b)
    b_bc = mol.GetBondBetweenAtoms(b, c)
    if b_ab is not None and b_bc is not None:
        common = set(ri.BondRingSizes(b_ab.GetIdx())) & set(ri.BondRingSizes(b_bc.GetIdx()))
        if common:
            n = min(common)
            poly = (n - 2) * 180.0 / n
            if hyb in (SP2, Chem.HybridizationType.SP) and n <= 6:
                return poly
            # SP3 in mixed rings (containing SP2 atoms): use polygon angle
            # so the ring angle sum is consistent with the SP2 atoms.
            if hyb == SP3 and n <= 6 and poly < 109.5:
                for ring in ri.AtomRings():
                    if len(ring) == n and b in ring:
                        if any(
                            mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in ring if x != b
                        ):
                            return poly
                        break
    return _HYBRID_ANGLES.get(hyb, 109.5)


def _bfs_subtree(mol, root, exclude):
    """BFS from root, excluding the given atom.  Returns set of reachable atoms."""
    visited = {root}
    queue = deque([root])
    while queue:
        curr = queue.popleft()
        for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
            b = nb.GetIdx()
            if b != exclude and b not in visited:
                visited.add(b)
                queue.append(b)
    return visited


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


def _find_core_atoms(mol, ring_atoms, n):
    """Minimal connected subgraph spanning all ring atoms.

    Core = ring atoms + bridging chain atoms between ring systems.
    Identified by starting with all atoms and iteratively pruning
    non-ring leaves.
    """
    core = set(range(n))
    changed = True
    while changed:
        changed = False
        for a in list(core):
            if a in ring_atoms:
                continue
            nbrs_in_core = sum(
                1 for nb in mol.GetAtomWithIdx(a).GetNeighbors() if nb.GetIdx() in core
            )
            if nbrs_in_core <= 1:
                core.discard(a)
                changed = True
    return core


def _fix_equivalent_terminals(mol, bond_dihedral, coords):
    """Fix dihedral reference-atom mismatches for equivalent terminals.

    When the reference atom is a degree-1 terminal with an equivalent sibling
    (same element, also degree-1), the encode/decode may pick different ones.
    Try each alternative and keep the better match.  Mutates bond_dihedral.
    """
    for (i, j), (mi, mj, angle) in list(bond_dihedral.items()):
        if i > j:
            continue
        for side, ref, bond_end, other_end in [(0, mi, i, j), (1, mj, j, i)]:
            a_ref = mol.GetAtomWithIdx(ref)
            if a_ref.GetDegree() != 1:
                continue
            siblings = [
                nb.GetIdx()
                for nb in mol.GetAtomWithIdx(bond_end).GetNeighbors()
                if nb.GetIdx() != other_end
                and nb.GetIdx() != ref
                and nb.GetDegree() == 1
                and nb.GetAtomicNum() == a_ref.GetAtomicNum()
            ]
            if not siblings:
                continue
            cur_mi, cur_mj = bond_dihedral[(i, j)][:2]
            best_ref = ref
            cur_diff = abs(
                (
                    measure_torsion(coords[cur_mi], coords[i], coords[j], coords[cur_mj])
                    - angle
                    + 180
                )
                % 360
                - 180
            )
            for alt in siblings:
                alt_mi = alt if side == 0 else cur_mi
                alt_mj = alt if side == 1 else cur_mj
                alt_diff = abs(
                    (
                        measure_torsion(coords[alt_mi], coords[i], coords[j], coords[alt_mj])
                        - angle
                        + 180
                    )
                    % 360
                    - 180
                )
                if alt_diff < cur_diff:
                    best_ref = alt
                    cur_diff = alt_diff
            if best_ref != ref:
                mi_new = best_ref if side == 0 else mi
                mj_new = best_ref if side == 1 else mj
                bond_dihedral[(i, j)] = (mi_new, mj_new, angle)
                bond_dihedral[(j, i)] = (mj_new, mi_new, angle)


def _correct_junction_dihedrals(mol, ring_systems, core_atoms, bond_dihedral, coords):
    """Rotate atoms in the embedding to match AMSR dihedrals at non-ring bonds.

    RDKit distance geometry often places separate ring systems at arbitrary
    torsion angles.  For each AMSR dihedral on a core bond that is not in
    a ring, BFS to find the downstream atoms and rotate them to match the
    encoded angle.  Mutates coords in place.
    """
    if len(ring_systems) <= 1:
        return

    for (i, j), (mi, mj, angle) in bond_dihedral.items():
        if i > j or i not in core_atoms or j not in core_atoms:
            continue
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is None or bond.IsInRing():
            continue
        # Skip bonds where an SP atom makes the torsion undefined.
        if mol.GetAtomWithIdx(i).GetHybridization() == SP:
            continue
        if mol.GetAtomWithIdx(j).GetHybridization() == SP:
            continue
        to_rotate = _bfs_subtree(mol, j, i)
        if mi in to_rotate or mj not in to_rotate:
            continue
        _set_dihedral(coords, mi, i, j, mj, angle, to_rotate)


def _build_outward_tree(mol, placed, n, ring_atoms):
    """BFS from placed (core) atoms outward to build parent tree and visit order.

    Returns (outward_parent, branch_order).  Mutates placed to include
    disconnected component seeds.
    """
    outward_parent: list[Optional[int]] = [None] * n
    branch_order: list[int] = []

    # BFS within the core to assign parents
    core_root = min(ring_atoms) if ring_atoms else (min(placed) if placed else 0)
    visited = {core_root}
    bfs = [core_root]
    qi = 0
    while qi < len(bfs):
        curr = bfs[qi]
        qi += 1
        for nb in mol.GetAtomWithIdx(curr).GetNeighbors():
            b = nb.GetIdx()
            if b in placed and b not in visited:
                outward_parent[b] = curr
                visited.add(b)
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

    return outward_parent, branch_order


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


def _to_array(lst, dtype=int, cols=None):
    """Convert list to numpy array, returning an appropriately shaped empty array if empty."""
    if lst:
        return np.array(lst, dtype=dtype)
    return np.empty((0, cols) if cols else (0,), dtype=dtype)


def _collect_ring_bonds(mol, sys_set, fixed):
    """Collect all bonds within the ring system and to fixed neighbors.

    Returns (pairs, ideal_lengths) — both numpy arrays.
    """
    pairs, ideals = [], []
    for a in sorted(sys_set):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            # b > a avoids double-counting intra-system bonds; fixed atoms
            # are external so the guard is unnecessary for them.
            if (b in sys_set and b > a) or b in fixed:
                pairs.append((a, b))
                ideals.append(_get_bond_length(mol, a, b))
    return np.array(pairs, dtype=int), np.array(ideals)


def _collect_ring_angles(mol, sys_set, fixed):
    """Collect all bond angle triples involving ring system atoms.

    Includes angles at ring atoms AND angles at fixed atoms that have
    two ring-system neighbors.  For SP2 atoms with exactly 3 angles,
    adjusts targets so they sum to 360° (important for fused ring junctions).
    Skips SP atoms (180° targets) — those are handled by _collect_linear_atoms.
    Returns (triples, ideal_angles) — numpy arrays.
    """

    available = sys_set | set(fixed)
    triples, ideals = [], []
    center_indices: dict[int, list[int]] = {}  # center atom → list of indices into triples
    # Angles at ring atoms (skip SP — handled by linearity constraint)
    for b in sorted(sys_set):
        if mol.GetAtomWithIdx(b).GetHybridization() == SP:
            continue
        nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(b).GetNeighbors() if nb.GetIdx() in available
        ]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                idx = len(triples)
                triples.append((a, b, c))
                ideals.append(_get_bond_angle(mol, a, b, c))
                center_indices.setdefault(b, []).append(idx)
    # Angles at fixed atoms with >=2 ring neighbors
    for fb in sorted(fixed):
        ring_nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(fb).GetNeighbors() if nb.GetIdx() in sys_set
        ]
        if len(ring_nbrs) >= 2:
            for ia in range(len(ring_nbrs)):
                for ic in range(ia + 1, len(ring_nbrs)):
                    a, c = ring_nbrs[ia], ring_nbrs[ic]
                    idx = len(triples)
                    triples.append((a, fb, c))
                    ideals.append(_get_bond_angle(mol, a, fb, c))
                    center_indices.setdefault(fb, []).append(idx)
    # For SP2 atoms with exactly 3 angle targets, ensure they sum to 360°.
    # At fused ring junctions the cross-ring angle has no common ring size
    # and defaults to 120°, but the correct value is 360° minus the two
    # in-ring angles.  Only apply when exactly one angle is cross-ring.
    ri = mol.GetRingInfo()
    for b, indices in center_indices.items():
        if len(indices) != 3:
            continue
        if mol.GetAtomWithIdx(b).GetHybridization() != SP2:
            continue
        total = sum(ideals[k] for k in indices)
        if abs(total - 360.0) < 1.0:
            continue
        # Find angles without a common ring size between their two bonds.
        # Distribute the deficit among cross-ring angles to make the sum 360°.
        # Skip if any in-ring angle uses a ring > 6 (polygon formula unreliable).
        no_common = []
        has_large_ring = False
        for k in indices:
            a, _, c = triples[k]
            b_ab = mol.GetBondBetweenAtoms(a, b)
            b_bc = mol.GetBondBetweenAtoms(b, c)
            if b_ab is None or b_bc is None:
                continue
            common = set(ri.BondRingSizes(b_ab.GetIdx())) & set(ri.BondRingSizes(b_bc.GetIdx()))
            if not common:
                no_common.append(k)
            elif min(common) > 6:
                has_large_ring = True
        if no_common and not has_large_ring:
            per = (360.0 - total) / len(no_common)
            for k in no_common:
                ideals[k] += per
    return _to_array(triples, cols=3), np.array(ideals) if ideals else np.empty(0)


def _collect_planar_atoms(mol, sys_set, fixed):
    """Collect planarity constraints for SP2 atoms in the ring system.

    For each SP2 atom with 3+ neighbors whose coords are available,
    returns (center, a, b, c) tuples — all four should be coplanar.
    """
    available = sys_set | set(fixed)
    groups = []
    for j in sorted(sys_set):
        atom = mol.GetAtomWithIdx(j)
        if atom.GetHybridization() != SP2:
            continue
        # Skip atoms with chirality tags — they are pyramidal, not planar.
        if atom.GetChiralTag() in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in available]
        if len(nbrs) >= 3:
            groups.append((j, nbrs[0], nbrs[1], nbrs[2]))
    return _to_array(groups, cols=4)


def _collect_linear_atoms(mol, sys_set, fixed):
    """Collect linearity constraints for SP atoms in the ring system.

    Returns Nx3 int array of (a, b, c) triples where b is SP and
    the angle a-b-c should be 180°.
    """

    available = sys_set | set(fixed)
    triples = []
    for j in sorted(sys_set):
        if mol.GetAtomWithIdx(j).GetHybridization() != SP:
            continue
        nbrs = [
            nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() in available
        ]
        if len(nbrs) == 2:
            triples.append((nbrs[0], j, nbrs[1]))
    return _to_array(triples, cols=3)


def _collect_chiral_atoms(mol, sys_set, fixed, coords=None):
    """Collect chirality constraints for SP3 chiral atoms in the ring system.

    Returns (Nx5 int array (center, a, b, c, sign), N float array of target volumes).
    When coords are available, target volumes are the embedding volumes;
    otherwise ±1 from the chirality tag.
    """
    available = sys_set | set(fixed)
    result = []
    target_vols = []

    def _get(a):
        return coords[a] if a not in fixed else fixed[a]

    for j in sorted(sys_set):
        atom = mol.GetAtomWithIdx(j)
        chiral = atom.GetChiralTag()
        if chiral not in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in available]
        if len(nbrs) < 3:
            continue
        if coords is not None:
            rj = _get(j)
            # For neighbors that are SP, the embedding position reflects
            # a bent alkyne.  Replace with the linearized direction at the
            # ideal bond length so the target volume matches the geometry
            # the optimizer will converge to.
            positions = []
            for nb in nbrs[:3]:
                if mol.GetAtomWithIdx(nb).GetHybridization() == SP:
                    # Find the other neighbor of the SP atom (the one that isn't j)
                    sp_other = [
                        x.GetIdx() for x in mol.GetAtomWithIdx(nb).GetNeighbors() if x.GetIdx() != j
                    ]
                    if sp_other:
                        # Direction: from center toward the SP chain
                        far = _get(sp_other[0])
                        direction = far - rj
                        n = _norm3(direction)
                        if n > 1e-10:
                            direction /= n
                            d = _get_bond_length(mol, j, nb)
                            positions.append(rj + direction * d)
                        else:
                            positions.append(_get(nb))
                    else:
                        positions.append(_get(nb))
                else:
                    positions.append(_get(nb))
            v1, v2, v3 = positions[0] - rj, positions[1] - rj, positions[2] - rj
            vol = np.dot(v1, np.cross(v2, v3))
            # When the embedding volume is near zero (e.g. pyramidal atom
            # embedded as planar) or an SP neighbor makes the linearized
            # volume unrealistically large, use sign-only target.
            denom = _norm3(v1) * _norm3(v2) * _norm3(v3)
            oop = abs(vol) / denom if denom > 1e-10 else 0.0
            has_sp_nbr = any(mol.GetAtomWithIdx(nb).GetHybridization() == SP for nb in nbrs[:3])
            if oop > 0.1 and not has_sp_nbr:
                sign = 1 if vol > 0 else -1
                target_vols.append(vol)
            else:
                sign = -1 if chiral == CW else 1
                target_vols.append(float(sign))
        else:
            sign = -1 if chiral == CW else 1
            target_vols.append(float(sign))
        result.append((j, nbrs[0], nbrs[1], nbrs[2], sign))
    idx = _to_array(result, cols=5)
    tvol = _to_array(target_vols, dtype=float)
    return idx, tvol


def _collect_ring_dihedrals(mol, sys_set, bond_dihedral, fixed):
    """Collect AMSR dihedral restraints for ring bonds and boundary bonds.

    Includes dihedrals for bonds within the ring system and bonds
    connecting ring atoms to fixed (placed) atoms.
    Skips dihedrals where an SP atom makes the torsion angle undefined.
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
                    # Skip if either bond endpoint is SP — the 180° angle
                    # makes the torsion undefined.
                    i, j = bond_key
                    if mol.GetAtomWithIdx(i).GetHybridization() == SP:
                        continue
                    if mol.GetAtomWithIdx(j).GetHybridization() == SP:
                        continue
                    quads.append((mi, i, j, mj))
                    targets.append(float(angle))
    return _to_array(quads, cols=4), _to_array(targets, dtype=float)


def _collect_ring_planarity_dihedrals(mol, sys_set, fixed):
    """Collect 0° torsion constraints for planar ring bonds.

    For each ring, constrains consecutive 4-atom sequences where all
    four atoms are SP2 to 0° torsion (planar).  This handles both
    all-SP2 rings and mixed SP2/SP3 rings — the SP2 portions stay flat
    while SP3 atoms are free to puck.  For large (>6) all-SP2 rings,
    only constrains Kekulized double bonds since single bonds can rotate.
    Returns (quads, targets) — Nx4 int, N float.
    """
    available = sys_set | set(fixed)
    # Kekulize a copy to identify single vs double bonds in aromatic rings.
    mol_k = Chem.RWMol(mol)
    try:
        Chem.Kekulize(mol_k, clearAromaticFlags=False)
    except Exception:
        mol_k = mol
    quads, targets = [], []
    seen = set()
    for ring in mol.GetRingInfo().AtomRings():
        if not all(a in available for a in ring):
            continue
        n = len(ring)
        all_sp2 = all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring)
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            # All four atoms must be SP2 for planarity.
            if not all(mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in (a, b, c, d)):
                continue
            # For large all-SP2 rings, skip Kekulized single bonds (they can rotate).
            if all_sp2 and n > 6:
                bond_bc = mol_k.GetBondBetweenAtoms(b, c)
                if bond_bc is not None and bond_bc.GetBondType() == Chem.BondType.SINGLE:
                    continue
            key = (min(b, c), max(b, c))
            if key not in seen:
                seen.add(key)
                quads.append((a, b, c, d))
                targets.append(0.0)
    return _to_array(quads, cols=4), _to_array(targets, dtype=float)


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
    return _to_array(quads, cols=4), _to_array(targets, dtype=float)


# ============================================================
# Cartesian ring geometry: RDKit embedding for initialization
# ============================================================


def _rdkit_embed(mol, n_confs=1, seed=42):
    """Embed molecule with RDKit distance geometry (adds/removes Hs internally).

    Returns list of Nx3 coordinate arrays (one per conformer), or empty list
    on failure.  Heavy-atom indices match the input mol.
    """
    from rdkit.Chem import AllChem

    # Strip E/Z stereo before embedding — pseudo-E/Z on bonds with
    # equivalent substituents causes pathological slowdowns in RDKit's
    # distance geometry.  The optimizer enforces E/Z via its own constraints.
    mol_e = Chem.RWMol(mol)
    for b in mol_e.GetBonds():
        if b.GetStereo() != Chem.BondStereo.STEREONONE:
            b.SetStereo(Chem.BondStereo.STEREONONE)
    mol_h = Chem.AddHs(mol_e)
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
_W_CHIRAL = 10.0
_W_DIHEDRAL = 0.1
_W_EZ = 0.3
_W_LINEAR = 20.0


def _optimize_ring_system(
    mol, system_atoms, all_rings, bond_dihedral, coords, placed, parent, ftol=1e-7, gtol=1e-5
):
    """Place and optimize ring system atoms in Cartesian space.

    Initializes ring atoms as regular polygons, then minimizes a cost
    function enforcing ideal bonds, angles, planarity, chirality, and
    AMSR dihedral restraints.
    """
    from scipy.optimize import minimize

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
    chiral_info, chiral_target_vols = _collect_chiral_atoms(mol, sys_set, fixed, coords=coords)
    # Store optimizer chirality signs and target volumes for z-matrix
    # chain placement.  Only core atoms get reliable values (embedding
    # chirality may be wrong for non-core atoms).
    if not hasattr(mol, "_optimized_chiral_sign"):
        mol._optimized_chiral_sign = {}
    if not hasattr(mol, "_optimized_chiral_vol"):
        mol._optimized_chiral_vol = {}
    for idx, row in enumerate(chiral_info):
        mol._optimized_chiral_sign[int(row[0])] = int(row[4])
        mol._optimized_chiral_vol[int(row[0])] = float(chiral_target_vols[idx])
    dih_quads, dih_targets = _collect_ring_dihedrals(mol, sys_set, bond_dihedral, fixed)

    linear_triples = _collect_linear_atoms(mol, sys_set, fixed)
    ez_quads, ez_targets = _collect_ez_constraints(mol, sys_set, fixed)
    rp_quads, rp_targets = _collect_ring_planarity_dihedrals(mol, sys_set, fixed)
    # Merge ring planarity dihedrals with AMSR dihedrals (similar weight)
    if len(rp_quads):
        dih_quads = np.concatenate([dih_quads, rp_quads]) if len(dih_quads) else rp_quads
        dih_targets = np.concatenate([dih_targets, rp_targets]) if len(dih_targets) else rp_targets

    # Build unified atom-index -> slot-index mapping for numpy fancy indexing.
    # Free atoms map to slots 0..n_sys-1; fixed atoms get slots n_sys..n_sys+n_fixed-1.
    max_atom = max(max(sys_list), max(fixed.keys()) if fixed else 0) + 1
    atom_to_slot = np.full(max_atom, -1, dtype=int)
    for atom, i in idx_map.items():
        atom_to_slot[atom] = i
    fixed_list = sorted(fixed.keys())
    fixed_coords = np.zeros((len(fixed_list), 3))
    for fi, fa in enumerate(fixed_list):
        atom_to_slot[fa] = n_sys + fi
        fixed_coords[fi] = fixed[fa]

    # Remap all constraint atom indices to slot indices for direct numpy indexing.
    def _remap(arr):
        return atom_to_slot[arr] if len(arr) else arr

    if len(bonds):
        bonds = _remap(bonds)
    if len(angle_triples):
        angle_triples = _remap(angle_triples)
    if len(planar_groups):
        planar_groups[:, :4] = _remap(planar_groups[:, :4])
    if len(chiral_info):
        chiral_info[:, :4] = _remap(chiral_info[:, :4])
    if len(dih_quads):
        dih_quads = _remap(dih_quads)
    if len(ez_quads):
        ez_quads = _remap(ez_quads)
    if len(linear_triples):
        linear_triples = _remap(linear_triples)

    from .cost_grad import CostGradProblem

    _lin = linear_triples if len(linear_triples) else None
    _objective = CostGradProblem(
        n_sys,
        fixed_coords,
        bonds,
        ideal_lengths,
        angle_triples,
        ideal_angles,
        planar_groups,
        chiral_info,
        chiral_target_vols,
        dih_quads,
        dih_targets,
        ez_quads,
        ez_targets,
        _W_BOND,
        _W_ANGLE,
        _W_PLANAR,
        _W_CHIRAL,
        _W_DIHEDRAL,
        _W_EZ,
        linear_triples=_lin,
        w_linear=_W_LINEAR,
    )

    # Initial coordinates come from the embedding already stored in coords.
    x0 = np.zeros(3 * n_sys)
    for a in sys_list:
        k = idx_map[a]
        x0[3 * k : 3 * k + 3] = coords[a]

    result = minimize(
        _objective, x0, method="L-BFGS-B", jac=True, options={"ftol": ftol, "gtol": gtol}
    )
    best_cost = result.fun
    best_x = result.x

    # If AMSR dihedral targets exist, try ring-inverted starting points.
    # For non-planar rings the optimizer can converge to the mirror-image
    # chair; inverting through the mean plane and re-optimizing often fixes it.
    if best_cost > 0.01 and len(dih_quads):
        # Global inversion (all atoms through overall mean plane)
        x_inv = best_x.copy().reshape(-1, 3)
        centroid = x_inv.mean(axis=0)
        centered = x_inv - centroid
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        normal = Vt[-1]
        for k in range(len(x_inv)):
            d = np.dot(x_inv[k] - centroid, normal)
            x_inv[k] -= 2.0 * d * normal
        r = minimize(
            _objective,
            x_inv.ravel(),
            method="L-BFGS-B",
            jac=True,
            options={"ftol": ftol, "gtol": gtol},
        )
        if r.fun < best_cost:
            best_cost = r.fun
            best_x = r.x

        # Per-ring inversions: invert each non-planar ring individually.
        # This handles cases where only one ring in a multi-ring system
        # needs to be flipped (e.g. cyclohexane chair in a mixed system).
        for ring in all_rings:
            ring_slots = [idx_map[a] for a in ring if a in idx_map]
            if len(ring_slots) < 4:
                continue
            if all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring):
                continue
            x_inv = best_x.copy().reshape(-1, 3)
            rcoords = x_inv[ring_slots]
            rc = rcoords.mean(axis=0)
            _, _, Vt = np.linalg.svd(rcoords - rc, full_matrices=False)
            rn = Vt[-1]
            for k in ring_slots:
                d = np.dot(x_inv[k] - rc, rn)
                x_inv[k] -= 2.0 * d * rn
            r = minimize(
                _objective,
                x_inv.ravel(),
                method="L-BFGS-B",
                jac=True,
                options={"ftol": ftol, "gtol": gtol},
            )
            if r.fun < best_cost:
                best_cost = r.fun
                best_x = r.x

    # Copy optimized coordinates back
    x_opt = best_x.reshape(-1, 3)
    for a in sys_list:
        coords[a] = x_opt[idx_map[a]]
    return best_cost


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
                        coords[i] = place_atom(
                            coords[k], coords[c1], coords[p], bond_len, ang, omega
                        )
                    else:
                        # Try +120 and -120 from ref k; pick the one farthest
                        # from all placed neighbors (avoids atom-ordering-dependent
                        # GetNeighbors() indexing for chirality).
                        best_pos = None
                        best_min_dist = -1.0
                        for omega in [120.0, -120.0]:
                            trial = place_atom(
                                coords[k], coords[c1], coords[p], bond_len, ang, omega
                            )
                            min_d = min(_norm3(trial - coords[nb]) for nb in other_placed)
                            if min_d > best_min_dist:
                                best_min_dist = min_d
                                best_pos = trial
                        coords[i] = best_pos
                else:
                    ref = _synthetic_ref(coords, p, c1)
                    omega = 180.0 if hyb == SP2 else 120.0 * child_count[p]
                    coords[i] = place_atom(ref, coords[c1], coords[p], bond_len, ang, omega)
    else:
        gg = parent[g]
        bond_angle = _get_bond_angle(mol, g, p, i)

        # When g is SP the ref chain (gg-g-p) is collinear and the
        # torsion is undefined.  Place using non-SP placed neighbors.
        if mol.GetAtomWithIdx(g).GetHybridization() == SP:
            non_sp = [
                nb.GetIdx()
                for nb in mol.GetAtomWithIdx(p).GetNeighbors()
                if nb.GetIdx() != i
                and _is_placed(coords, nb.GetIdx())
                and mol.GetAtomWithIdx(nb.GetIdx()).GetHybridization() != SP
            ]
            if len(non_sp) >= 2:
                k1, k2 = non_sp[0], non_sp[1]
                hyb = mol.GetAtomWithIdx(p).GetHybridization()
                if hyb == SP2:
                    coords[i] = place_atom(
                        coords[k2], coords[k1], coords[p], bond_len, bond_angle, 180.0
                    )
                else:
                    best_pos, best_min = None, -1.0
                    all_placed = [
                        nb.GetIdx()
                        for nb in mol.GetAtomWithIdx(p).GetNeighbors()
                        if nb.GetIdx() != i and _is_placed(coords, nb.GetIdx())
                    ]
                    for omega in (120.0, -120.0):
                        trial = place_atom(
                            coords[k2], coords[k1], coords[p], bond_len, bond_angle, omega
                        )
                        min_d = min(_norm3(trial - coords[nb]) for nb in all_placed)
                        if min_d > best_min:
                            best_min = min_d
                            best_pos = trial
                    coords[i] = best_pos
            elif len(non_sp) == 1:
                k = non_sp[0]
                ref = _synthetic_ref(coords, p, k)
                coords[i] = place_atom(ref, coords[k], coords[p], bond_len, bond_angle, 120.0)
            else:
                # All neighbors are SP — just use synthetic ref from g
                ref = _synthetic_ref(coords, g, p)
                coords[i] = place_atom(ref, coords[g], coords[p], bond_len, bond_angle, 120.0)
        else:
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

    # After placing a child of a chiral atom, check volume to catch sign
    # convention mismatches.  If wrong, reflect the new atom across the
    # plane of the parent's other neighbors.
    atom_p = mol.GetAtomWithIdx(p)
    if atom_p.GetChiralTag() in (CW, CCW):
        nbrs = [nb.GetIdx() for nb in atom_p.GetNeighbors()]
        if all(_is_placed(coords, nb) for nb in nbrs) and len(nbrs) >= 3:
            rp = coords[p]
            vs = [coords[nb] - rp for nb in nbrs[:3]]
            vol = np.dot(vs[0], np.cross(vs[1], vs[2]))
            expected_sign = getattr(mol, "_optimized_chiral_sign", {}).get(p)
            if expected_sign is None:
                expected_sign = -1.0 if atom_p.GetChiralTag() == CW else 1.0
            if np.sign(vol) != expected_sign:
                # Reflect atom i across the plane of p's other placed neighbors
                others = [nb for nb in nbrs if nb != i and _is_placed(coords, nb)]
                if len(others) >= 2:
                    v1 = coords[others[0]] - rp
                    v2 = coords[others[1]] - rp
                    normal = np.cross(v1, v2)
                    nn = _norm3(normal)
                    if nn > 1e-10:
                        normal /= nn
                        d = coords[i] - rp
                        coords[i] = rp + d - 2.0 * np.dot(d, normal) * normal
                        # Update first_child_torsion so subsequent children
                        # use the post-reflection base angle.
                        if first_child_idx[p] == i and g is not None:
                            std_ref = coords[gg] if gg is not None else _synthetic_ref(coords, g, p)
                            first_child_torsion[p] = measure_torsion(
                                std_ref, coords[g], coords[p], coords[i]
                            )

    return alternatives


# ============================================================
# Public API
# ============================================================


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
    ftol: float = 1e-7,
    gtol: float = 1e-5,
    max_confs: int = 10,
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

    # Convert quad-keyed dihedrals to bond-keyed for fast lookup.
    bond_dihedral: dict[tuple[int, int], tuple[int, int, int]] = {}
    if dihedral:
        for (mi, i, j, mj), angle in dihedral.items():
            bond_dihedral[(i, j)] = (mi, mj, angle)
            bond_dihedral[(j, i)] = (mj, mi, angle)

    # Default parent tree (used for z-matrix chain placement).
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

    # --- Phase 1: Cartesian-optimize the ring core ---
    core_atoms: set[int] = set()
    if ring_atoms:
        core_atoms = _find_core_atoms(mol, ring_atoms, n)
        # Generate embeddings lazily (one at a time).  Stop when the cost
        # is low and additional embeddings aren't improving it.
        best_opt_cost = float("inf")
        best_coords = None
        best_bd: dict[tuple[int, int], tuple[int, int, int]] = {}
        best_chiral_sign = {}
        best_chiral_vol = {}
        bd_saved = dict(bond_dihedral)
        stale = 0  # consecutive embeddings without >10% improvement
        for attempt in range(max_confs):
            ec_list = _rdkit_embed(mol, n_confs=1, seed=42 + attempt)
            if not ec_list:
                continue
            ec = ec_list[0]
            bond_dihedral.update(bd_saved)
            coords[:] = ec
            _fix_equivalent_terminals(mol, bond_dihedral, coords)
            _correct_junction_dihedrals(mol, ring_systems, core_atoms, bond_dihedral, coords)
            fixed_for_opt = {i: coords[i].copy() for i in range(n) if i not in core_atoms}
            oc = _optimize_ring_system(
                mol,
                core_atoms,
                all_rings,
                bond_dihedral,
                coords,
                fixed_for_opt,
                parent,
                ftol=ftol,
                gtol=gtol,
            )
            improved = oc < best_opt_cost * 0.9
            if oc < best_opt_cost:
                best_opt_cost = oc
                best_coords = coords.copy()
                best_bd = dict(bond_dihedral)
                best_chiral_sign = dict(getattr(mol, "_optimized_chiral_sign", {}))
                best_chiral_vol = dict(getattr(mol, "_optimized_chiral_vol", {}))
            stale = 0 if improved else stale + 1
            # Stop early only when cost has converged to a low value.
            if stale >= 2 and best_opt_cost < 1.0:
                break
        if best_coords is not None:
            coords[:] = best_coords
            bond_dihedral.update(best_bd)
            mol._optimized_chiral_sign = best_chiral_sign
            mol._optimized_chiral_vol = best_chiral_vol

    # --- Phase 2: z-matrix place branch atoms outward from core ---
    placed: set[int] = set(core_atoms)
    if not ring_atoms:
        placed.add(0)
    outward_parent, branch_order = _build_outward_tree(mol, placed, n, ring_atoms)

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
        if p is not None and _has_collision(mol, i, coords, placed, threshold=1.0):
            resolved = False
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
                    resolved = True
                    break
            # When torsion alternatives fail and the parent is SP3 with
            # 3+ placed neighbors, the remaining tetrahedral vertex is
            # uniquely determined.
            if not resolved and mol.GetAtomWithIdx(p).GetHybridization() == SP3:
                nbrs_p = [nb.GetIdx() for nb in mol.GetAtomWithIdx(p).GetNeighbors()]
                placed_nbrs = [nb for nb in nbrs_p if nb != i and _is_placed(coords, nb)]
                if len(placed_nbrs) >= 3:
                    rp = coords[p]
                    vsum = np.zeros(3)
                    for nb in placed_nbrs[:3]:
                        v = coords[nb] - rp
                        nv = _norm3(v)
                        if nv > 1e-10:
                            vsum += v / nv
                    direction = -vsum
                    n_dir = _norm3(direction)
                    if n_dir > 1e-10:
                        coords[i] = rp + _get_bond_length(mol, p, i) * direction / n_dir

    # --- Phase 3: correct any unsatisfied AMSR dihedrals ---
    # Some dihedrals (e.g. on forward bonds not consumed during z-matrix
    # placement) may not have been applied.  Rotate subtrees to fix them.

    for (mi, i, j, mj), angle in (dihedral or {}).items():
        # Skip bonds where an SP atom makes the torsion undefined.
        if mol.GetAtomWithIdx(i).GetHybridization() == SP:
            continue
        if mol.GetAtomWithIdx(j).GetHybridization() == SP:
            continue
        actual = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
        diff = abs((actual - angle + 180) % 360 - 180)
        if diff < 5.0:
            continue
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is None or bond.IsInRing():
            continue
        to_rotate = _bfs_subtree(mol, j, i)
        if mi in to_rotate or mj not in to_rotate:
            continue
        _set_dihedral(coords, mi, i, j, mj, angle, to_rotate)

    # Build RDKit conformer
    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()
