"""Conformer generation from AMSR dihedrals.

Embed the entire molecule with RDKit distance geometry, then optimize
all atoms with a Cartesian-space cost function enforcing ideal bond
lengths, bond angles, planarity at SP2 centers, chirality at SP3
centers, E/Z constraints, and AMSR dihedral restraints.
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

SP = Chem.HybridizationType.SP
SP2 = Chem.HybridizationType.SP2
SP3 = Chem.HybridizationType.SP3
CW = Chem.ChiralType.CHI_TETRAHEDRAL_CW
CCW = Chem.ChiralType.CHI_TETRAHEDRAL_CCW


# ---------------------------------------------------------------------------
# Geometry primitives
# ---------------------------------------------------------------------------


def _cross3(a, b):
    """Cross product for 3-element arrays (avoids numpy.cross overhead)."""
    return np.array(
        [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
    )


def _norm3(v):
    """Euclidean norm for 3-element array (avoids numpy.linalg.norm overhead)."""
    return np.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def measure_torsion(p0, p1, p2, p3):
    """Torsion angle (degrees) for four 3-D points."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = _cross3(b1, b2), _cross3(b2, b3)
    n1n, n2n = _norm3(n1), _norm3(n2)
    if n1n < 1e-10 or n2n < 1e-10:
        return 0.0
    n1, n2 = n1 / n1n, n2 / n2n
    return np.degrees(np.arctan2(np.dot(_cross3(n1, n2), b2 / _norm3(b2)), np.dot(n1, n2)))


# ---------------------------------------------------------------------------
# RDKit helpers
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
    # Single bonds between two SP2 atoms have resonance character
    # (e.g. amide C-N, ester/carbamate C-O).  Use the average of
    # single and aromatic bond lengths.
    if bo == 1:
        ai, aj = mol.GetAtomWithIdx(i), mol.GetAtomWithIdx(j)
        if ai.GetHybridization() == SP2 and aj.GetHybridization() == SP2:
            pair = (min(s1, s2), max(s1, s2))
            single = _BOND_LENGTHS.get((*pair, 1))
            arom = _BOND_LENGTHS.get((*pair, 1.5))
            if single is not None and arom is not None:
                return 0.5 * (single + arom)
    key = (min(s1, s2), max(s1, s2), bo)
    if key in _BOND_LENGTHS:
        return _BOND_LENGTHS[key]
    return _COVALENT_RADII.get(s1, 1.5) + _COVALENT_RADII.get(s2, 1.5)


_ELEMENT_ANGLES = {
    "S": {SP2: 102.0, SP3: 96.0},
    "Se": {SP2: 100.0, SP3: 95.0},
}


def _get_bond_angle(mol, a, b, c):
    """Ideal bond angle a-b-c in degrees."""
    atom_b = mol.GetAtomWithIdx(b)
    hyb = atom_b.GetHybridization()
    sym = atom_b.GetSymbol()
    if sym in _ELEMENT_ANGLES and hyb in _ELEMENT_ANGLES[sym]:
        # Reduced angles only apply to low-coordination heteroatoms
        # (e.g. thioether R-S-R).  High-coordination centers like
        # sulfonyl S(=O)2 are approximately tetrahedral.
        if atom_b.GetDegree() <= 2:
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
            if hyb == SP3 and n <= 6 and poly < 109.5:
                for ring in ri.AtomRings():
                    if len(ring) == n and b in ring:
                        if any(
                            mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in ring if x != b
                        ):
                            return poly
                        break
    return _HYBRID_ANGLES.get(hyb, 109.5)


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


def _subtree(mol, root, exclude):
    """BFS to find all atoms reachable from root without crossing exclude."""
    visited = {root}
    queue = [root]
    while queue:
        a = queue.pop(0)
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in visited and ni != exclude:
                visited.add(ni)
                queue.append(ni)
    return visited


def _rotate_subtree(coords, atoms, origin, axis_dir, angle_deg):
    """Rotate atoms around an axis through origin by angle_deg (Rodrigues)."""
    rad = np.radians(angle_deg)
    axis = axis_dir / _norm3(axis_dir)
    cos_a, sin_a = np.cos(rad), np.sin(rad)
    for k in atoms:
        v = coords[k] - origin
        coords[k] = (
            origin + v * cos_a + _cross3(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)
        )


def _set_dihedrals(mol, bond_dihedral, coords):
    """Rotate subtrees around non-ring bonds to match AMSR dihedral targets."""
    seen = set()
    for (i, j), (mi, mj, target) in bond_dihedral.items():
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is not None and bond.IsInRing():
            continue
        current = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
        delta = (target - current + 180) % 360 - 180
        if abs(delta) < 1.0:
            continue
        tree_j = _subtree(mol, j, i)
        tree_i = _subtree(mol, i, j)
        if len(tree_j) <= len(tree_i):
            _rotate_subtree(coords, tree_j - {j}, coords[i], coords[j] - coords[i], delta)
        else:
            _rotate_subtree(coords, tree_i - {i}, coords[j], coords[i] - coords[j], -delta)


def _fix_ring_puckers(mol, bond_dihedral, coords):
    """Reconstruct ring atom positions from AMSR ring dihedrals.

    For each ring with encoded dihedral targets, build coordinates from
    internal geometry (bond lengths, angles, dihedrals), Kabsch-align
    onto the embedding, and translate attached substituents.
    """
    all_rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    # Only fix isolated rings (no atoms shared with other rings).
    # Fused ring systems have coupled geometry — let the optimizer handle them.
    all_ring_atoms = {}
    for idx, ring in enumerate(all_rings):
        for a in ring:
            all_ring_atoms.setdefault(a, set()).add(idx)
    shared_atoms = {a for a, rings in all_ring_atoms.items() if len(rings) > 1}

    for ring in all_rings:
        ring_set = set(ring)
        if ring_set & shared_atoms:
            continue
        n_ring = len(ring)
        if n_ring < 4:
            continue

        # Collect target dihedrals for consecutive ring quadruples
        ring_dihedrals: dict[int, float] = {}
        for idx in range(n_ring):
            b, c = ring[(idx + 1) % n_ring], ring[(idx + 2) % n_ring]
            key_bc = (min(b, c), max(b, c))
            for (bi, bj), (mi, mj, target) in bond_dihedral.items():
                if (min(bi, bj), max(bi, bj)) == key_bc and mi in ring_set and mj in ring_set:
                    ring_dihedrals[idx] = target
                    break
        if len(ring_dihedrals) < n_ring:
            continue

        # Bond lengths and angles around the ring
        lengths = [_get_bond_length(mol, ring[k], ring[(k + 1) % n_ring]) for k in range(n_ring)]
        angles = [
            _get_bond_angle(mol, ring[(k - 1) % n_ring], ring[k], ring[(k + 1) % n_ring])
            for k in range(n_ring)
        ]

        # Build ring from internal coordinates (z-matrix style)
        new_pos = np.zeros((n_ring, 3))
        new_pos[0] = [0, 0, 0]
        new_pos[1] = [lengths[0], 0, 0]
        if n_ring >= 3:
            a = np.radians(180 - angles[1])
            new_pos[2] = [lengths[0] - lengths[1] * np.cos(a), lengths[1] * np.sin(a), 0]
        for idx in range(3, n_ring):
            p1, p2, p3 = new_pos[idx - 3], new_pos[idx - 2], new_pos[idx - 1]
            d_len = lengths[idx - 1]
            theta = np.radians(180 - angles[idx])
            tau = np.radians(ring_dihedrals.get((idx - 3) % n_ring, 0))
            bc = p3 - p2
            bc_n = bc / _norm3(bc)
            ab = p2 - p1
            n_vec = _cross3(ab, bc)
            nn = _norm3(n_vec)
            n_vec = n_vec / nn if nn > 1e-10 else np.array([0.0, 0.0, 1.0])
            m_vec = _cross3(n_vec, bc_n)
            new_pos[idx] = p3 + d_len * (
                bc_n * np.cos(theta)
                + m_vec * np.sin(theta) * np.cos(tau)
                + n_vec * np.sin(theta) * np.sin(tau)
            )

        # Kabsch alignment onto embedded ring positions
        old_pos = coords[ring]
        old_center = old_pos.mean(axis=0)
        new_center = new_pos.mean(axis=0)
        H = (new_pos - new_center).T @ (old_pos - old_center)
        U, _, Vt = np.linalg.svd(H)
        d = np.linalg.det(Vt.T @ U.T)
        R = Vt.T @ np.diag([1, 1, np.sign(d)]) @ U.T
        aligned = (new_pos - new_center) @ R.T + old_center

        # Move ring atoms and translate their substituent subtrees
        for idx, ai in enumerate(ring):
            disp = aligned[idx] - coords[ai]
            for nb in mol.GetAtomWithIdx(ai).GetNeighbors():
                ni = nb.GetIdx()
                if ni in ring_set:
                    continue
                for k in _subtree(mol, ni, ai):
                    coords[k] += disp
            coords[ai] = aligned[idx]


# ============================================================
# Constraint collection
# ============================================================


def _to_array(lst, dtype=int, cols=None):
    """Convert list to numpy array, returning an appropriately shaped empty array if empty."""
    if lst:
        return np.array(lst, dtype=dtype)
    return np.empty((0, cols) if cols else (0,), dtype=dtype)


def _collect_bonds(mol, atoms):
    """Collect all bonds between atoms in the set."""
    pairs, ideals = [], []
    for a in sorted(atoms):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b in atoms and b > a:
                pairs.append((a, b))
                ideals.append(_get_bond_length(mol, a, b))
    return _to_array(pairs, cols=2), np.array(ideals) if ideals else np.empty(0)


def _collect_angles(mol, atoms):
    """Collect all bond angle triples.

    For SP2 atoms with exactly 3 angles, adjusts targets so they sum to
    360° (important for fused ring junctions).
    """
    triples, ideals = [], []
    center_indices: dict[int, list[int]] = {}
    for b in sorted(atoms):
        if mol.GetAtomWithIdx(b).GetHybridization() == SP:
            continue
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(b).GetNeighbors() if nb.GetIdx() in atoms]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                idx = len(triples)
                triples.append((a, b, c))
                ideals.append(_get_bond_angle(mol, a, b, c))
                center_indices.setdefault(b, []).append(idx)
    ri = mol.GetRingInfo()
    for b, indices in center_indices.items():
        if len(indices) != 3:
            continue
        if mol.GetAtomWithIdx(b).GetHybridization() != SP2:
            continue
        total = sum(ideals[k] for k in indices)
        if abs(total - 360.0) < 1.0:
            continue
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


def _collect_planar_atoms(mol, atoms):
    """Collect planarity constraints for SP2 atoms with 3+ neighbors."""
    groups = []
    for j in sorted(atoms):
        atom = mol.GetAtomWithIdx(j)
        if atom.GetHybridization() != SP2:
            continue
        if atom.GetChiralTag() in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in atoms]
        if len(nbrs) >= 3:
            groups.append((j, nbrs[0], nbrs[1], nbrs[2]))
    return _to_array(groups, cols=4)


def _collect_linear_atoms(mol, atoms):
    """Collect linearity constraints for SP atoms."""
    triples = []
    for j in sorted(atoms):
        if mol.GetAtomWithIdx(j).GetHybridization() != SP:
            continue
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() in atoms]
        if len(nbrs) == 2:
            triples.append((nbrs[0], j, nbrs[1]))
    return _to_array(triples, cols=3)


def _collect_chiral_atoms(mol, atoms, coords):
    """Collect chirality constraints for SP3 chiral atoms.

    Uses the embedding volume sign when clearly pyramidal; otherwise
    falls back to the CW/CCW tag.
    """
    result = []
    target_vols = []
    for j in sorted(atoms):
        atom = mol.GetAtomWithIdx(j)
        chiral = atom.GetChiralTag()
        if chiral not in (CW, CCW):
            continue
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in atoms]
        if len(nbrs) < 3:
            continue
        rj = coords[j]
        positions = []
        for nb in nbrs[:3]:
            if mol.GetAtomWithIdx(nb).GetHybridization() == SP:
                sp_other = [
                    x.GetIdx() for x in mol.GetAtomWithIdx(nb).GetNeighbors() if x.GetIdx() != j
                ]
                if sp_other:
                    far = coords[sp_other[0]]
                    direction = far - rj
                    dn = _norm3(direction)
                    if dn > 1e-10:
                        direction /= dn
                        d = _get_bond_length(mol, j, nb)
                        positions.append(rj + direction * d)
                    else:
                        positions.append(coords[nb])
                else:
                    positions.append(coords[nb])
            else:
                positions.append(coords[nb])
        v1, v2, v3 = positions[0] - rj, positions[1] - rj, positions[2] - rj
        vol = np.dot(v1, np.cross(v2, v3))
        denom = _norm3(v1) * _norm3(v2) * _norm3(v3)
        oop = abs(vol) / denom if denom > 1e-10 else 0.0
        has_sp_nbr = any(mol.GetAtomWithIdx(nb).GetHybridization() == SP for nb in nbrs[:3])
        if oop > 0.1 and not has_sp_nbr:
            sign = 1 if vol > 0 else -1
        else:
            sign = -1 if chiral == CW else 1
        target_vols.append(sign * 2.5)
        result.append((j, nbrs[0], nbrs[1], nbrs[2], sign))
    return _to_array(result, cols=5), _to_array(target_vols, dtype=float)


def _collect_dihedral_restraints(mol, atoms, bond_dihedral):
    """Collect AMSR dihedral restraints for all bonds."""
    quads, targets = [], []
    seen = set()
    for a in sorted(atoms):
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            b = nb.GetIdx()
            if b not in atoms:
                continue
            key = (min(a, b), max(a, b))
            if key in seen:
                continue
            seen.add(key)
            bond_key = (a, b) if (a, b) in bond_dihedral else (b, a)
            if bond_key in bond_dihedral:
                mi, mj, angle = bond_dihedral[bond_key]
                if mi in atoms and mj in atoms:
                    i, j = bond_key
                    if mol.GetAtomWithIdx(i).GetHybridization() == SP:
                        continue
                    if mol.GetAtomWithIdx(j).GetHybridization() == SP:
                        continue
                    quads.append((mi, i, j, mj))
                    targets.append(float(angle))
    return _to_array(quads, cols=4), _to_array(targets, dtype=float)


def _collect_planarity_dihedrals(mol, atoms):
    """Collect 0° torsion constraints for planar ring bonds."""
    mol_k = Chem.RWMol(mol)
    try:
        Chem.Kekulize(mol_k, clearAromaticFlags=False)
    except Exception:
        mol_k = mol
    quads, targets = [], []
    seen = set()
    for ring in mol.GetRingInfo().AtomRings():
        if not all(a in atoms for a in ring):
            continue
        n = len(ring)
        all_sp2 = all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring)
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            if not all(mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in (a, b, c, d)):
                continue
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


def _collect_ez_constraints(mol, atoms):
    """Collect cis/trans dihedral constraints for E/Z double bonds."""
    quads, targets = [], []
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i not in atoms and j not in atoms:
            continue
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) < 2:
            continue
        si, sj = stereo_atoms[0], stereo_atoms[1]
        if si not in atoms or sj not in atoms:
            continue
        target = 0.0 if stereo == Chem.BondStereo.STEREOZ else 180.0
        quads.append((si, i, j, sj))
        targets.append(target)
    return _to_array(quads, cols=4), _to_array(targets, dtype=float)


# ============================================================
# RDKit embedding
# ============================================================


def _rdkit_embed(mol, n_confs=1, seed=42):
    """Embed molecule with RDKit distance geometry.

    Returns list of Nx3 coordinate arrays (one per conformer), or empty list
    on failure.  Heavy-atom indices match the input mol.
    """
    from rdkit.Chem import AllChem

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
# Optimizer
# ============================================================

_W_BOND = 5.0
_W_ANGLE = 2.0
_W_PLANAR = 3.0
_W_CHIRAL = 10.0
_W_DIHEDRAL = 0.1
_W_EZ = 0.3
_W_LINEAR = 20.0


def _optimize(mol, bond_dihedral, coords, ftol=1e-3, gtol=1e-1):
    """Optimize all atom positions to satisfy geometry constraints.

    Returns the optimizer cost.  Mutates coords in place.
    """
    from scipy.optimize import minimize

    n = mol.GetNumAtoms()
    atoms = set(range(n))

    bonds, ideal_lengths = _collect_bonds(mol, atoms)
    angle_triples, ideal_angles = _collect_angles(mol, atoms)
    planar_groups = _collect_planar_atoms(mol, atoms)
    chiral_info, chiral_target_vols = _collect_chiral_atoms(mol, atoms, coords)
    dih_quads, dih_targets = _collect_dihedral_restraints(mol, atoms, bond_dihedral)
    linear_triples = _collect_linear_atoms(mol, atoms)
    ez_quads, ez_targets = _collect_ez_constraints(mol, atoms)
    rp_quads, rp_targets = _collect_planarity_dihedrals(mol, atoms)

    if len(rp_quads):
        dih_quads = np.concatenate([dih_quads, rp_quads]) if len(dih_quads) else rp_quads
        dih_targets = np.concatenate([dih_targets, rp_targets]) if len(dih_targets) else rp_targets

    from .cost_grad import CostGradProblem

    _objective = CostGradProblem(
        n,
        np.zeros((0, 3)),
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
        linear_triples=linear_triples if len(linear_triples) else None,
        w_linear=_W_LINEAR,
    )

    result = minimize(
        _objective,
        coords.ravel(),
        method="L-BFGS-B",
        jac=True,
        options={"ftol": ftol, "gtol": gtol},
    )
    best_cost = result.fun
    best_x = result.x

    # Try ring-inverted starting points.  For non-planar rings the
    # optimizer can converge to the wrong chair; inverting through the
    # mean plane and re-optimizing often fixes it.
    all_rings = [tuple(r) for r in mol.GetRingInfo().AtomRings()]
    if best_cost > 1.0 and len(dih_quads):
        # Global inversion
        x_inv = best_x.copy().reshape(-1, 3)
        centroid = x_inv.mean(axis=0)
        centered = x_inv - centroid
        _, _, Vt = np.linalg.svd(centered, full_matrices=False)
        normal = Vt[-1]
        for k in range(n):
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

        # Per-ring inversions
        for ring in all_rings:
            if len(ring) < 4:
                continue
            if all(mol.GetAtomWithIdx(a).GetHybridization() == SP2 for a in ring):
                continue
            x_inv = best_x.copy().reshape(-1, 3)
            rcoords = x_inv[list(ring)]
            rc = rcoords.mean(axis=0)
            _, _, Vt = np.linalg.svd(rcoords - rc, full_matrices=False)
            rn = Vt[-1]
            for k in ring:
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

        # Per-chiral-center inversions: reflect all neighbors of a
        # chiral atom through the plane of its three placed neighbors.
        # The embedding may assign the wrong handedness to pseudo-chiral
        # centers (e.g. C(OH)(Ph)(Ph) with two identical phenyls).
        for row in chiral_info:
            center = int(row[0])
            nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(center).GetNeighbors()]
            if len(nbrs) < 3:
                continue
            x_inv = best_x.copy().reshape(-1, 3)
            rp = x_inv[center]
            v1 = x_inv[nbrs[0]] - rp
            v2 = x_inv[nbrs[1]] - rp
            normal = np.cross(v1, v2)
            nn = _norm3(normal)
            if nn < 1e-10:
                continue
            normal /= nn
            # Reflect all atoms bonded on the "far" side through the plane.
            # Simplest: reflect atom center's non-plane neighbors.
            for nb in nbrs[2:]:
                d = np.dot(x_inv[nb] - rp, normal)
                x_inv[nb] -= 2.0 * d * normal
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

    coords[:] = best_x.reshape(-1, 3)
    return best_cost


# ============================================================
# Public API
# ============================================================


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
    ftol: float = 1e-3,
    gtol: float = 1e-1,
    max_confs: int = 30,
) -> Chem.Mol:
    """Generate 3D conformer.

    1. Embed entire molecule with RDKit distance geometry.
    2. Optimize all atom positions with a cost function enforcing ideal
       bond lengths, angles, planarity, chirality, and AMSR dihedrals.
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

    best_cost = float("inf")
    best_coords = None

    for attempt in range(max_confs):
        ec_list = _rdkit_embed(mol, n_confs=1, seed=42 + attempt)
        if not ec_list:
            continue
        bd = dict(bond_dihedral)  # copy — _fix_equivalent_terminals mutates
        coords[:] = ec_list[0]
        _fix_equivalent_terminals(mol, bd, coords)
        _fix_ring_puckers(mol, bd, coords)
        _set_dihedrals(mol, bd, coords)
        oc = _optimize(mol, bd, coords, ftol=ftol, gtol=gtol)
        if oc < best_cost:
            best_cost = oc
            best_coords = coords.copy()
            if best_cost < 1.0:
                break

    if best_coords is not None:
        coords[:] = best_coords

    conf = Chem.Conformer(n)
    conf.Set3D(True)
    for i in range(n):
        conf.SetAtomPosition(i, coords[i].tolist())
    mol = Chem.RWMol(mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol()


# ============================================================
# Legacy MMFF-based conformer generation
# ============================================================


def GetConformerAndEnergy(
    mol: Chem.Mol, dihedral: Optional[dict[tuple[int, int, int, int], int]] = None
) -> Chem.Mol:
    """Return a conformer using MMFF94 force field.

    :param mol: RDKit Mol
    :param dihedral: optional dictionary of dihedral angle constraints
    :return: (RDKit Mol, energy in kcal/mol)
    """
    from rdkit.Chem import rdMolTransforms as rdMT

    mol = Chem.AddHs(mol)
    Chem.AllChem.EmbedMolecule(mol)
    mp = Chem.AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = Chem.AllChem.MMFFGetMoleculeForceField(mol, mp)
    if dihedral is not None:
        for (i, j, k, last), v in dihedral.items():
            if not mol.GetBondBetweenAtoms(j, k).IsInRing():
                rdMT.SetDihedralDeg(mol.GetConformer(0), i, j, k, last, v)
                ff.MMFFAddTorsionConstraint(i, j, k, last, False, v, v, 1e3)
    ff.Minimize(maxIts=100000)
    rdMT.CanonicalizeConformer(mol.GetConformer())
    ener = ff.CalcEnergy()  # kcal/mol
    return Chem.RemoveHs(mol), ener


def GetRoundedDihedral(mol: Chem.Mol, dihedral: tuple[int, int, int, int], ndeg: int) -> int:
    """Return dihedral angle rounded to nearest ndeg degrees.

    :param mol: RDKit Mol
    :param dihedral: tuple of four atom indices
    :param ndeg: int
    :return: rounded dihedral angle
    """
    from rdkit.Chem import rdMolTransforms as rdMT

    return round(rdMT.GetDihedralDeg(mol.GetConformer(0), *dihedral) / ndeg) * ndeg
