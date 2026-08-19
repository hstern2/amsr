"""Conformer generation from AMSR dihedrals.

Embed with distance geometry using AMSR-encoded dihedrals, then optimize
all atoms with a Cartesian-space cost function enforcing ideal bond
lengths, bond angles, planarity at SP2 centers, chirality at SP3
centers, E/Z constraints, and AMSR dihedral restraints.
"""

import ctypes
from collections import deque
from importlib.machinery import EXTENSION_SUFFIXES
from pathlib import Path
from typing import Optional

import numpy as np
from rdkit import Chem

# ---------------------------------------------------------------------------
# C shared library (conf_util.c): cost_and_grad + embed
# ---------------------------------------------------------------------------

_LIB_DIR = Path(__file__).parent / "src"
_LIB_CANDIDATES = [
    _LIB_DIR / "conf_util.so",
    _LIB_DIR / "conf_util.dylib",
    *(_LIB_DIR / f"conf_util{suffix}" for suffix in EXTENSION_SUFFIXES),
]
_LIB_PATH = next((path for path in _LIB_CANDIDATES if path.is_file()), None)
if _LIB_PATH is None:
    raise ImportError(
        "AMSR conformer library is missing; reinstall AMSR from a source or platform wheel"
    )
_lib = ctypes.CDLL(str(_LIB_PATH))

_c_lbfgs = _lib.lbfgs_optimize
_c_lbfgs.restype = ctypes.c_double
_c_lbfgs.argtypes = [
    ctypes.c_void_p,
    ctypes.c_int,  # x, ndim
    ctypes.c_int,  # n_atoms
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # bonds
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # angles
    ctypes.c_void_p,
    ctypes.c_int,  # planar
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # chiral
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # dihedral
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # ez
    ctypes.c_void_p,
    ctypes.c_int,  # linear
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,  # w_bond, w_angle, w_planar
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_double,  # w_chiral, w_dih, w_ez
    ctypes.c_double,  # w_linear
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,  # ftol, gtol, max_iter
]

_c_embed = _lib.embed
_c_embed.restype = None
_c_embed.argtypes = [
    ctypes.c_int,
    ctypes.c_void_p,  # n, coords_out
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # bonds
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # angles
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # dihedrals
    ctypes.c_uint,  # seed
]

_c_prepare_embed_bounds = _lib.prepare_embed_bounds
_c_prepare_embed_bounds.restype = None
_c_prepare_embed_bounds.argtypes = [
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,  # n, lower, upper, angle_distances
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # bonds
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # angles
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # dihedrals
]

_c_embed_bounds = _lib.embed_bounds
_c_embed_bounds.restype = None
_c_embed_bounds.argtypes = [
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_void_p,  # n, coords_out, lower, upper
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # bonds
    ctypes.c_int,
    ctypes.c_void_p,
    ctypes.c_void_p,  # angles
    ctypes.c_uint,  # seed
]


def _to_int32(arr):
    if isinstance(arr, np.ndarray) and arr.dtype == np.int32 and arr.flags["C_CONTIGUOUS"]:
        return arr
    a = np.asarray(arr, dtype=np.int32)
    return np.ascontiguousarray(a)


def _to_f64(arr):
    if isinstance(arr, np.ndarray) and arr.dtype == np.float64 and arr.flags["C_CONTIGUOUS"]:
        return arr
    a = np.asarray(arr, dtype=np.float64)
    return np.ascontiguousarray(a)


def _dptr(arr):
    """Data pointer for a numpy array, or NULL."""
    return arr.ctypes.data if arr.size else 0


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


def _build_ring_index(mol):
    """Pre-compute per-atom ring membership: atom_idx → list of (ring_size, ring_set, has_sp2)."""
    ri = mol.GetRingInfo()
    ring_data = []
    for ring in ri.AtomRings():
        rs = frozenset(ring)
        n = len(ring)
        has_sp2 = any(mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in ring)
        ring_data.append((n, rs, has_sp2))
    # Map each atom to the rings it belongs to
    atom_rings: dict[int, list[tuple[int, frozenset[int], bool]]] = {}
    for rd in ring_data:
        for a in rd[1]:
            atom_rings.setdefault(a, []).append(rd)
    return atom_rings


def _get_bond_angle(mol, a, b, c, atom_rings=None):
    """Ideal bond angle a-b-c in degrees."""
    atom_b = mol.GetAtomWithIdx(b)
    hyb = atom_b.GetHybridization()
    sym = atom_b.GetSymbol()
    if sym in _ELEMENT_ANGLES and hyb in _ELEMENT_ANGLES[sym]:
        if atom_b.GetDegree() <= 2:
            return _ELEMENT_ANGLES[sym][hyb]
    b_ab = mol.GetBondBetweenAtoms(a, b)
    b_bc = mol.GetBondBetweenAtoms(b, c)
    if b_ab is not None and b_bc is not None:
        # Find the smallest ring containing all three atoms a-b-c.
        best_n = 0
        best_ring = None
        best_has_sp2 = False
        if atom_rings is not None:
            b_rings = atom_rings.get(b)
            if b_rings:
                for n, rs, has_sp2 in b_rings:
                    if a in rs and c in rs and (best_ring is None or n < best_n):
                        best_n = n
                        best_ring = rs
                        best_has_sp2 = has_sp2
        else:
            for ring in mol.GetRingInfo().AtomRings():
                if a in ring and b in ring and c in ring:
                    n = len(ring)
                    if best_ring is None or n < best_n:
                        best_n = n
                        best_ring = ring
                        best_has_sp2 = any(
                            mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in ring if x != b
                        )
        if best_ring is not None:
            poly = (best_n - 2) * 180.0 / best_n
            if hyb in (SP2, Chem.HybridizationType.SP) and best_n <= 6:
                return poly
            if hyb == SP3 and best_n <= 4 and poly < 109.5:
                return poly
            if hyb == SP3 and best_n <= 6 and poly < 109.5:
                if best_has_sp2:
                    return poly
    default = _HYBRID_ANGLES.get(hyb, 109.5)
    # SP2 with degree 2 and one implicit H: the two heavy-atom
    # neighbors occupy the wider angle (~126°) while the H takes ~117°.
    # Exception: if both neighbors are SP2 (conjugated system), use 120°.
    if hyb == SP2 and atom_b.GetDegree() == 2 and atom_b.GetTotalNumHs() > 0:
        na, nc = mol.GetAtomWithIdx(a), mol.GetAtomWithIdx(c)
        if na.GetHybridization() != SP2 or nc.GetHybridization() != SP2:
            return 126.0
    # SP3 with degree 2 and implicit H: heavy-atom angle opens up
    # because H atoms occupy the smaller angular sectors (Bent's rule).
    if hyb == SP3 and atom_b.GetDegree() == 2 and atom_b.GetTotalNumHs() > 0:
        return 112.0
    return default


def _build_pseudo_ez_ops(mol, subtree_cache=None):
    """Pre-compute pseudo-E/Z reflections for _fix_pseudo_ez."""
    ops = []
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) < 2:
            continue
        si, sj = stereo_atoms[0], stereo_atoms[1]
        ni = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() != j]
        nj = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() != i]
        ni_ranks = [ranks[n] for n in ni]
        nj_ranks = [ranks[n] for n in nj]
        if len(set(ni_ranks)) == len(ni_ranks) and len(set(nj_ranks)) == len(nj_ranks):
            continue
        atoms_to_reflect: set[int] = set()
        for nb_idx in nj:
            sub = subtree_cache[(nb_idx, j)] if subtree_cache else _subtree(mol, nb_idx, j)
            atoms_to_reflect.symmetric_difference_update(sub)
        if atoms_to_reflect:
            target = 0.0 if stereo == Chem.BondStereo.STEREOZ else 180.0
            ops.append((si, i, j, sj, target, np.fromiter(atoms_to_reflect, dtype=np.intp)))
    return ops


def _fix_pseudo_ez(mol, coords, subtree_cache=None, pseudo_ez_ops=None):
    """Fix pseudo-E/Z double bonds after embedding.

    For double bonds with graph-equivalent substituents on at least one
    side, check if the E/Z geometry matches the tag and reflect if wrong.
    """
    if pseudo_ez_ops is None:
        pseudo_ez_ops = _build_pseudo_ez_ops(mol, subtree_cache)
    for si, i, j, sj, target, atoms in pseudo_ez_ops:
        actual = measure_torsion(coords[si], coords[i], coords[j], coords[sj])
        diff = abs((actual - target + 180) % 360 - 180)
        if diff < 45:
            continue
        bond_vec = coords[j] - coords[i]
        v_si = coords[si] - coords[i]
        normal = np.cross(bond_vec, v_si)
        nn = _norm3(normal)
        if nn < 1e-10:
            continue
        normal /= nn
        mid = 0.5 * (coords[i] + coords[j])
        _reflect_atoms(coords, atoms, mid, normal)


def _subtree(mol, root, exclude):
    """BFS to find all atoms reachable from root without crossing exclude."""
    visited = {root}
    queue = deque([root])
    while queue:
        a = queue.popleft()
        for nb in mol.GetAtomWithIdx(a).GetNeighbors():
            ni = nb.GetIdx()
            if ni not in visited and ni != exclude:
                visited.add(ni)
                queue.append(ni)
    return visited


def _rotate_subtree(coords, atoms, origin, axis_dir, angle_deg):
    """Rotate atoms around an axis through origin by angle_deg (Rodrigues)."""
    if len(atoms) == 0:
        return
    rad = np.radians(angle_deg)
    axis = axis_dir / _norm3(axis_dir)
    cos_a, sin_a = np.cos(rad), np.sin(rad)
    if len(atoms) >= 4:
        v = coords[atoms] - origin
        coords[atoms] = (
            origin
            + v * cos_a
            + np.cross(axis, v) * sin_a
            + axis * (v @ axis)[:, None] * (1 - cos_a)
        )
        return
    for k in atoms:
        v = coords[k] - origin
        coords[k] = (
            origin + v * cos_a + _cross3(axis, v) * sin_a + axis * np.dot(axis, v) * (1 - cos_a)
        )


def _reflect_atoms(coords, atoms, center, normal):
    """Reflect atoms through the plane defined by center and normal."""
    if len(atoms) >= 4:
        group = coords[atoms]
        coords[atoms] = group - 2.0 * ((group - center) @ normal)[:, None] * normal
        return
    for k in atoms:
        d = np.dot(coords[k] - center, normal)
        coords[k] -= 2.0 * d * normal


def _build_subtree_cache(mol):
    """Pre-compute subtrees for all bonds: cache[(root, exclude)] = set of atoms."""
    cache = {}
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        cache[(i, j)] = _subtree(mol, i, j)
        cache[(j, i)] = _subtree(mol, j, i)
    return cache


def _build_dihedral_ops(mol, bond_dihedral, subtree_cache=None):
    """Pre-compute non-ring dihedral rotations for _set_dihedrals."""
    ops = []
    seen = set()
    for (i, j), (mi, mj, target) in bond_dihedral.items():
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        bond = mol.GetBondBetweenAtoms(i, j)
        if bond is not None and bond.IsInRing():
            continue
        if mol.GetAtomWithIdx(i).GetHybridization() == SP:
            continue
        if mol.GetAtomWithIdx(j).GetHybridization() == SP:
            continue

        tree_j = subtree_cache[(j, i)] if subtree_cache else _subtree(mol, j, i)
        tree_i = subtree_cache[(i, j)] if subtree_cache else _subtree(mol, i, j)
        if len(tree_j) <= len(tree_i) or mi in tree_i:
            atoms = np.fromiter(tree_j - {j}, dtype=np.intp)
            ops.append((mi, i, j, mj, target, atoms, i, j, 1.0))
        else:
            atoms = np.fromiter(tree_i - {i}, dtype=np.intp)
            ops.append((mi, i, j, mj, target, atoms, j, i, -1.0))
    return ops


def _build_chirality_ops(mol, subtree_cache=None):
    """Pre-compute chiral-center reflections for _fix_chirality."""
    ops = []
    n = mol.GetNumAtoms()
    for atom in mol.GetAtoms():
        chiral = atom.GetChiralTag()
        if chiral not in (CW, CCW):
            continue
        j = atom.GetIdx()
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
        if len(nbrs) < 3:
            continue
        nbrs3 = tuple(nbrs[:3])
        expected = -1 if chiral == CW else 1
        best = None
        for ni in nbrs3:
            sub = subtree_cache[(ni, j)] if subtree_cache else _subtree(mol, ni, j)
            if len(sub) <= n // 2 and (best is None or len(sub) < len(best[1])):
                best = (ni, sub)
        if best is None:
            continue
        reflect_nb, reflect_sub = best
        others = tuple(ni for ni in nbrs3 if ni != reflect_nb)
        ops.append(
            (
                j,
                nbrs3,
                expected,
                np.fromiter(reflect_sub, dtype=np.intp),
                others,
            )
        )
    return ops


def _fix_chirality(mol, coords, subtree_cache=None, chirality_ops=None):
    """Fix chiral centers whose signed volume has the wrong sign.

    For each wrong CW/CCW center, find the smallest neighbor subtree
    that is a true branch (< half the molecule) and reflect it through
    the plane of the other two neighbors.  Processes from interior
    to exterior and iterates until stable.
    """
    if chirality_ops is None:
        chirality_ops = _build_chirality_ops(mol, subtree_cache)
    for _iteration in range(10):
        fixed_any = False
        for j, nbrs, expected, reflect_atoms, others in chirality_ops:
            rj = coords[j]
            v = [coords[ni] - rj for ni in nbrs[:3]]
            vol = np.dot(v[0], _cross3(v[1], v[2]))
            if (vol > 0) == (expected > 0):
                continue
            va, vb = coords[others[0]] - rj, coords[others[1]] - rj
            normal = _cross3(va, vb)
            nn = _norm3(normal)
            if nn < 1e-10:
                continue
            normal /= nn
            _reflect_atoms(coords, reflect_atoms, rj, normal)
            fixed_any = True
        if not fixed_any:
            break


def _set_dihedrals(mol, bond_dihedral, coords, subtree_cache=None, dihedral_ops=None):
    """Rotate subtrees around non-ring bonds to match AMSR dihedral targets."""
    if dihedral_ops is None:
        dihedral_ops = _build_dihedral_ops(mol, bond_dihedral, subtree_cache)
    for mi, i, j, mj, target, atoms, origin_idx, axis_end_idx, delta_sign in dihedral_ops:
        current = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
        delta = (target - current + 180) % 360 - 180
        if abs(delta) < 1.0:
            continue
        _rotate_subtree(
            coords,
            atoms,
            coords[origin_idx],
            coords[axis_end_idx] - coords[origin_idx],
            delta_sign * delta,
        )


def _build_ring_topology(mol, subtree_cache=None):
    """Pre-compute ring topology data (invariant across conformer attempts)."""
    all_rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    ring_atom_arrays = [np.asarray(r, dtype=np.intp) for r in all_rings]
    all_ring_atoms = {}
    for idx, ring in enumerate(all_rings):
        for a in ring:
            all_ring_atoms.setdefault(a, set()).add(idx)
    shared_atoms = {a for a, rings in all_ring_atoms.items() if len(rings) > 1}

    ring_sub_by_atom: list[dict[int, list[set[int]]]] = []
    for ring in all_rings:
        ring_set = set(ring)
        per_atom: dict[int, list[set[int]]] = {}
        for a in ring:
            subs: list[set[int]] = []
            for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                ni = nb.GetIdx()
                if ni not in ring_set:
                    sub = subtree_cache[(ni, a)] if subtree_cache else _subtree(mol, ni, a)
                    subs.append(sub)
            if subs:
                per_atom[a] = subs
        ring_sub_by_atom.append(per_atom)

    # Reflection atom sets for pucker flips.  The original operation reflects
    # the ring atoms and then each attached subtree through the same plane.
    # Since a reflection is its own inverse, only atoms reached an odd number
    # of times need to be transformed.
    ring_reflect_atoms: list[np.ndarray] = []
    for ring_idx, per_atom in enumerate(ring_sub_by_atom):
        atoms_to_reflect: set[int] = set()
        for a in all_rings[ring_idx]:
            atoms_to_reflect.symmetric_difference_update((a,))
        for a in all_rings[ring_idx]:
            for sub in per_atom.get(a, ()):
                atoms_to_reflect.symmetric_difference_update(sub)
        ring_reflect_atoms.append(np.fromiter(atoms_to_reflect, dtype=np.intp))

    # Pre-compute chiral info per ring atom (avoid repeated GetChiralTag calls)
    ring_chiral: list[list[tuple[int, list[int], int]]] = []
    for ring in all_rings:
        chiral_atoms = []
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            chiral = atom.GetChiralTag()
            if chiral in (CW, CCW):
                nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
                if len(nbrs) >= 3:
                    expected = -1 if chiral == CW else 1
                    chiral_atoms.append((a, nbrs[:3], expected))
        ring_chiral.append(chiral_atoms)

    return (
        all_rings,
        ring_atom_arrays,
        shared_atoms,
        ring_sub_by_atom,
        ring_reflect_atoms,
        ring_chiral,
    )


def _fix_ring_puckers(mol, bond_dihedral, coords, atom_rings=None, ring_topo=None):
    """Fix ring conformations.

    Two strategies:
    1. For isolated rings with encoded dihedral targets, rebuild from
       internal coordinates (z-matrix style) and Kabsch-align.
    2. For any ring (including fused) with a chiral atom whose volume
       has the wrong sign, reflect the ring through its mean plane.
       Iterate until stable.
    """
    if ring_topo is None:
        ring_topo = _build_ring_topology(mol)
    (
        all_rings,
        ring_atom_arrays,
        shared_atoms,
        ring_sub_by_atom,
        ring_reflect_atoms,
        ring_chiral,
    ) = ring_topo
    if atom_rings is None:
        atom_rings = _build_ring_index(mol)

    # Strategy 2: flip ring puckers where chiral atoms have wrong sign.
    # Only applies to small rings (≤8); reflecting a macrocycle through its
    # mean plane is geometrically ill-defined and destroys fused-ring geometry.
    for _iteration in range(10):
        flipped_any = False
        for ri_idx, ring in enumerate(all_rings):
            if len(ring) > 8:
                continue
            wrong = 0
            for a, nbrs, expected in ring_chiral[ri_idx]:
                rj = coords[a]
                v = [coords[ni] - rj for ni in nbrs]
                vol = np.dot(v[0], _cross3(v[1], v[2]))
                if (vol > 0) != (expected > 0):
                    wrong += 1
            if wrong == 0:
                continue
            # Reflect ring atoms through mean plane
            ring_atoms = ring_atom_arrays[ri_idx]
            ring_coords = coords[ring_atoms]
            center = ring_coords.mean(axis=0)
            centered = ring_coords - center
            _, vecs = np.linalg.eigh(centered.T @ centered)
            normal = vecs[:, 0]
            atoms = ring_reflect_atoms[ri_idx]
            group = coords[atoms]
            coords[atoms] = group - 2.0 * ((group - center) @ normal)[:, None] * normal
            flipped_any = True
        if not flipped_any:
            break

    # Pre-index bond_dihedral by canonical bond key for fast lookup.
    dih_by_bond: dict[tuple[int, int], list[tuple[int, int, float]]] = {}
    for (bi, bj), (mi, mj, target) in bond_dihedral.items():
        key = (min(bi, bj), max(bi, bj))
        dih_by_bond.setdefault(key, []).append((mi, mj, target))

    # Strategy 1: rebuild isolated rings from internal coordinates.
    for ri_idx, ring in enumerate(all_rings):
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
            for mi, mj, target in dih_by_bond.get(key_bc, ()):
                if mi in ring_set and mj in ring_set:
                    ring_dihedrals[idx] = target
                    break
        if len(ring_dihedrals) < n_ring:
            continue

        # Bond lengths and angles around the ring
        lengths = [_get_bond_length(mol, ring[k], ring[(k + 1) % n_ring]) for k in range(n_ring)]
        angles = [
            _get_bond_angle(
                mol, ring[(k - 1) % n_ring], ring[k], ring[(k + 1) % n_ring], atom_rings
            )
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
        sub_map = ring_sub_by_atom[ri_idx]
        for idx, ai in enumerate(ring):
            disp = aligned[idx] - coords[ai]
            for sub in sub_map.get(ai, ()):
                for k in sub:
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


def _collect_angles(mol, atoms, atom_rings=None):
    """Collect all bond angle triples.

    For SP2 atoms with exactly 3 angles, adjusts targets so they sum to
    360° (important for fused ring junctions).
    """
    ri = mol.GetRingInfo()
    if atom_rings is None:
        atom_rings = _build_ring_index(mol)
    # For angles where all three atoms are in the same 3-membered ring,
    # bond lengths alone determine the geometry — skip those angles.
    three_ring_atoms: set[frozenset[int]] = set()
    for n, rs, _ in {rd for rds in atom_rings.values() for rd in rds}:
        if n == 3:
            three_ring_atoms.add(rs)
    triples, ideals = [], []
    center_indices: dict[int, list[int]] = {}
    for b in sorted(atoms):
        if mol.GetAtomWithIdx(b).GetHybridization() == SP:
            continue
        nbrs = [nb.GetIdx() for nb in mol.GetAtomWithIdx(b).GetNeighbors() if nb.GetIdx() in atoms]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                if any(a in r and b in r and c in r for r in three_ring_atoms):
                    continue
                idx = len(triples)
                triples.append((a, b, c))
                ideals.append(_get_bond_angle(mol, a, b, c, atom_rings))
                center_indices.setdefault(b, []).append(idx)
    # For degree-3 SP2 centers, redistribute angles so they sum to 360°
    # (important for fused ring junctions).
    for b, indices in center_indices.items():
        if mol.GetAtomWithIdx(b).GetHybridization() != SP2:
            continue
        if len(indices) != 3:
            continue
        total = sum(ideals[k] for k in indices)
        if abs(total - 360.0) < 1.0:
            continue
        # Angles not constrained by any small ring (≤6) are adjustable.
        # This includes angles with no common ring AND angles whose only
        # common ring is a large ring (macrocycle), so that SP2 junction
        # atoms in small-ring + macrocycle systems sum to 360°.
        adjustable = []
        for k in indices:
            a, _, c = triples[k]
            b_ab = mol.GetBondBetweenAtoms(a, b)
            b_bc = mol.GetBondBetweenAtoms(b, c)
            if b_ab is None or b_bc is None:
                adjustable.append(k)
                continue
            common = set(ri.BondRingSizes(b_ab.GetIdx())) & set(ri.BondRingSizes(b_bc.GetIdx()))
            if not common or min(common) > 6:
                adjustable.append(k)
        if adjustable:
            per = (360.0 - total) / len(adjustable)
            for k in adjustable:
                ideals[k] += per
    return _to_array(triples, cols=3), np.array(ideals) if ideals else np.empty(0)


def _collect_planar_atoms(mol, atoms):
    """Collect planarity constraints for SP2 atoms with 3+ neighbors.

    Also includes sulfonamide N (N bonded to S), which is experimentally
    near-planar regardless of RDKit's hybridization assignment.
    """
    groups = []
    for j in sorted(atoms):
        atom = mol.GetAtomWithIdx(j)
        nbrs = [nb.GetIdx() for nb in atom.GetNeighbors() if nb.GetIdx() in atoms]
        if len(nbrs) < 3:
            continue
        if atom.GetHybridization() == SP2:
            groups.append((j, nbrs[0], nbrs[1], nbrs[2]))
        elif atom.GetSymbol() == "N" and any(
            mol.GetAtomWithIdx(ni).GetSymbol() == "S" for ni in nbrs
        ):
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


_CHIRAL_VOL_REF = 2.5  # target volume for C(SP3) with three C neighbors
_CHIRAL_BL_REF = 1.54**3  # product of C-C bond lengths


def _collect_chiral_atoms(mol, atoms):
    """Collect chirality constraints for chiral atoms.

    Target volume is scaled by the product of ideal bond lengths
    relative to a C(C,C,C) reference.  SP2 atoms with chiral tags
    (e.g. nitrogen pseudo-stereocenters) use a small target volume
    to encode the face preference while remaining nearly planar.
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
        sign = -1 if chiral == CW else 1
        if atom.GetHybridization() == SP2:
            # Sulfonamide N: partially pyramidal, use measured volume.
            # All other SP2 chiral tags (amide N, etc.): skip — planarity
            # constraint handles their geometry.
            if atom.GetSymbol() == "N" and any(
                mol.GetAtomWithIdx(ni).GetSymbol() == "S" for ni in nbrs[:3]
            ):
                vol = 1.7
            else:
                continue
        else:
            d = (
                _get_bond_length(mol, j, nbrs[0])
                * _get_bond_length(mol, j, nbrs[1])
                * _get_bond_length(mol, j, nbrs[2])
            )
            vol = _CHIRAL_VOL_REF * d / _CHIRAL_BL_REF
        target_vols.append(sign * vol)
        result.append((j, nbrs[0], nbrs[1], nbrs[2], sign))
        # For 4-neighbor sp3 centers, constrain the 4th neighbor too.
        # V(n0,n1,n3) has the opposite sign from V(n0,n1,n2) in a valid tetrahedron.
        if len(nbrs) >= 4 and atom.GetHybridization() != SP2:
            target_vols.append(-sign * vol)
            result.append((j, nbrs[0], nbrs[1], nbrs[3], -sign))
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
    """Collect 0° torsion constraints for planar ring bonds.

    Two sources:
    1. Consecutive SP2 quads in aromatic/conjugated rings.
    2. Double bonds in rings whose neighbors aren't all SP2 (e.g. alkenes
       in partially saturated rings).  Each SP2 end gets a constraint
       using its SP2 neighbor on the other side of the double bond.
    """
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
        if n > 7:
            continue
        for i in range(n):
            a, b, c, d = ring[i], ring[(i + 1) % n], ring[(i + 2) % n], ring[(i + 3) % n]
            if not all(mol.GetAtomWithIdx(x).GetHybridization() == SP2 for x in (a, b, c, d)):
                continue
            key = (min(b, c), max(b, c))
            if key not in seen:
                seen.add(key)
                quads.append((a, b, c, d))
                targets.append(0.0)
    # Double bonds in rings: ensure planarity even when neighbors are SP3.
    # Use ring neighbors (cis across the double bond) to form a 0° torsion.
    ri = mol.GetRingInfo()
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE or not bond.IsInRing():
            continue
        b, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if b not in atoms or c not in atoms:
            continue
        if mol.GetAtomWithIdx(b).GetHybridization() != SP2:
            continue
        if mol.GetAtomWithIdx(c).GetHybridization() != SP2:
            continue
        key = (min(b, c), max(b, c))
        if key in seen:
            continue
        # Find a ring containing both b and c, and pick the ring neighbors
        # (which are cis across the double bond).
        for ring in ri.AtomRings():
            if b not in ring or c not in ring:
                continue
            ring = list(ring)
            ib, ic = ring.index(b), ring.index(c)
            n = len(ring)
            # Ring neighbor of b on the far side from c
            nb_b = ring[(ib - 1) % n] if ring[(ib + 1) % n] == c else ring[(ib + 1) % n]
            # Ring neighbor of c on the far side from b
            nb_c = ring[(ic + 1) % n] if ring[(ic - 1) % n] == b else ring[(ic - 1) % n]
            if nb_b in atoms and nb_c in atoms:
                seen.add(key)
                quads.append((nb_b, b, c, nb_c))
                targets.append(0.0)
                break
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


def _build_embed_data(mol, bond_dihedral, atom_rings):
    """Build static arrays for C distance-geometry embedding."""
    n = mol.GetNumAtoms()
    atoms = set(range(n))

    bond_pairs = []
    bond_lengths = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond_pairs.append((i, j))
        bond_lengths.append(_get_bond_length(mol, i, j))

    angle_triples = []
    angle_values = []
    for b_idx in sorted(atoms):
        atom_b = mol.GetAtomWithIdx(b_idx)
        nbrs = [nb.GetIdx() for nb in atom_b.GetNeighbors() if nb.GetIdx() in atoms]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                angle_triples.append((a, b_idx, c))
                angle_values.append(_get_bond_angle(mol, a, b_idx, c, atom_rings))

    dihedral_quads = []
    dihedral_values = []
    seen = set()
    for (i, j), (mi, mj, angle) in bond_dihedral.items():
        key = (min(i, j), max(i, j))
        if key in seen:
            continue
        seen.add(key)
        dihedral_quads.append((mi, i, j, mj))
        dihedral_values.append(float(angle))

    bp = _to_int32(bond_pairs) if bond_pairs else np.empty((0, 2), dtype=np.int32)
    bl = _to_f64(bond_lengths) if bond_lengths else np.empty(0)
    at = _to_int32(angle_triples) if angle_triples else np.empty((0, 3), dtype=np.int32)
    av = _to_f64(angle_values) if angle_values else np.empty(0)
    dq = _to_int32(dihedral_quads) if dihedral_quads else np.empty((0, 4), dtype=np.int32)
    dv = _to_f64(dihedral_values) if dihedral_values else np.empty(0)

    lower = np.empty((n, n), dtype=np.float64)
    upper = np.empty((n, n), dtype=np.float64)
    angle_distances = np.empty(len(at), dtype=np.float64)
    _c_prepare_embed_bounds(
        n,
        _dptr(lower),
        _dptr(upper),
        _dptr(angle_distances),
        len(bp),
        _dptr(bp),
        _dptr(bl),
        len(at),
        _dptr(at),
        _dptr(av),
        len(dq),
        _dptr(dq),
        _dptr(dv),
    )
    return bp, bl, at, av, dq, dv, lower, upper, angle_distances


def _dg_embed(mol, bond_dihedral, atom_rings, seed=42, embed_data=None, coords=None):
    """Embed molecule using C distance-geometry with AMSR dihedrals.

    Builds distance bounds from bond lengths, bond angles, and AMSR-encoded
    dihedrals, then embeds via metric matrix eigendecomposition.
    Returns Nx3 coordinate array.
    """
    n = mol.GetNumAtoms()
    if embed_data is None:
        embed_data = _build_embed_data(mol, bond_dihedral, atom_rings)
    if coords is None:
        coords = np.zeros((n, 3), dtype=np.float64)
    if len(embed_data) == 9:
        bp, bl, at, _, _, _, lower, upper, angle_distances = embed_data
        _c_embed_bounds(
            n,
            _dptr(coords),
            _dptr(lower),
            _dptr(upper),
            len(bp),
            _dptr(bp),
            _dptr(bl),
            len(at),
            _dptr(at),
            _dptr(angle_distances),
            seed,
        )
    else:
        bp, bl, at, av, dq, dv = embed_data
        _c_embed(
            n,
            _dptr(coords),
            len(bp),
            _dptr(bp),
            _dptr(bl),
            len(at),
            _dptr(at),
            _dptr(av),
            len(dq),
            _dptr(dq),
            _dptr(dv),
            seed,
        )
    return coords


# ============================================================
# Optimizer
# ============================================================

_W_BOND = 5.0
_W_ANGLE = 2.0
_W_PLANAR = 5.0
_W_CHIRAL = 10.0
_W_DIHEDRAL = 0.1
_W_EZ = 0.3
_W_LINEAR = 20.0

DEFAULT_MAX_CONFS = 500


def _build_optimizer_data(mol, bond_dihedral, atom_rings=None):
    """Build static arrays for the C optimizer."""
    n = mol.GetNumAtoms()
    atoms = set(range(n))

    bonds, ideal_lengths = _collect_bonds(mol, atoms)
    angle_triples, ideal_angles = _collect_angles(mol, atoms, atom_rings)
    planar_groups = _collect_planar_atoms(mol, atoms)
    chiral_info, chiral_target_vols = _collect_chiral_atoms(mol, atoms)
    dih_quads, dih_targets = _collect_dihedral_restraints(mol, atoms, bond_dihedral)
    linear_triples = _collect_linear_atoms(mol, atoms)
    ez_quads, ez_targets = _collect_ez_constraints(mol, atoms)
    rp_quads, rp_targets = _collect_planarity_dihedrals(mol, atoms)

    if len(rp_quads):
        dih_quads = np.concatenate([dih_quads, rp_quads]) if len(dih_quads) else rp_quads
        dih_targets = np.concatenate([dih_targets, rp_targets]) if len(dih_targets) else rp_targets

    bp = _to_int32(bonds) if len(bonds) else np.empty((0, 2), dtype=np.int32)
    il = _to_f64(ideal_lengths) if len(ideal_lengths) else np.empty(0, dtype=np.float64)
    at = _to_int32(angle_triples) if len(angle_triples) else np.empty((0, 3), dtype=np.int32)
    ia = _to_f64(ideal_angles) if len(ideal_angles) else np.empty(0, dtype=np.float64)
    pg = _to_int32(planar_groups[:, :4]) if len(planar_groups) else np.empty((0, 4), dtype=np.int32)
    ci = _to_int32(chiral_info) if len(chiral_info) else np.empty((0, 5), dtype=np.int32)
    ctv = _to_f64(chiral_target_vols) if len(chiral_target_vols) else np.empty(0, dtype=np.float64)
    dq = _to_int32(dih_quads) if len(dih_quads) else np.empty((0, 4), dtype=np.int32)
    dt = _to_f64(dih_targets) if len(dih_targets) else np.empty(0, dtype=np.float64)
    eq = _to_int32(ez_quads) if len(ez_quads) else np.empty((0, 4), dtype=np.int32)
    et = _to_f64(ez_targets) if len(ez_targets) else np.empty(0, dtype=np.float64)
    lt = _to_int32(linear_triples) if len(linear_triples) else np.empty((0, 3), dtype=np.int32)

    return bp, il, at, ia, pg, ci, ctv, dq, dt, eq, et, lt


def _optimize(
    mol,
    bond_dihedral,
    coords,
    atom_rings=None,
    ftol=1e-3,
    gtol=1e-1,
    optimizer_data=None,
):
    """Optimize all atom positions to satisfy geometry constraints.

    Returns the optimizer cost.  Mutates coords in place.
    """
    n = mol.GetNumAtoms()
    if optimizer_data is None:
        optimizer_data = _build_optimizer_data(mol, bond_dihedral, atom_rings)
    bp, il, at, ia, pg, ci, ctv, dq, dt, eq, et, lt = optimizer_data

    x = coords.ravel()
    copied = x.dtype != np.float64 or not x.flags["C_CONTIGUOUS"]
    if copied:
        x = np.ascontiguousarray(x, dtype=np.float64)

    cost = _c_lbfgs(
        _dptr(x),
        len(x),
        n,
        _dptr(bp),
        _dptr(il),
        len(bp),
        _dptr(at),
        _dptr(ia),
        len(at),
        _dptr(pg),
        len(pg),
        _dptr(ci),
        _dptr(ctv),
        len(ci),
        _dptr(dq),
        _dptr(dt),
        len(dq),
        _dptr(eq),
        _dptr(et),
        len(eq),
        _dptr(lt),
        len(lt),
        _W_BOND,
        _W_ANGLE,
        _W_PLANAR,
        _W_CHIRAL,
        _W_DIHEDRAL,
        _W_EZ,
        _W_LINEAR,
        ftol,
        gtol,
        2000,
    )
    if copied:
        coords[:] = x.reshape(-1, 3)
    return cost


# ============================================================
# Public API
# ============================================================


def GetConformer(
    mol: Chem.Mol,
    dihedral: Optional[dict[tuple[int, int, int, int], int]] = None,
    ftol: float = 1e-3,
    gtol: float = 1e-1,
    max_confs: int = DEFAULT_MAX_CONFS,
) -> Chem.Mol:
    """Generate 3D conformer.

    1. Embed with distance geometry using AMSR-encoded dihedrals.
    2. Fix chirality, pseudo-E/Z, ring puckers, and acyclic dihedrals.
    3. Optimize all atom positions with a cost function enforcing ideal
       bond lengths, angles, planarity, chirality, and AMSR dihedrals.

    Stops early if any start reaches a very low internal cost.
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
    atom_rings = _build_ring_index(mol)
    subtree_cache = _build_subtree_cache(mol)
    chirality_ops = _build_chirality_ops(mol, subtree_cache)
    dihedral_ops = _build_dihedral_ops(mol, bond_dihedral, subtree_cache)
    pseudo_ez_ops = _build_pseudo_ez_ops(mol, subtree_cache)
    ring_topo = _build_ring_topology(mol, subtree_cache)
    embed_data = _build_embed_data(mol, bond_dihedral, atom_rings)
    optimizer_data = _build_optimizer_data(mol, bond_dihedral, atom_rings)

    for attempt in range(max_confs):
        _dg_embed(
            mol,
            bond_dihedral,
            atom_rings,
            seed=42 + attempt,
            embed_data=embed_data,
            coords=coords,
        )
        if chirality_ops:
            _fix_chirality(mol, coords, subtree_cache, chirality_ops)
        if dihedral_ops:
            _set_dihedrals(mol, bond_dihedral, coords, subtree_cache, dihedral_ops)
        if pseudo_ez_ops:
            _fix_pseudo_ez(mol, coords, subtree_cache, pseudo_ez_ops)
        if ring_topo[0]:
            _fix_ring_puckers(mol, bond_dihedral, coords, atom_rings, ring_topo)
        if dihedral_ops:
            _set_dihedrals(mol, bond_dihedral, coords, subtree_cache, dihedral_ops)
        oc = _optimize(
            mol,
            bond_dihedral,
            coords,
            atom_rings,
            ftol=ftol,
            gtol=gtol,
            optimizer_data=optimizer_data,
        )
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
