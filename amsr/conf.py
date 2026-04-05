"""Conformer generation from AMSR dihedrals.

Embed with distance geometry using AMSR-encoded dihedrals, then optimize
all atoms with a Cartesian-space cost function enforcing ideal bond
lengths, bond angles, planarity at SP2 centers, chirality at SP3
centers, E/Z constraints, and AMSR dihedral restraints.
"""

import ctypes
import os
import sys
from typing import Optional

import numpy as np
from rdkit import Chem

# ---------------------------------------------------------------------------
# C shared library (conf_util.c): cost_and_grad + embed
# ---------------------------------------------------------------------------

_LIB_NAME = "conf_util.dylib" if sys.platform == "darwin" else "conf_util.so"
_LIB_PATH = os.path.join(os.path.dirname(__file__), _LIB_NAME)
_lib = ctypes.CDLL(_LIB_PATH)

_cost_and_grad = _lib.cost_and_grad
_cost_and_grad.restype = ctypes.c_double
_cost_and_grad.argtypes = [
    ctypes.c_void_p,
    ctypes.c_void_p,
    ctypes.c_int,  # x, grad, n_free
    ctypes.c_void_p,
    ctypes.c_int,  # fixed, n_fixed
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


def _to_int32(arr, cols=None):
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
    return arr.ctypes.data_as(ctypes.c_void_p) if arr.size else ctypes.c_void_p(0)


class _CostGradProblem:
    """Pre-cached C cost_and_grad arguments for repeated evaluation."""

    __slots__ = (
        "_grad",
        "_n_free",
        "_fc",
        "_n_fixed",
        "_p_fc",
        "_bp",
        "_il",
        "_n_bonds",
        "_p_bp",
        "_p_il",
        "_at",
        "_ia",
        "_n_angles",
        "_p_at",
        "_p_ia",
        "_pg",
        "_n_planar",
        "_p_pg",
        "_ci",
        "_ctv",
        "_n_chiral",
        "_p_ci",
        "_p_ctv",
        "_dq",
        "_dt",
        "_n_dih",
        "_p_dq",
        "_p_dt",
        "_eq",
        "_et",
        "_n_ez",
        "_p_eq",
        "_p_et",
        "_lt",
        "_n_linear",
        "_p_lt",
        "_w_bond",
        "_w_angle",
        "_w_planar",
        "_w_chiral",
        "_w_dih",
        "_w_ez",
        "_w_linear",
    )

    def __init__(
        self,
        n_free,
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
        w_bond,
        w_angle,
        w_planar,
        w_chiral,
        w_dih,
        w_ez,
        linear_triples=None,
        w_linear=10.0,
    ):
        self._n_free = n_free
        self._grad = np.zeros(n_free * 3)
        self._fc = _to_f64(fixed_coords) if len(fixed_coords) else np.empty(0, dtype=np.float64)
        self._n_fixed = len(fixed_coords)
        self._p_fc = _dptr(self._fc)
        self._bp = _to_int32(bonds) if len(bonds) else np.empty((0, 2), dtype=np.int32)
        self._il = _to_f64(ideal_lengths) if len(ideal_lengths) else np.empty(0, dtype=np.float64)
        self._n_bonds = len(self._bp)
        self._p_bp = _dptr(self._bp)
        self._p_il = _dptr(self._il)
        self._at = (
            _to_int32(angle_triples) if len(angle_triples) else np.empty((0, 3), dtype=np.int32)
        )
        self._ia = _to_f64(ideal_angles) if len(ideal_angles) else np.empty(0, dtype=np.float64)
        self._n_angles = len(self._at)
        self._p_at = _dptr(self._at)
        self._p_ia = _dptr(self._ia)
        pg = planar_groups[:, :4] if len(planar_groups) else np.empty((0, 4), dtype=np.int32)
        self._pg = _to_int32(pg)
        self._n_planar = len(self._pg)
        self._p_pg = _dptr(self._pg)
        self._ci = _to_int32(chiral_info) if len(chiral_info) else np.empty((0, 5), dtype=np.int32)
        self._ctv = (
            _to_f64(chiral_target_vols)
            if len(chiral_target_vols)
            else np.empty(0, dtype=np.float64)
        )
        self._n_chiral = len(self._ci)
        self._p_ci = _dptr(self._ci)
        self._p_ctv = _dptr(self._ctv)
        self._dq = _to_int32(dih_quads) if len(dih_quads) else np.empty((0, 4), dtype=np.int32)
        self._dt = _to_f64(dih_targets) if len(dih_targets) else np.empty(0, dtype=np.float64)
        self._n_dih = len(self._dq)
        self._p_dq = _dptr(self._dq)
        self._p_dt = _dptr(self._dt)
        self._eq = _to_int32(ez_quads) if len(ez_quads) else np.empty((0, 4), dtype=np.int32)
        self._et = _to_f64(ez_targets) if len(ez_targets) else np.empty(0, dtype=np.float64)
        self._n_ez = len(self._eq)
        self._p_eq = _dptr(self._eq)
        self._p_et = _dptr(self._et)
        self._w_bond = w_bond
        self._w_angle = w_angle
        self._w_planar = w_planar
        self._w_chiral = w_chiral
        self._w_dih = w_dih
        self._w_ez = w_ez
        lt = linear_triples if linear_triples is not None and len(linear_triples) else None
        self._lt = _to_int32(lt) if lt is not None else np.empty((0, 3), dtype=np.int32)
        self._n_linear = len(self._lt)
        self._p_lt = _dptr(self._lt)
        self._w_linear = w_linear

    def __call__(self, x):
        x = np.ascontiguousarray(x, dtype=np.float64)
        grad = self._grad
        cost = _cost_and_grad(
            x.ctypes.data_as(ctypes.c_void_p),
            grad.ctypes.data_as(ctypes.c_void_p),
            self._n_free,
            self._p_fc,
            self._n_fixed,
            self._p_bp,
            self._p_il,
            self._n_bonds,
            self._p_at,
            self._p_ia,
            self._n_angles,
            self._p_pg,
            self._n_planar,
            self._p_ci,
            self._p_ctv,
            self._n_chiral,
            self._p_dq,
            self._p_dt,
            self._n_dih,
            self._p_eq,
            self._p_et,
            self._n_ez,
            self._p_lt,
            self._n_linear,
            self._w_bond,
            self._w_angle,
            self._w_planar,
            self._w_chiral,
            self._w_dih,
            self._w_ez,
            self._w_linear,
        )
        return cost, grad.copy()


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


def _fix_pseudo_ez(mol, coords):
    """Fix pseudo-E/Z double bonds after embedding.

    For double bonds with graph-equivalent substituents on at least one
    side, check if the E/Z geometry matches the tag and reflect if wrong.
    """
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    # Pseudo-E/Z double bonds
    for bond in mol.GetBonds():
        stereo = bond.GetStereo()
        if stereo not in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
            continue
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        stereo_atoms = list(bond.GetStereoAtoms())
        if len(stereo_atoms) < 2:
            continue
        si, sj = stereo_atoms[0], stereo_atoms[1]
        # Only fix pseudo-E/Z: substituents on at least one side share a rank
        ni = [nb.GetIdx() for nb in mol.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() != j]
        nj = [nb.GetIdx() for nb in mol.GetAtomWithIdx(j).GetNeighbors() if nb.GetIdx() != i]
        ni_ranks = [ranks[n] for n in ni]
        nj_ranks = [ranks[n] for n in nj]
        if len(set(ni_ranks)) == len(ni_ranks) and len(set(nj_ranks)) == len(nj_ranks):
            continue
        # Check current torsion
        actual = measure_torsion(coords[si], coords[i], coords[j], coords[sj])
        target = 0.0 if stereo == Chem.BondStereo.STEREOZ else 180.0
        diff = abs((actual - target + 180) % 360 - 180)
        if diff < 45:
            continue
        # Reflect the j-side substituents through the double bond plane
        bond_vec = coords[j] - coords[i]
        # Normal to the plane containing the double bond and si
        v_si = coords[si] - coords[i]
        normal = np.cross(bond_vec, v_si)
        nn = _norm3(normal)
        if nn < 1e-10:
            continue
        normal /= nn
        mid = 0.5 * (coords[i] + coords[j])
        for nb_idx in nj:
            sub = _subtree(mol, nb_idx, j)
            for k in sub:
                d = np.dot(coords[k] - mid, normal)
                coords[k] -= 2.0 * d * normal


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


def _fix_chirality(mol, coords):
    """Fix chiral centers whose signed volume has the wrong sign.

    For each wrong CW/CCW center, find the smallest neighbor subtree
    that is a true branch (< half the molecule) and reflect it through
    the plane of the other two neighbors.  Processes from interior
    to exterior and iterates until stable.
    """
    n = mol.GetNumAtoms()
    for _iteration in range(10):
        fixed_any = False
        for atom in mol.GetAtoms():
            chiral = atom.GetChiralTag()
            if chiral not in (CW, CCW):
                continue
            j = atom.GetIdx()
            nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
            if len(nbrs) < 3:
                continue
            rj = coords[j]
            v = [coords[ni] - rj for ni in nbrs[:3]]
            vol = np.dot(v[0], _cross3(v[1], v[2]))
            expected = -1 if chiral == CW else 1
            if (vol > 0) == (expected > 0):
                continue
            # Find smallest branch to reflect
            best = None
            for ni in nbrs[:3]:
                sub = _subtree(mol, ni, j)
                if len(sub) <= n // 2 and (best is None or len(sub) < len(best[1])):
                    best = (ni, sub)
            if best is None:
                continue
            reflect_nb, reflect_sub = best
            others = [ni for ni in nbrs[:3] if ni != reflect_nb]
            va, vb = coords[others[0]] - rj, coords[others[1]] - rj
            normal = _cross3(va, vb)
            nn = _norm3(normal)
            if nn < 1e-10:
                continue
            normal /= nn
            for k in reflect_sub:
                d = np.dot(coords[k] - rj, normal)
                coords[k] -= 2.0 * d * normal
            fixed_any = True
        if not fixed_any:
            break


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
        # Skip bonds involving SP (linear) atoms — torsion is undefined
        if mol.GetAtomWithIdx(i).GetHybridization() == SP:
            continue
        if mol.GetAtomWithIdx(j).GetHybridization() == SP:
            continue
        current = measure_torsion(coords[mi], coords[i], coords[j], coords[mj])
        delta = (target - current + 180) % 360 - 180
        if abs(delta) < 1.0:
            continue
        tree_j = _subtree(mol, j, i)
        tree_i = _subtree(mol, i, j)
        # Prefer rotating the smaller subtree to minimize cascading
        # disruption, but avoid rotating the side containing mi (the
        # reference atom), which would change the dihedral by 2*delta.
        if len(tree_j) <= len(tree_i) or mi in tree_i:
            _rotate_subtree(coords, tree_j - {j}, coords[i], coords[j] - coords[i], delta)
        else:
            _rotate_subtree(coords, tree_i - {i}, coords[j], coords[i] - coords[j], -delta)


def _fix_ring_puckers(mol, bond_dihedral, coords):
    """Fix ring conformations.

    Two strategies:
    1. For isolated rings with encoded dihedral targets, rebuild from
       internal coordinates (z-matrix style) and Kabsch-align.
    2. For any ring (including fused) with a chiral atom whose volume
       has the wrong sign, reflect the ring through its mean plane.
       Iterate until stable.
    """
    all_rings = [list(r) for r in mol.GetRingInfo().AtomRings()]
    all_ring_atoms = {}
    for idx, ring in enumerate(all_rings):
        for a in ring:
            all_ring_atoms.setdefault(a, set()).add(idx)
    shared_atoms = {a for a, rings in all_ring_atoms.items() if len(rings) > 1}

    # Strategy 2: flip ring puckers where chiral atoms have wrong sign.
    for _iteration in range(10):
        flipped_any = False
        for ring in all_rings:
            ring_set = set(ring)
            wrong = 0
            for a in ring:
                atom = mol.GetAtomWithIdx(a)
                chiral = atom.GetChiralTag()
                if chiral not in (CW, CCW):
                    continue
                nbrs = [nb.GetIdx() for nb in atom.GetNeighbors()]
                if len(nbrs) < 3:
                    continue
                rj = coords[a]
                v = [coords[ni] - rj for ni in nbrs[:3]]
                vol = np.dot(v[0], _cross3(v[1], v[2]))
                expected = -1 if chiral == CW else 1
                if (vol > 0) != (expected > 0):
                    wrong += 1
            if wrong == 0:
                continue
            # Reflect ring atoms through mean plane
            ring_coords = coords[ring]
            center = ring_coords.mean(axis=0)
            centered = ring_coords - center
            _, _, Vt = np.linalg.svd(centered)
            normal = Vt[2]
            for a in ring:
                d = np.dot(coords[a] - center, normal)
                coords[a] -= 2.0 * d * normal
                for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                    ni = nb.GetIdx()
                    if ni in ring_set:
                        continue
                    for k in _subtree(mol, ni, a):
                        d = np.dot(coords[k] - center, normal)
                        coords[k] -= 2.0 * d * normal
            flipped_any = True
        if not flipped_any:
            break

    # Strategy 1: rebuild isolated rings from internal coordinates.
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


def _dg_embed(mol, bond_dihedral, seed=42):
    """Embed molecule using C distance-geometry with AMSR dihedrals.

    Builds distance bounds from bond lengths, bond angles, and AMSR-encoded
    dihedrals, then embeds via metric matrix eigendecomposition.
    Returns Nx3 coordinate array.
    """
    n = mol.GetNumAtoms()
    atoms = set(range(n))

    # Collect bonds
    bond_pairs = []
    bond_lengths = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bond_pairs.append((i, j))
        bond_lengths.append(_get_bond_length(mol, i, j))

    # Collect angles
    angle_triples = []
    angle_values = []
    for b_idx in sorted(atoms):
        atom_b = mol.GetAtomWithIdx(b_idx)
        nbrs = [nb.GetIdx() for nb in atom_b.GetNeighbors() if nb.GetIdx() in atoms]
        for ia in range(len(nbrs)):
            for ic in range(ia + 1, len(nbrs)):
                a, c = nbrs[ia], nbrs[ic]
                angle_triples.append((a, b_idx, c))
                angle_values.append(_get_bond_angle(mol, a, b_idx, c))

    # Collect AMSR dihedrals
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
    coords = np.zeros((n, 3), dtype=np.float64)
    _c_embed(
        n,
        coords.ctypes.data,
        len(bl),
        bp.ctypes.data,
        bl.ctypes.data,
        len(av),
        at.ctypes.data,
        av.ctypes.data,
        len(dv),
        dq.ctypes.data,
        dv.ctypes.data,
        seed,
    )
    return coords


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

    _objective = _CostGradProblem(
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
    coords[:] = result.x.reshape(-1, 3)
    return result.fun


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

    1. Embed with distance geometry using AMSR-encoded dihedrals.
    2. Fix chirality, pseudo-E/Z, ring puckers, and acyclic dihedrals.
    3. Optimize all atom positions with a cost function enforcing ideal
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
        coords[:] = _dg_embed(mol, bond_dihedral, seed=42 + attempt)
        _fix_chirality(mol, coords)
        _set_dihedrals(mol, bond_dihedral, coords)
        _fix_pseudo_ez(mol, coords)
        _fix_ring_puckers(mol, bond_dihedral, coords)
        _set_dihedrals(mol, bond_dihedral, coords)
        oc = _optimize(mol, bond_dihedral, coords, ftol=ftol, gtol=gtol)
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
