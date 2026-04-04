"""Python wrapper for the C cost_and_grad shared library."""

import ctypes
import os
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Load the C shared library
# ---------------------------------------------------------------------------

_lib = None
_c_func = None


def _load_lib():
    global _lib, _c_func
    d = os.path.dirname(__file__)
    name = "_cost_grad_c.dylib" if sys.platform == "darwin" else "_cost_grad_c.so"
    path = os.path.join(d, name)
    if not os.path.exists(path):
        raise RuntimeError(
            f"C extension {name} not found in {d}. Build it with 'make' in the project root."
        )
    try:
        _lib = ctypes.CDLL(path)
    except OSError as e:
        raise RuntimeError(f"Failed to load C extension {path}: {e}. Rebuild with 'make'.") from e
    fn = _lib.cost_and_grad
    fn.restype = ctypes.c_double
    fn.argtypes = [
        ctypes.c_void_p,  # x
        ctypes.c_void_p,  # grad
        ctypes.c_int,  # n_free
        ctypes.c_void_p,  # fixed
        ctypes.c_int,  # n_fixed
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int,  # bonds
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int,  # angles
        ctypes.c_void_p,
        ctypes.c_int,  # planar
        ctypes.c_void_p,
        ctypes.c_void_p,  # chiral (info + target_vols)
        ctypes.c_int,
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int,  # dih
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int,  # ez
        ctypes.c_void_p,
        ctypes.c_int,  # linear
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,  # weights (7: bond, angle, planar, chiral, dih, ez, linear)
    ]
    _c_func = fn


_load_lib()


def is_available():
    """True if the C library is loaded."""
    return _c_func is not None


def _to_int32(arr):
    return np.ascontiguousarray(arr, dtype=np.int32)


def _to_f64(arr):
    return np.ascontiguousarray(arr, dtype=np.float64)


def _dptr(arr):
    """Data pointer for a numpy array, or NULL."""
    return arr.ctypes.data_as(ctypes.c_void_p) if arr.size else ctypes.c_void_p(0)


class CostGradProblem:
    """Pre-cached C cost_and_grad arguments for repeated evaluation.

    All constraint arrays are converted to the right dtypes and their
    ctypes pointers cached once.  Only x and grad change per call.
    """

    __slots__ = (
        "_fn",
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
        "_w_bond",
        "_w_angle",
        "_w_planar",
        "_w_chiral",
        "_w_dih",
        "_w_ez",
        "_lt",
        "_n_linear",
        "_p_lt",
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
        self._fn = _c_func
        self._n_free = n_free
        self._grad = np.zeros(n_free * 3)

        # Pre-convert and cache pointers
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
        """Evaluate cost and gradient. Returns (cost, grad_flat)."""
        x = np.ascontiguousarray(x, dtype=np.float64)
        grad = self._grad
        cost = self._fn(
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
