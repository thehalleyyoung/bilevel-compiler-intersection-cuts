"""Big-M KKT reformulation of a bilevel program with an LP follower, built in HiGHS.

Bilevel program (optimistic):

    min_{x,y}  d^T x + e^T y
    s.t.       C x + D y <= h
               x_lb <= x <= x_ub   (x integer if x_integer)
               y in argmin_{y'} { c^T y' : A y' <= b + B x, 0 <= y' <= y_ub }

The follower LP is replaced by its KKT conditions. With A', b', B' the
follower rows extended by the finite upper bounds y <= y_ub:

    A' y - B' x <= b'                       (primal feasibility)
    A'^T lam + c >= 0,  lam >= 0            (dual feasibility, reduced cost r = A'^T lam + c)
    lam_k <= M z_k,  (b' + B' x - A' y)_k <= M (1 - z_k)
    y_j   <= M t_j,  r_j                  <= M (1 - t_j)

z and t are binary. The formulation is exact whenever M bounds every
multiplier, every slack, every y_j and every reduced cost at some bilevel
optimum; a too-small M can cut off the optimum or return a point that is
not follower-optimal.

Column layout: x [0, n_x), y [n_x, n_x+n_y), lam (m'), z (m'), t (n_y).
"""

from typing import Optional

import numpy as np
import highspy

KINTEGER = highspy.HighsVarType.kInteger
INF = 1e30


def build_kkt_highs(n_x: int, n_y: int,
                    d: np.ndarray, e: np.ndarray,
                    C: np.ndarray, D: np.ndarray, h_vec: np.ndarray,
                    c: np.ndarray, A: np.ndarray, b: np.ndarray, B: np.ndarray,
                    x_lb: np.ndarray, x_ub: np.ndarray, x_integer: bool,
                    big_m: float,
                    y_ub: Optional[np.ndarray] = None,
                    relax_integrality: bool = False) -> highspy.Highs:
    """Return a HiGHS model of the big-M KKT reformulation described above."""
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)

    # Extend follower rows with finite upper bounds on y.
    rows_A, rows_B, rows_b = [A], [B], [b]
    if y_ub is not None:
        for j in range(n_y):
            if y_ub[j] < INF / 10:
                a = np.zeros((1, n_y)); a[0, j] = 1.0
                rows_A.append(a)
                rows_B.append(np.zeros((1, n_x)))
                rows_b.append(np.array([float(y_ub[j])]))
    Ap = np.vstack(rows_A)
    Bp = np.vstack(rows_B)
    bp = np.concatenate(rows_b)
    m = Ap.shape[0]

    off_y = n_x
    off_lam = n_x + n_y
    off_z = off_lam + m
    off_t = off_z + m

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)

    for i in range(n_x):
        h.addVar(float(x_lb[i]), float(x_ub[i]))
        h.changeColCost(i, float(d[i]))
        if x_integer and not relax_integrality:
            h.changeColIntegrality(i, KINTEGER)
    for j in range(n_y):
        h.addVar(0.0, INF)
        h.changeColCost(off_y + j, float(e[j]))
    for k in range(m):
        h.addVar(0.0, INF)
    for k in range(m):
        h.addVar(0.0, 1.0)
        if not relax_integrality:
            h.changeColIntegrality(off_z + k, KINTEGER)
    for j in range(n_y):
        h.addVar(0.0, 1.0)
        if not relax_integrality:
            h.changeColIntegrality(off_t + j, KINTEGER)

    def add_row(lo, hi, idx, val):
        h.addRow(lo, hi, len(idx), np.array(idx, dtype=np.int32),
                 np.array(val, dtype=float))

    # Upper level: C x + D y <= h
    for i in range(C.shape[0]):
        idx, val = [], []
        for j in range(n_x):
            if abs(C[i, j]) > 1e-12: idx.append(j); val.append(C[i, j])
        for j in range(n_y):
            if abs(D[i, j]) > 1e-12: idx.append(off_y + j); val.append(D[i, j])
        if idx:
            add_row(-INF, float(h_vec[i]), idx, val)

    # Primal feasibility: A' y - B' x <= b'
    for k in range(m):
        idx, val = [], []
        for j in range(n_y):
            if abs(Ap[k, j]) > 1e-12: idx.append(off_y + j); val.append(Ap[k, j])
        for j in range(n_x):
            if abs(Bp[k, j]) > 1e-12: idx.append(j); val.append(-Bp[k, j])
        add_row(-INF, float(bp[k]), idx, val)

    # Dual feasibility: A'^T lam >= -c
    for j in range(n_y):
        idx, val = [], []
        for k in range(m):
            if abs(Ap[k, j]) > 1e-12: idx.append(off_lam + k); val.append(Ap[k, j])
        add_row(float(-c[j]), INF, idx, val)

    # Complementarity between lam and the primal slack.
    for k in range(m):
        add_row(-INF, 0.0, [off_lam + k, off_z + k], [1.0, -big_m])
        idx, val = [], []
        for j in range(n_y):
            if abs(Ap[k, j]) > 1e-12: idx.append(off_y + j); val.append(-Ap[k, j])
        for j in range(n_x):
            if abs(Bp[k, j]) > 1e-12: idx.append(j); val.append(Bp[k, j])
        idx.append(off_z + k); val.append(big_m)
        add_row(-INF, big_m - float(bp[k]), idx, val)

    # Complementarity between y and its reduced cost r = A'^T lam + c.
    for j in range(n_y):
        add_row(-INF, 0.0, [off_y + j, off_t + j], [1.0, -big_m])
        idx, val = [], []
        for k in range(m):
            if abs(Ap[k, j]) > 1e-12: idx.append(off_lam + k); val.append(Ap[k, j])
        idx.append(off_t + j); val.append(big_m)
        add_row(-INF, big_m - float(c[j]), idx, val)

    return h


def num_kkt_columns(n_x: int, n_y: int, m_lower: int, n_y_bounds: int = 0) -> int:
    m = m_lower + n_y_bounds
    return n_x + n_y + 2 * m + n_y
