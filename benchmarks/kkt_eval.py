#!/usr/bin/env python3
"""Evaluate the big-M KKT reformulation on 122 generated bilevel instances.

For every instance (integer leader, LP follower; see generators.py) this
script

  1. solves the KKT MILP (kkt_highs.py) with HiGHS at M = 1000,
  2. checks the returned (x, y) with one follower LP: (x, y) is
     bilevel feasible iff it satisfies all constraints and
     c^T y <= phi(x) + tol, where phi(x) is the follower optimum at x,
  3. computes the exact bilevel optimum by enumerating every integer x
     when there are at most MAX_ENUM leader vectors (for fixed x the
     optimistic problem is a single LP),
  4. records the LP relaxation of the KKT MILP and the high-point
     relaxation (follower optimality dropped),
  5. repeats the KKT solve for a range of M values and classifies each
     answer against the reference optimum and the follower-LP check.

Results go to kkt_eval_output/kkt_eval.json. Dependencies: numpy, highspy.
"""

import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import highspy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generators import (BilevelInstance, generate_hard_knapsack_interdiction,
                        generate_dense_bilevel,
                        generate_bilevel_with_integer_linking,
                        generate_stackelberg_game)
from kkt_highs import build_kkt_highs, KINTEGER

INF = 1e30
MAX_ENUM = 4096
MAIN_M = 1000.0
SWEEP_M = [1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 1000.0, 10000.0]
MAIN_TL = 60.0
SWEEP_TL = 30.0
REL_TOL = 1e-6


def instances():
    out = []
    for n in [8, 12, 16, 20, 25, 30, 40, 50]:
        for seed in [42, 123, 456]:
            for corr in ["weakly", "strongly", "uncorrelated"]:
                out.append(generate_hard_knapsack_interdiction(n, seed, corr))
    for nx, ny in [(5, 5), (8, 8), (10, 10), (10, 15), (15, 15), (20, 20), (20, 30)]:
        for seed in [42, 123]:
            for density in [0.5, 0.8]:
                out.append(generate_dense_bilevel(nx, ny, density, seed))
    for n in [5, 8, 10, 15, 20, 25]:
        for seed in [42, 123]:
            out.append(generate_bilevel_with_integer_linking(n, seed))
    for n in [5, 8, 10, 15, 20]:
        for seed in [42, 123]:
            out.append(generate_stackelberg_game(n, seed))
    return out


def _new_highs(time_limit=None):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    if time_limit is not None:
        h.setOptionValue("time_limit", float(time_limit))
    return h


def _add_rows(h, M, lo, hi, offset=0):
    for k in range(M.shape[0]):
        idx = [offset + j for j in range(M.shape[1]) if abs(M[k, j]) > 1e-12]
        val = [float(M[k, j - offset]) for j in idx]
        h.addRow(float(lo[k]), float(hi[k]), len(idx),
                 np.array(idx, dtype=np.int32), np.array(val, dtype=float))


def follower_value(inst: BilevelInstance, x: np.ndarray):
    """phi(x) = min { c^T y : A y <= b + B x, y >= 0 }, or None."""
    h = _new_highs()
    for j in range(inst.n_y):
        h.addVar(0.0, INF)
        h.changeColCost(j, float(inst.c[j]))
    rhs = inst.b + inst.B @ x
    _add_rows(h, inst.A, np.full(inst.m_lower, -INF), rhs)
    h.run()
    if h.getModelStatus() == highspy.HighsModelStatus.kOptimal:
        return h.getInfoValue("objective_function_value")[1]
    return None


def check_bilevel_feasible(inst, x, y, tol=REL_TOL):
    """One-LP check that (x, y) is feasible and follower-optimal."""
    feas_tol = 1e-6
    if np.any(y < -feas_tol):
        return False, None
    if np.any(inst.C @ x + inst.D @ y > inst.h + feas_tol * (1 + np.abs(inst.h))):
        return False, None
    rhs = inst.b + inst.B @ x
    if np.any(inst.A @ y > rhs + feas_tol * (1 + np.abs(rhs))):
        return False, None
    phi = follower_value(inst, x)
    if phi is None:
        return False, None
    return bool(inst.c @ y <= phi + tol * (1 + abs(phi))), phi


def optimistic_value_at(inst, x):
    """Leader objective of the optimistic bilevel problem with x fixed, or None."""
    phi = follower_value(inst, x)
    if phi is None:
        return None
    if np.any(inst.x_lb - 1e-9 > x) or np.any(x > inst.x_ub + 1e-9):
        return None
    h = _new_highs()
    for j in range(inst.n_y):
        h.addVar(0.0, INF)
        h.changeColCost(j, float(inst.e[j]))
    _add_rows(h, inst.D, np.full(inst.m_upper, -INF), inst.h - inst.C @ x)
    _add_rows(h, inst.A, np.full(inst.m_lower, -INF), inst.b + inst.B @ x)
    _add_rows(h, inst.c.reshape(1, -1), [-INF], [phi + 1e-9 * (1 + abs(phi))])
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        return None
    return float(inst.d @ x) + h.getInfoValue("objective_function_value")[1]


def exact_optimum(inst):
    """Enumerate integer x; return (optimum or None, number of x vectors) or (None, n) if too many."""
    ranges = [range(int(round(inst.x_lb[i])), int(round(inst.x_ub[i])) + 1)
              for i in range(inst.n_x)]
    count = int(np.prod([len(r) for r in ranges], dtype=float))
    if count > MAX_ENUM:
        return None, count, False
    best = None
    for xs in itertools.product(*ranges):
        v = optimistic_value_at(inst, np.array(xs, dtype=float))
        if v is not None and (best is None or v < best):
            best = v
    return best, count, True


def solve_kkt(inst, big_m, time_limit, relax=False):
    h = build_kkt_highs(inst.n_x, inst.n_y, inst.d, inst.e, inst.C, inst.D,
                        inst.h, inst.c, inst.A, inst.b, inst.B,
                        inst.x_lb, inst.x_ub, inst.x_binary, big_m,
                        relax_integrality=relax)
    h.setOptionValue("time_limit", float(time_limit))
    t0 = time.time()
    h.run()
    elapsed = time.time() - t0
    ms = h.getModelStatus()
    out = {"status": "other", "objective": None, "time_s": elapsed, "nodes": 0}
    if not relax:
        try:
            out["nodes"] = int(h.getInfoValue("mip_node_count")[1])
        except Exception:
            pass
    if ms == highspy.HighsModelStatus.kOptimal:
        out["status"] = "optimal"
    elif ms == highspy.HighsModelStatus.kInfeasible:
        out["status"] = "infeasible"
        return out, None, None
    elif ms == highspy.HighsModelStatus.kTimeLimit:
        out["status"] = "time_limit"
    if h.getInfoValue("primal_solution_status")[1] == 2:
        sol = h.getSolution()
        x = np.array(sol.col_value[:inst.n_x])
        y = np.array(sol.col_value[inst.n_x:inst.n_x + inst.n_y])
        out["objective"] = h.getInfoValue("objective_function_value")[1]
        return out, x, y
    return out, None, None


def solve_hpr(inst, time_limit):
    """High-point relaxation: all constraints, follower optimality dropped."""
    h = _new_highs(time_limit)
    for i in range(inst.n_x):
        h.addVar(float(inst.x_lb[i]), float(inst.x_ub[i]))
        h.changeColCost(i, float(inst.d[i]))
        if inst.x_binary:
            h.changeColIntegrality(i, KINTEGER)
    for j in range(inst.n_y):
        h.addVar(0.0, INF)
        h.changeColCost(inst.n_x + j, float(inst.e[j]))
    _add_rows(h, np.hstack([inst.C, inst.D]), np.full(inst.m_upper, -INF), inst.h)
    _add_rows(h, np.hstack([-inst.B, inst.A]), np.full(inst.m_lower, -INF), inst.b)
    h.run()
    if h.getModelStatus() == highspy.HighsModelStatus.kOptimal:
        return h.getInfoValue("objective_function_value")[1]
    return None


def same(a, b):
    return a is not None and b is not None and abs(a - b) <= 1e-6 * max(1.0, abs(b))


def main():
    out_dir = Path(__file__).parent / "kkt_eval_output"
    out_dir.mkdir(exist_ok=True)
    insts = instances()
    print(f"{len(insts)} instances")
    rows = []
    for inst in insts:
        t0 = time.time()
        exact, n_x_vectors, enumerated = exact_optimum(inst)
        t_enum = time.time() - t0
        lp, _, _ = solve_kkt(inst, MAIN_M, MAIN_TL, relax=True)
        hpr = solve_hpr(inst, MAIN_TL)
        main_r, x, y = solve_kkt(inst, MAIN_M, MAIN_TL)
        verified = None
        if x is not None:
            verified, _ = check_bilevel_feasible(inst, x, y)
        sweep = []
        for M in SWEEP_M:
            r, xs, ys = solve_kkt(inst, M, SWEEP_TL)
            r["big_m"] = M
            r["verified"] = None
            if xs is not None:
                r["verified"], _ = check_bilevel_feasible(inst, xs, ys)
            sweep.append(r)
        row = {
            "instance": inst.name, "category": inst.category,
            "n_x": inst.n_x, "n_y": inst.n_y, "m_lower": inst.m_lower,
            "n_x_vectors": n_x_vectors, "enumerated": enumerated,
            "exact_optimum": exact, "enum_time_s": t_enum if enumerated else None,
            "lp_relaxation": lp["objective"] if lp["status"] == "optimal" else None,
            "hpr": hpr,
            "kkt": {**main_r, "verified": verified},
            "sweep": sweep,
        }
        rows.append(row)
        print(f"{inst.name:<32} kkt={main_r['status']:<10} obj={main_r['objective']} "
              f"verified={verified} exact={exact} t={main_r['time_s']:.3f}s nodes={main_r['nodes']}",
              flush=True)

    summary = summarize(rows)
    print(json.dumps(summary, indent=2))
    with open(out_dir / "kkt_eval.json", "w") as f:
        json.dump({"date": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "highspy_version": getattr(highspy, "__version__", None),
                   "main_big_m": MAIN_M, "sweep_big_m": SWEEP_M,
                   "max_enum": MAX_ENUM, "summary": summary, "instances": rows},
                  f, indent=2, default=float)
    print(f"saved {out_dir / 'kkt_eval.json'}")


def reference(row):
    if row["enumerated"]:
        return row["exact_optimum"]
    big = [s for s in row["sweep"] if s["big_m"] == max(SWEEP_M)][0]
    return big["objective"] if big["status"] == "optimal" and big["verified"] else None


def summarize(rows):
    cats = {}
    for r in rows:
        c = cats.setdefault(r["category"], {"count": 0, "kkt_optimal": 0, "kkt_verified": 0,
                                            "enumerated": 0, "kkt_matches_exact": 0,
                                            "times": [], "nodes": [], "lp_gaps": [], "hpr_gaps": []})
        c["count"] += 1
        k = r["kkt"]
        if k["status"] == "optimal":
            c["kkt_optimal"] += 1
            c["times"].append(k["time_s"]); c["nodes"].append(k["nodes"])
            if k["verified"]:
                c["kkt_verified"] += 1
            if r["lp_relaxation"] is not None and abs(k["objective"]) > 1e-8:
                c["lp_gaps"].append(abs(k["objective"] - r["lp_relaxation"]) / abs(k["objective"]) * 100)
            if r["hpr"] is not None and abs(k["objective"]) > 1e-8:
                c["hpr_gaps"].append(abs(k["objective"] - r["hpr"]) / abs(k["objective"]) * 100)
        if r["enumerated"]:
            c["enumerated"] += 1
            if k["status"] == "optimal" and same(k["objective"], r["exact_optimum"]):
                c["kkt_matches_exact"] += 1
    per_cat = {}
    for name, c in cats.items():
        per_cat[name] = {
            "count": c["count"], "kkt_optimal": c["kkt_optimal"],
            "kkt_verified_bilevel_feasible": c["kkt_verified"],
            "enumerated": c["enumerated"], "kkt_matches_exact": c["kkt_matches_exact"],
            "mean_time_s": float(np.mean(c["times"])) if c["times"] else None,
            "max_time_s": float(np.max(c["times"])) if c["times"] else None,
            "mean_nodes": float(np.mean(c["nodes"])) if c["nodes"] else None,
            "max_nodes": int(np.max(c["nodes"])) if c["nodes"] else None,
            "mean_lp_gap_pct": float(np.mean(c["lp_gaps"])) if c["lp_gaps"] else None,
            "max_lp_gap_pct": float(np.max(c["lp_gaps"])) if c["lp_gaps"] else None,
            "zero_lp_gap": int(sum(g < 1e-6 for g in c["lp_gaps"])),
            "mean_hpr_gap_pct": float(np.mean(c["hpr_gaps"])) if c["hpr_gaps"] else None,
            "instances_hpr_differs": int(sum(g > 1e-6 for g in c["hpr_gaps"])),
        }

    sweep = {}
    for M in SWEEP_M:
        s = {"runs": 0, "correct": 0, "wrong_flagged_by_check": 0,
             "wrong_passed_check": 0, "reported_infeasible": 0, "no_answer": 0,
             "correct_failed_check": 0}
        for r in rows:
            ref = reference(r)
            if ref is None:
                continue
            run = [x for x in r["sweep"] if x["big_m"] == M][0]
            s["runs"] += 1
            if run["status"] == "infeasible":
                s["reported_infeasible"] += 1
            elif run["objective"] is None or run["status"] != "optimal":
                s["no_answer"] += 1
            elif same(run["objective"], ref):
                s["correct"] += 1
                if not run["verified"]:
                    s["correct_failed_check"] += 1
            elif run["verified"]:
                s["wrong_passed_check"] += 1
            else:
                s["wrong_flagged_by_check"] += 1
        sweep[str(M)] = s
    return {"categories": per_cat, "big_m_sweep": sweep,
            "instances": len(rows),
            "kkt_optimal": sum(c["kkt_optimal"] for c in per_cat.values()),
            "kkt_verified": sum(c["kkt_verified_bilevel_feasible"] for c in per_cat.values()),
            "enumerated": sum(c["enumerated"] for c in per_cat.values()),
            "kkt_matches_exact": sum(c["kkt_matches_exact"] for c in per_cat.values())}


if __name__ == "__main__":
    main()
