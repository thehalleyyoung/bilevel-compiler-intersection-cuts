#!/usr/bin/env python3
"""BOBILib evaluation: KKT on the LP-relaxed follower versus an exact
leader-enumeration loop, on the first 20 entries of bobilib/instance_list.json.

Methods (all solved with SCIP through PySCIPOpt):

  KKT          KKT conditions of the follower's LP relaxation, complementarity
               encoded with indicator constraints. Exact when the follower is
               continuous; for an integer follower the answer is then checked
               by solving the follower MIP at the returned leader decision.

  Enumeration  High-point relaxation (HPR) as a lower bound. At each HPR
               optimum x^, the follower MIP gives phi(x^) and one more MIP gives
               the optimistic leader value at x^ (minimise the leader
               objective over follower-optimal responses). A no-good
               constraint then removes x^ from the HPR. The loop stops when the
               HPR bound reaches the incumbent or the HPR becomes infeasible,
               which proves optimality. It requires a pure-integer leader
               with finite bounds; all 20 instances satisfy this.

Every objective is compared with the optimum BOBILib reports for the
instance. Results go to bobilib_results/bobilib_eval.json.
"""

import gzip, json, os, sys, time, traceback
import numpy as np

# ---------------------------------------------------------------------------
# Parsers (same as v3)
# ---------------------------------------------------------------------------

def parse_mps_gz(path):
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt', errors='replace') as f:
        lines = f.readlines()
    rows, row_sense, obj_name = [], {}, ''
    cols, seen_cols = [], set()
    A, rhs, lo, hi = {}, {}, {}, {}
    integers = set()
    section, in_integer = None, False
    for raw in lines:
        line = raw.rstrip('\n')
        if not line or line[0] == '*': continue
        if line[0] != ' ':
            tok = line.split()[0]
            if tok in ('NAME','ROWS','COLUMNS','RHS','RANGES','BOUNDS','ENDATA'):
                section = tok
            continue
        parts = line.split()
        if section == 'ROWS' and len(parts) >= 2:
            rows.append(parts[1]); row_sense[parts[1]] = parts[0]
            if parts[0] == 'N': obj_name = parts[1]
        elif section == 'COLUMNS':
            if "'MARKER'" in line:
                in_integer = "'INTORG'" in line; continue
            if len(parts) < 3: continue
            col = parts[0]
            if col not in seen_cols:
                seen_cols.add(col); cols.append(col); lo[col] = 0.0; hi[col] = 1e30
            if in_integer: integers.add(col)
            i = 1
            while i+1 < len(parts):
                A[(parts[i], col)] = float(parts[i+1]); i += 2
        elif section == 'RHS' and len(parts) >= 3:
            i = 1
            while i+1 < len(parts):
                rhs[parts[i]] = float(parts[i+1]); i += 2
        elif section == 'BOUNDS' and len(parts) >= 3:
            bt, _, col = parts[0], parts[1], parts[2]
            val = float(parts[3]) if len(parts) > 3 else 0.0
            if bt == 'UP': hi[col] = val
            elif bt == 'LO': lo[col] = val
            elif bt == 'FX': lo[col] = hi[col] = val
            elif bt == 'FR': lo[col] = -1e30; hi[col] = 1e30
            elif bt == 'MI': lo[col] = -1e30
            elif bt == 'PL': hi[col] = 1e30
            elif bt == 'BV': lo[col] = 0.0; hi[col] = 1.0; integers.add(col)
    return {'obj': obj_name, 'rows': rows, 'sense': row_sense,
            'cols': cols, 'A': A, 'rhs': rhs, 'lo': lo, 'hi': hi, 'int': integers}

def parse_aux(path):
    with open(path) as f:
        lines = f.read().strip().split('\n')
    lower_vars, lower_obj, lower_constrs = [], {}, []
    sec = None
    for line in lines:
        l = line.strip()
        if not l: continue
        if l in ('@VARSBEGIN','@VARSEND','@CONSTRSBEGIN','@CONSTRSEND',
                 '@NUMVARS','@NUMCONSTRS','@NAME','@MPS'):
            sec = l; continue
        if sec == '@VARSBEGIN':
            p = l.split(); lower_vars.append(p[0])
            lower_obj[p[0]] = float(p[1]) if len(p) > 1 else 0.0
        elif sec == '@CONSTRSBEGIN':
            lower_constrs.append(l.split()[0])
        elif sec in ('@NUMVARS','@NUMCONSTRS','@NAME','@MPS'):
            sec = None
    return lower_vars, lower_obj, lower_constrs

def _extract_scip(m, t0):
    elapsed = time.time() - t0
    r = {'time': elapsed, 'nodes': 0}
    try: r['nodes'] = int(m.getNNodes())
    except: pass
    st = m.getStatus()
    if st == 'optimal':
        r['status'] = 'optimal'; r['objective'] = m.getObjVal()
    elif st == 'infeasible':
        r['status'] = 'infeasible'; r['objective'] = None
    elif st in ('timelimit', 'gaplimit'):
        r['status'] = 'time_limit'
        try: r['objective'] = m.getObjVal()
        except: r['objective'] = None
    else:
        r['status'] = f'other:{st}'
        try: r['objective'] = m.getObjVal()
        except: r['objective'] = None
    return r

# ---------------------------------------------------------------------------
# Check feasibility of a full solution (x_val dict) against all constraints
# ---------------------------------------------------------------------------

def check_feasibility(mps, x_val, tol=1e-5):
    """Check if x_val satisfies all constraints and bounds."""
    for v in mps['cols']:
        val = x_val.get(v, 0.0)
        lb = mps['lo'].get(v, 0.0)
        ub = mps['hi'].get(v, 1e30)
        if val < lb - tol or val > ub + tol:
            return False
    for row in mps['rows']:
        s = mps['sense'].get(row)
        if s == 'N': continue
        lhs = sum(mps['A'].get((row, v), 0.0) * x_val.get(v, 0.0) for v in mps['cols'])
        rval = mps['rhs'].get(row, 0.0)
        if s == 'L' and lhs > rval + tol: return False
        if s == 'G' and lhs < rval - tol: return False
        if s == 'E' and abs(lhs - rval) > tol: return False
    return True

# ---------------------------------------------------------------------------
# Solve integer follower for fixed x
# ---------------------------------------------------------------------------

def solve_follower_mip(mps, lower_vars, lower_obj, lower_constrs, x_val, time_limit=10,
                       return_status=False):
    from pyscipopt import Model, quicksum
    lower_constr_set = set(lower_constrs)
    all_vars = mps['cols']
    lower_rows = [r for r in mps['rows'] if mps['sense'].get(r) != 'N' and r in lower_constr_set]

    m = Model(); m.hideOutput()
    m.setParam("limits/time", time_limit); m.setParam("limits/gap", 1e-8)

    yv = {}
    for v in lower_vars:
        lb = max(mps['lo'].get(v, 0.0), -1e20)
        ub = min(mps['hi'].get(v, 1e30), 1e20)
        yv[v] = m.addVar(name=v, lb=lb if lb > -1e20 else None,
                          ub=ub if ub < 1e20 else None,
                          vtype="I" if v in mps['int'] else "C")

    c_lower = {v: lower_obj.get(v, 0.0) for v in lower_vars}
    m.setObjective(quicksum(c_lower.get(v, 0) * yv[v] for v in lower_vars
                             if abs(c_lower.get(v, 0)) > 1e-15), "minimize")

    for row in lower_rows:
        s = mps['sense'][row]; rval = mps['rhs'].get(row, 0.0)
        rhs_adj = rval
        for v in all_vars:
            if v in yv: continue
            a = mps['A'].get((row, v), 0.0)
            if abs(a) > 1e-15: rhs_adj -= a * x_val.get(v, 0.0)
        lhs = quicksum(mps['A'].get((row, v), 0) * yv[v] for v in lower_vars
                        if abs(mps['A'].get((row, v), 0)) > 1e-15)
        if s == 'L': m.addCons(lhs <= rhs_adj)
        elif s == 'G': m.addCons(lhs >= rhs_adj)
        elif s == 'E': m.addCons(lhs == rhs_adj)

    m.optimize()
    st = m.getStatus()
    if st in ('optimal', 'gaplimit'):
        out = (m.getObjVal(), {v: m.getVal(yv[v]) for v in lower_vars})
    else:
        out = (None, None)
    return out + (st,) if return_status else out

# ---------------------------------------------------------------------------
# Method 1: Indicator KKT (baseline)
# ---------------------------------------------------------------------------

def solve_indicator(mps, lower_vars, lower_obj, lower_constrs, time_limit, return_vals=False):
    from pyscipopt import Model, quicksum
    lower_var_set = set(lower_vars)
    lower_constr_set = set(lower_constrs)
    all_vars = mps['cols']
    constraint_rows = [r for r in mps['rows'] if mps['sense'].get(r) != 'N']
    lower_rows = [r for r in constraint_rows if r in lower_constr_set]
    c_upper = {v: mps['A'].get((mps['obj'], v), 0.0) for v in all_vars}
    c_lower = {v: lower_obj.get(v, 0.0) for v in lower_vars}

    m = Model(); m.hideOutput()
    m.setParam("limits/time", time_limit); m.setParam("limits/gap", 1e-6)

    x = {}
    for v in all_vars:
        lb = max(mps['lo'].get(v, 0.0), -1e20)
        ub = min(mps['hi'].get(v, 1e30), 1e20)
        x[v] = m.addVar(name=v, lb=lb if lb > -1e20 else None,
                         ub=ub if ub < 1e20 else None,
                         vtype="I" if v in mps['int'] else "C")

    m.setObjective(quicksum(c_upper.get(v, 0) * x[v] for v in all_vars
                             if abs(c_upper.get(v, 0)) > 1e-15), "minimize")
    for row in constraint_rows:
        s = mps['sense'][row]; rval = mps['rhs'].get(row, 0.0)
        lhs = quicksum(mps['A'].get((row, v), 0) * x[v] for v in all_vars
                        if abs(mps['A'].get((row, v), 0)) > 1e-15)
        if s == 'L': m.addCons(lhs <= rval, name=row)
        elif s == 'G': m.addCons(lhs >= rval, name=row)
        elif s == 'E': m.addCons(lhs == rval, name=row)

    if not lower_rows or not lower_vars:
        t0 = time.time(); m.optimize()
        r = _extract_scip(m, t0)
        if return_vals and r.get('objective') is not None:
            r['x_val'] = {v: m.getVal(x[v]) for v in all_vars}
        return r

    std_rows = []
    for row in lower_rows:
        s = mps['sense'][row]; sign = -1.0 if s == 'G' else 1.0
        std_rows.append((row, sign, s))

    lam = {}
    for k, (row, sign, s) in enumerate(std_rows):
        lam[k] = m.addVar(f"lam_{k}", lb=None if s == 'E' else 0, ub=None, vtype="C")

    mu_lo, mu_hi = {}, {}
    for v in lower_vars:
        if mps['lo'].get(v, 0.0) > -1e20:
            mu_lo[v] = m.addVar(f"mu_lo_{v}", lb=0, vtype="C")
        if mps['hi'].get(v, 1e30) < 1e20:
            mu_hi[v] = m.addVar(f"mu_hi_{v}", lb=0, vtype="C")

    for v in lower_vars:
        grad = c_lower.get(v, 0.0); terms = []
        for k, (row, sign, _) in enumerate(std_rows):
            a = mps['A'].get((row, v), 0.0)
            if abs(a) < 1e-15: continue
            terms.append(lam[k] * (sign * a))
        stat_expr = grad + quicksum(terms) if terms else grad
        if v in mu_lo: stat_expr -= mu_lo[v]
        if v in mu_hi: stat_expr += mu_hi[v]
        m.addCons(stat_expr == 0, name=f"stat_{v}")

    for k, (row, sign, s) in enumerate(std_rows):
        if s == 'E': continue
        rval = mps['rhs'].get(row, 0.0)
        z = m.addVar(f"z_{k}", vtype="B")
        m.addConsIndicator(lam[k] <= 0, z, activeone=False)
        slack_expr = sign * (rval - quicksum(
            mps['A'].get((row, v), 0) * x[v] for v in all_vars
            if abs(mps['A'].get((row, v), 0)) > 1e-15))
        m.addConsIndicator(slack_expr <= 0, z, activeone=True)

    for v in lower_vars:
        lb_v = mps['lo'].get(v, 0.0); ub_v = mps['hi'].get(v, 1e30)
        if v in mu_lo:
            zl = m.addVar(f"zl_{v}", vtype="B")
            m.addConsIndicator(mu_lo[v] <= 0, zl, activeone=False)
            m.addConsIndicator(x[v] - lb_v <= 0, zl, activeone=True)
        if v in mu_hi and ub_v < 1e20:
            zu = m.addVar(f"zu_{v}", vtype="B")
            m.addConsIndicator(mu_hi[v] <= 0, zu, activeone=False)
            m.addConsIndicator(ub_v - x[v] <= 0, zu, activeone=True)

    t0 = time.time(); m.optimize()
    r = _extract_scip(m, t0)
    if return_vals and r.get('objective') is not None:
        r['x_val'] = {v: m.getVal(x[v]) for v in all_vars}
    return r


# ---------------------------------------------------------------------------
# Method 2: exact leader enumeration with HPR bound
# ---------------------------------------------------------------------------

def optimistic_value_at(mps, lower_vars, lower_obj, lower_constrs, x_lead, time_limit):
    """Optimistic bilevel objective with the leader fixed at x_lead.

    Returns (value or None, phi or None, exact). value is None when the
    follower problem is infeasible at x_lead or no follower-optimal response
    satisfies the leader's constraints; exact is False when a sub-MIP stopped
    before proving optimality or infeasibility.
    """
    from pyscipopt import Model, quicksum
    phi, _, fst = solve_follower_mip(mps, lower_vars, lower_obj, lower_constrs,
                                     x_lead, time_limit=max(1.0, time_limit / 2),
                                     return_status=True)
    if phi is None:
        return None, None, fst == 'infeasible'

    lower_set = set(lower_vars)
    all_vars = mps['cols']
    m = Model(); m.hideOutput()
    m.setParam("limits/time", max(1.0, time_limit / 2)); m.setParam("limits/gap", 1e-9)
    y = {}
    for v in lower_vars:
        lb = max(mps['lo'].get(v, 0.0), -1e20)
        ub = min(mps['hi'].get(v, 1e30), 1e20)
        y[v] = m.addVar(name=v, lb=lb if lb > -1e20 else None,
                        ub=ub if ub < 1e20 else None,
                        vtype="I" if v in mps['int'] else "C")
    const = 0.0
    for v in all_vars:
        if v not in lower_set:
            const += mps['A'].get((mps['obj'], v), 0.0) * x_lead[v]
    m.setObjective(quicksum(mps['A'].get((mps['obj'], v), 0.0) * y[v] for v in lower_vars
                            if abs(mps['A'].get((mps['obj'], v), 0.0)) > 1e-15), "minimize")
    for row in mps['rows']:
        s = mps['sense'].get(row)
        if s == 'N':
            continue
        rhs = mps['rhs'].get(row, 0.0)
        for v in all_vars:
            if v not in lower_set:
                rhs -= mps['A'].get((row, v), 0.0) * x_lead[v]
        terms = [mps['A'].get((row, v), 0.0) * y[v] for v in lower_vars
                 if abs(mps['A'].get((row, v), 0.0)) > 1e-15]
        lhs = quicksum(terms) if terms else 0.0
        if not terms:
            ok = (s == 'L' and 0.0 <= rhs + 1e-9) or (s == 'G' and 0.0 >= rhs - 1e-9) \
                or (s == 'E' and abs(rhs) <= 1e-9)
            if not ok:
                return None, phi, True
            continue
        if s == 'L': m.addCons(lhs <= rhs)
        elif s == 'G': m.addCons(lhs >= rhs)
        elif s == 'E': m.addCons(lhs == rhs)
    m.addCons(quicksum(lower_obj.get(v, 0.0) * y[v] for v in lower_vars
                       if abs(lower_obj.get(v, 0.0)) > 1e-15)
              <= phi + 1e-6 * max(1.0, abs(phi)))
    m.optimize()
    st = m.getStatus()
    if st in ('optimal', 'gaplimit'):
        return const + m.getObjVal(), phi, True
    return None, phi, st == 'infeasible'


def solve_enumeration(mps, lower_vars, lower_obj, lower_constrs, time_limit, max_iters=500):
    from pyscipopt import Model, quicksum
    t0 = time.time()
    lower_set = set(lower_vars)
    all_vars = mps['cols']
    leader = [v for v in all_vars if v not in lower_set]
    if any(v not in mps['int'] for v in leader):
        return {'status': 'unsupported', 'objective': None, 'time': 0.0, 'nodes': 0}
    bounds = {}
    for v in leader:
        lo, hi = mps['lo'].get(v, 0.0), mps['hi'].get(v, 1e30)
        if lo < -1e19 or hi > 1e19:
            return {'status': 'unsupported', 'objective': None, 'time': 0.0, 'nodes': 0}
        bounds[v] = (int(round(lo)), int(round(hi)))

    m = Model(); m.hideOutput()
    m.setParam("limits/gap", 1e-9)
    x = {}
    for v in all_vars:
        lb = max(mps['lo'].get(v, 0.0), -1e20)
        ub = min(mps['hi'].get(v, 1e30), 1e20)
        x[v] = m.addVar(name=v, lb=lb if lb > -1e20 else None,
                        ub=ub if ub < 1e20 else None,
                        vtype="I" if v in mps['int'] else "C")
    m.setObjective(quicksum(mps['A'].get((mps['obj'], v), 0.0) * x[v] for v in all_vars
                            if abs(mps['A'].get((mps['obj'], v), 0.0)) > 1e-15), "minimize")
    for row in mps['rows']:
        s = mps['sense'].get(row)
        if s == 'N':
            continue
        rval = mps['rhs'].get(row, 0.0)
        lhs = quicksum(mps['A'].get((row, v), 0) * x[v] for v in all_vars
                       if abs(mps['A'].get((row, v), 0)) > 1e-15)
        if s == 'L': m.addCons(lhs <= rval)
        elif s == 'G': m.addCons(lhs >= rval)
        elif s == 'E': m.addCons(lhs == rval)

    best, best_x, proved, iters, lb, all_exact = None, None, False, 0, None, True
    while iters < max_iters:
        remaining = time_limit - (time.time() - t0)
        if remaining < 1.0:
            break
        m.setParam("limits/time", remaining)
        m.optimize()
        st = m.getStatus()
        if st == 'infeasible':
            proved = True
            break
        if st not in ('optimal', 'gaplimit'):
            break
        lb = m.getObjVal()
        if best is not None and lb >= best - 1e-6 * max(1.0, abs(best)):
            proved = True
            break
        iters += 1
        x_hat = {v: float(int(round(m.getVal(x[v])))) for v in leader}
        val, _, exact = optimistic_value_at(mps, lower_vars, lower_obj, lower_constrs, x_hat,
                                            max(1.0, time_limit - (time.time() - t0)))
        all_exact = all_exact and exact
        if val is not None and (best is None or val < best):
            best, best_x = val, x_hat
        # No-good: exclude the integer leader vector x_hat.
        m.freeTransform()
        picks = []
        for v in leader:
            lo, hi = bounds[v]
            xv = int(x_hat[v])
            if xv > lo:
                d = m.addVar(vtype="B")
                m.addCons(x[v] <= xv - 1 + (hi - xv + 1) * (1 - d))
                picks.append(d)
            if xv < hi:
                d = m.addVar(vtype="B")
                m.addCons(x[v] >= xv + 1 - (xv + 1 - lo) * (1 - d))
                picks.append(d)
        if not picks:
            proved = True
            break
        m.addCons(quicksum(picks) >= 1)

    proved = proved and all_exact
    return {'status': 'optimal' if (proved and best is not None) else
                      ('infeasible' if proved else ('feasible' if best is not None else 'no_solution')),
            'objective': best, 'proved_optimal': proved, 'lower_bound': lb,
            'time': time.time() - t0, 'nodes': 0, 'iterations': iters}


def solve_kkt_checked(mps, lower_vars, lower_obj, lower_constrs, time_limit):
    """KKT of the LP-relaxed follower, then check the answer with the follower MIP."""
    r = solve_indicator(mps, lower_vars, lower_obj, lower_constrs, time_limit, return_vals=True)
    x_val = r.pop('x_val', None)
    r['bilevel_feasible'] = None
    if x_val is not None:
        phi, _ = solve_follower_mip(mps, lower_vars, lower_obj, lower_constrs, x_val, time_limit=10)
        fy = sum(lower_obj.get(v, 0.0) * x_val[v] for v in lower_vars)
        y_int = all(abs(x_val[v] - round(x_val[v])) <= 1e-6 for v in lower_vars if v in mps['int'])
        r['bilevel_feasible'] = bool(phi is not None and y_int
                                     and fy <= phi + 1e-6 * max(1.0, abs(phi))
                                     and check_feasibility(mps, x_val))
    return r


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    list_file = os.path.join(script_dir, 'bobilib', 'instance_list.json')
    if not os.path.exists(list_file):
        print("ERROR: instance_list.json not found"); sys.exit(1)

    with open(list_file) as f:
        instance_list = json.load(f)

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=20, help="number of instances from instance_list.json")
    ap.add_argument("--time-limit", type=float, default=60.0, help="seconds per method and instance")
    ap.add_argument("--out", default="bobilib_eval.json")
    args = ap.parse_args()
    os.chdir(script_dir)
    MAX = args.max
    instance_list = instance_list[:MAX]
    TL = args.time_limit

    methods = [
        ("KKT",         lambda mps, lv, lo, lc: solve_kkt_checked(mps, lv, lo, lc, TL)),
        ("Enumeration", lambda mps, lv, lo, lc: solve_enumeration(mps, lv, lo, lc, TL)),
    ]

    print("=" * 100)
    print("BOBILib evaluation: KKT (LP-relaxed follower) vs exact leader enumeration")
    print(f"Instances: {len(instance_list)} | Time limit: {TL}s | All SCIP")
    print("=" * 100)

    all_results = []
    for idx, entry in enumerate(instance_list):
        name = entry['name']
        aux_path = entry['aux_path']
        known_obj = entry['obj']

        mps_gz = aux_path.replace('.aux', '.mps.gz')
        try:
            mps = parse_mps_gz(mps_gz)
            lower_vars, lower_obj, lower_constrs = parse_aux(aux_path)
        except Exception as e:
            print(f"\n[{idx+1}] {name}: PARSE ERROR: {e}"); continue

        n_lint = len([v for v in lower_vars if v in mps['int']])
        print(f"\n[{idx+1}/{len(instance_list)}] {name} "
              f"(v={len(mps['cols'])} lv={len(lower_vars)} lint={n_lint} opt={known_obj})")

        row = {'name': name, 'known_obj': known_obj, 'methods': {}}

        for mname, fn in methods:
            try:
                r = fn(mps, lower_vars, lower_obj, lower_constrs)
            except Exception as e:
                r = {'status': 'crash', 'objective': None, 'time': 0}
                traceback.print_exc()

            obj = r.get('objective')
            gap = None
            if obj is not None and known_obj is not None:
                gap = abs(obj - known_obj) / max(1.0, abs(known_obj)) if abs(known_obj) > 1e-8 else abs(obj - known_obj)
            r['gap'] = gap
            row['methods'][mname] = r

            obj_s = f"{obj:.2f}" if obj is not None else "None"
            gap_s = f"{gap:.1e}" if gap is not None else "---"
            ok = "✓" if gap is not None and gap < 0.01 else ("✗" if gap is not None else "?")
            extra = ""
            for k in ('bilevel_feasible', 'proved_optimal', 'iterations'):
                if k in r: extra += f" {k}={r[k]}"
            print(f"  {mname:12s} {r['status']:12s} obj={obj_s:>12s} gap={gap_s:>9s} "
                  f"{ok} {r['time']:.2f}s{extra}")

        all_results.append(row)

    # Summary
    print("\n" + "=" * 100)
    mnames = [m[0] for m in methods]
    N = len(all_results)
    summary = {"instances": N, "time_limit_s": TL, "methods": {}}
    print(f"\n{'Method':12s} {'=opt':>6s} {'<=1%':>6s} {'proved':>7s} {'feasible':>9s} {'AvgTime':>8s}")
    for mn in mnames:
        rs = [r['methods'].get(mn, {}) for r in all_results]
        exact = sum(1 for x in rs if x.get('gap') is not None and x['gap'] <= 1e-6)
        within = sum(1 for x in rs if x.get('gap') is not None and x['gap'] <= 0.01)
        proved = sum(1 for x in rs if x.get('proved_optimal'))
        feas = sum(1 for x in rs if x.get('bilevel_feasible'))
        times = [x['time'] for x in rs if x.get('time') is not None]
        summary["methods"][mn] = {"matches_bobilib_optimum": exact, "within_1pct": within,
                                  "proved_optimal": proved, "checked_bilevel_feasible": feas,
                                  "mean_time_s": float(np.mean(times)) if times else None}
        print(f"{mn:12s} {exact:>6d} {within:>6d} {proved:>7d} {feas:>9d} {np.mean(times):>8.2f}")
    print(json.dumps(summary, indent=2))

    out_dir = os.path.join(script_dir, 'bobilib_results')
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, args.out), 'w') as f:
        json.dump({"summary": summary, "instances": all_results}, f, indent=2, default=str)
    print(f"\nSaved to {out_dir}/{args.out}")


if __name__ == '__main__':
    main()
