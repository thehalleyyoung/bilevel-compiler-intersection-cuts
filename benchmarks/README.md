# Benchmarks

Requirements: Python 3.10 or later with `numpy`, `highspy` and `pyscipopt`.

| Script | What it does | Output |
|---|---|---|
| `kkt_eval.py` | 122 generated instances with integer leader and LP follower (`generators.py`). Big-M KKT model (`kkt_highs.py`) solved with HiGHS at M = 1000, follower-LP check of every answer, exact optimum by leader enumeration where there are at most 4096 leader vectors, and a sweep over M = 1, 2, 5, 10, 50, 100, 1000, 10^4. | `kkt_eval_output/kkt_eval.json` |
| `bobilib_eval.py` | BOBILib instances listed in `bobilib/instance_list.json`. KKT model of the follower's LP relaxation (SCIP, indicator constraints) checked with the follower MIP, and the HPR-bounded leader enumeration warm-started from the checked KKT value. Every objective is compared with the optimum BOBILib reports. | `bobilib_results/bobilib_eval.json` |

```bash
python3 benchmarks/kkt_eval.py
python3 benchmarks/bobilib_eval.py --max 80 --time-limit 60
```

`bobilib/` holds the BOBILib `.mps.gz` and `.aux` files and the solution
files with the reported optima, taken from
<https://gitlab.uni-trier.de/nlopt/bobilib-data>.
