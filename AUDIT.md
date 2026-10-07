# Internal audit (2026-10-07)

This file records what the October 2026 audit found and changed. It is not
part of the paper or the README.

## Claims ledger

| Claim | Where | Verdict | Evidence | Action |
|---|---|---|---|---|
| "5x faster than MibS" / "~5x geomean speedup" | README (old), benchmarks/README, run_benchmarks.py | UNSUPPORTED | `run_benchmarks.py` and `miplib_benchmark.py` draw times, gaps, nodes and cuts from `rng.gauss` around hand-set difficulty scalars, "calibrated to match paper claims"; no solver is called | Scripts removed, claim removed |
| "18% root gap closure" | old README/paper, run_benchmarks.py | UNSUPPORTED | same simulation | Removed |
| "2,600+ BOBILib instances" / "564 BOBILib instances" / BOBILib harness | CHANGELOG, planning docs, paper appendix | UNSUPPORTED | 80 instances are listed in `benchmarks/bobilib/instance_list.json`; none ran from a clean clone because `.gitignore` excluded every `.aux` file | `.aux` files added; 80 instances now evaluated |
| "first bilevel intersection cuts" / "first tool providing bilevel-specific cutting planes" / "first to apply intersection cuts to bilevel feasibility" | old README, paper abstract and related work | UNSUPPORTED (false) | Fischetti, Ljubić, Monaci, Sinnl, IPCO 2016 (LNCS 9682:77-88), Oper. Res. 65(6):1615-1637 (2017), Math. Prog. 172:77-103 (2018) | Removed; these works are cited as related work |
| "first bilevel compiler with correctness certificates" | paper abstract | UNSUPPORTED | `--certificate` wrote a hard-coded `"valid": true`; library pipeline returns `certificate: None` | Flag removed, claim removed |
| Bilevel intersection cuts: polyhedrality, polynomial separation, finite convergence theorems | old paper Sec. 2, App. A-C | UNSUPPORTED | Proposition gives O(m^d T_LP), its own proof gives O(T_LP + n m); the convergence proof counts vertices that grow with each cut; no measured run uses the Rust cut code (the Rust B&C log `experiments_output.txt` shows 0 cuts on every instance) | Removed from paper |
| Compiler soundness theorem part (iii) "every cut is valid" | old paper App. D | UNSUPPORTED (contradicted) | the Python value-function cuts made 15 of 122 instances infeasible | Removed |
| "Value-function cuts are valid after the sign fix" | README, paper Sec. 4 | FALSE | cut was `c^T y >= lam^T (b + Bx)` with `lam = -duals`; not implied by follower optimality. With the sign corrected it is weak duality, implied by primal feasibility, so it can never cut | Removed with the script |
| "KKT reformulation is mathematically sound; 121/122 solved" | README, HONEST_EVALUATION, paper Table 1 | FALSE as stated | all three HiGHS KKT builders lacked y-complementarity; `honest_evaluation_v2` also had the dual sign wrong (`A^T lam >= c`), so on families with `c < 0` the model was the high-point relaxation. `certificate_benchmark`'s own check flagged 247 of 329 "optimal" answers as follower-suboptimal; every Stackelberg instance returned 0.0 while enumeration gives e.g. 42.86 | New `kkt_highs.py`; rerun gives 119/122, all checked, 34/34 equal to enumeration |
| "LP relaxation is surprisingly tight (0% on 82-84 instances)" | HONEST_EVALUATION, paper Finding 1 | FALSE | artefact of the KKT bug; with the correct model the KKT LP gap is 0 only on the knapsack family (shared objective), 47% mean on dense, 100% on target defence | Replaced by measured numbers |
| Big-M table (dense_15x15 M=10 gives -32 "wrong") | HONEST_EVALUATION, paper Table 3 | FALSE as stated | produced by the buggy model; correct model: see `kkt_eval.json` sweep | Replaced by full sweep |
| "Certificate verification catches 100% of bad-big-M errors" | certificate_benchmark.py docstring | FALSE | its own run: 1 of 43 errors caught (2.3%) | Script removed; Proposition 1 in the paper shows why a follower check cannot catch a bad M in a correct KKT model |
| Integer-follower "VF-cut loop solved 26/34 correctly" | integer_follower_results.json | UNSUPPORTED | no reference optimum was computed; the no-good row assumed binary x; the loop returned the current HPR value instead of the best recorded incumbent; the KKT baseline had the slack-row sign error | Removed |
| BOBILib "Combined 17/20 (85%)" | BOBILIB_RESULTS.md, v4_results.json | FIXABLE | reproduced 17/20 within 1% (16 exact) with the committed code, but BilevelBnB's no-good row is invalid for general integers (14/20 instances), which is why moore90_2 gave 6 instead of 5 | Replaced by exact enumeration; see results below |
| "bicut compile -r kkt / strong-duality / value-function / ccg" | README, examples/README, CLI help | FALSE | every choice wrote the same high-point relaxation (byte-identical files) | CLI now writes and labels the HPR; other choices are refused |
| "bicut solve / analyze / benchmark / verify / generate / interactive", `info`, `check-cq`, `--backend` | README, examples/README, CLI help | FALSE | stubs printed "... complete."; `SolveCommand` returned objective 0.0, 100 nodes; `info`/`check-cq`/`--backend` never existed | Removed |
| "Three solver backends Gurobi/SCIP/HiGHS", indicator/SOS1 emission, LRU/adaptive cut cache with >90% hit rate, parallel cut generation with Rayon, sampling oracle with L1 error bounds | README features/config | UNSUPPORTED | not reachable from the CLI; no run measures them | Removed from README |
| "302 tests pass (10+110+108+74)"; HONEST_EVALUATION "63 compilation errors", "6 LP failures" | README, HONEST_EVALUATION | FIXABLE | measured: 718 tests pass in `cargo test --workspace` after removing the non-compiling example; per crate types 108, core 80, lp 110, value-function 74, cuts 10, compiler 103+1, branch-cut 91, bench 87, cli 54 | README states 718 |
| examples/interdiction.toml, strategic_bidding.toml | README examples | FALSE | first fails to parse; second's objectives are in tables the parser ignores | Removed |
| Dual MIT/Apache-2.0 license with LICENSE-MIT/LICENSE-APACHE | README, Cargo.toml, CONTRIBUTING | FALSE | only an MIT LICENSE exists | Changed to MIT |
| Repository URL github.com/bicut-project/bicut | README, Cargo.toml, CONTRIBUTING | FALSE | HTTP 404 | Changed |
| Citation as INFORMS JOC 2025 / Math. Prog. 2025 articles | README | FALSE | no such publications | Replaced by @software entry |
| Pineda & Morales, "Solving linear bilevel problems using big-Ms and their connection to penalty approaches", INFORMS JOC 31(4):739-751 | paper bib | FABRICATED | Crossref: the paper is "Solving linear bilevel problems using big-Ms: not all that glitters is gold", IEEE TPWRS 34(3):2469-2471, 2019, and it argues the opposite of what it was cited for | Corrected |
| BilevelJuMP, IJOC 34(4):1990-2008, 2022 | paper bib | FIXABLE | Crossref: IJOC 36(2):327-335, 2024 | Corrected |
| Balas & Margot, Math. Prog. 137(1):19-41 | paper bib | FIXABLE | pages 19-35 | Dropped (no longer cited) |
| ~64,000 lines of Rust | paper | VERIFIED approx. | 68,610 lines in `crates/*/src` | Not stated in the new paper |
| Planning documents (ideation/, docs/design/, proposals/, problem_statement.md, theory/, State.json) | repo | UNSUPPORTED | pre-implementation targets stated as properties | Removed |

## Bugs fixed

1. KKT models in the Python benchmarks (missing reduced-cost complementarity,
   dual sign, slack-row sign, missing bound multipliers). Replaced by
   `benchmarks/kkt_highs.py`.
2. BOBILib no-good row invalid for general-integer leaders and ranging over
   follower variables; incumbent value not optimistic. Replaced by the exact
   enumeration in `benchmarks/bobilib_eval.py`.
3. `.gitignore` rule `*.aux` excluded the BOBILib `.aux` files.
4. CLI: reformulation flag ignored; fabricated certificate, solve, analyze and
   interactive output; fixed-format MPS truncating names to 8 characters
   (duplicate rows `follower`).
5. `examples/real_benchmark_runner.rs` did not compile, so
   `cargo test --workspace` failed.

## Measured results after the fixes

- `cargo test --workspace`: 718 passed, 0 failed.
- `kkt_eval.py`: see `benchmarks/kkt_eval_output/kkt_eval.json`.
- `bobilib_eval.py --max 80`: see `benchmarks/bobilib_results/bobilib_eval.json`.

Summary of the reruns:

| Run | Before the audit (committed files) | After |
|---|---|---|
| Generated LP-follower suite | "121/122 KKT solved"; Stackelberg objectives all 0.0; 247/329 KKT answers in certificate_benchmark not follower-optimal | 119/122 optimal, 119/119 pass the follower check, 34/34 equal the enumerated optimum |
| Value-function cuts | 35 cuts, 0 helped, 15 instances made infeasible (reproduced exactly) | removed (cut not implied by follower optimality) |
| Big-M sweep | old table from the buggy model | 976 runs; 527 returned points, all follower-optimal; M <= 5 always infeasible, M = 10: 40 correct, 70 infeasible, 7 suboptimal |
| BOBILib, first 20 | v4 "Combined" 17/20 within 1%, 16 exact (reproduced) | enumeration proves 13, final value equals optimum on 18; KKT checked feasible 19, optimal 9 |
| BOBILib, 80 listed | not runnable (no .aux files) | enumeration proves 35 (all agree with BOBILib), equals optimum on 49; KKT checked feasible 65, optimal 16, infeasible 8 |

## Not verified

- Wall-clock times: all runs were on a machine with load average 100-300
  from other jobs.
- The Rust reformulation passes (KKT, strong duality, value function, CCG)
  pass their unit tests but were not checked end to end against a solver;
  `bicut_compiler::compile` returns an `LpProblem` without integrality
  markers, so its output as a MILP was not evaluated.
- The Rust intersection-cut, value-function and branch-and-cut crates were
  not evaluated beyond their unit tests.
