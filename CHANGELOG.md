# Changelog

## [0.2.0]

### Added
- `benchmarks/kkt_highs.py`: big-M KKT model of a bilevel program with an
  LP follower, with complementarity on constraint slacks and on reduced
  costs.
- `benchmarks/kkt_eval.py`: evaluation on 122 generated instances with a
  follower-LP answer check, leader enumeration for small instances and a
  big-M sweep.
- `benchmarks/bobilib_eval.py`: integer-follower evaluation on BOBILib with
  a checked KKT model of the follower's LP relaxation and an HPR-bounded
  leader enumeration that proves optimality.
- The BOBILib `.aux` files needed to read the instances.

### Changed
- `bicut compile` writes the high-point relaxation and says so; other
  `--reformulation` values are rejected.
- MPS output uses free format with full row and column names.

### Removed
- CLI subcommands that did not perform their task (`solve`, `analyze`,
  `benchmark`, `verify`, `generate`, `interactive`) and the `--certificate`
  flag.
- Simulated benchmark scripts and their outputs.

## [0.1.0]

- Rust workspace: problem types, LP/MPS readers and writers, reformulation
  passes, CLI.
