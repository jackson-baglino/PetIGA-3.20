# k_eff corrector solver benchmark

**Why.** In batch 3a the k_eff solves took 55–75% of every run's wall time.
Measured cost per sample and per iteration:

| run | ranks | s/sample (median) | s per CG iteration |
|---|---|---|---|
| batch 2 | 241 | 65 | 0.49 |
| batch 3a | 401 | 108 | 0.85 |
| rev64 (2026-09-18) | 615 | 4 | 0.034 |

- **Iteration counts are similar everywhere** (~110–140), so the problem is
  cost per iteration, not convergence.
- **It gets worse with MORE ranks**, the signature of a solve that is
  communication- or setup-bound rather than work-bound. One candidate: a
  multigrid coarse level that stops shrinking, solved on one rank every
  iteration.

Shakedown cost: $195 against a ~$60–100 estimate. Projected over the
remaining 118 runs at this rate: ~$4k.

**What.** Replay 5 snapshots of the 3a φ = 0.475 run under each
configuration, with `-log_view`. No new simulation; each job is ~10–20 min.

| label | change from production |
|---|---|
| `base` | none: CG + GAMG, rtol 1e-10 (the reference) |
| `rtol7` | rtol 1e-7. Symmetry already shows k_eff exact to ~1e-11; seed scatter is ~3% |
| `gamgthr` | rtol 1e-7 + GAMG drop threshold 0.02 (helps high-contrast/anisotropic coarsening) |
| `hypre` | rtol 1e-7 + BoomerAMG (linked into this PETSc build) |
| `freeze` | rtol 1e-7 + hold the multigrid hierarchy, rebuilding every 10 samples |
| `base100k` | `base` at 100k DoF/core (241 ranks): the rank-count question on its own |

The submit commands are in the headers of `bench.txt` and `bench_100k.txt`.

**Read.**

```bash
venv_enceladus/bin/python studies/keff_sintering/solver_benchmark/analyze_bench.py \
    <local copy of the 0.475 run dir> <local copy of the bench batch dirs>
```

Fetch the run dir's `bench_*.csv` and each bench batch's `*.o*` logs first.

**Decision rule.** Adopt the fastest configuration whose
max |Δk|/k_iso against `base` is below 1e-5. Record it in
`scripts/HPC/submit_keff_production.sh` before batch 3b, and note in `TODO.md`
that 3a used `base`. The difference is far below anything reported, and
replaying the 3a snapshots with the new setting can show it. If nothing
clears the bar, keep `base` and just revert the allocation.
