# LOG.md — running session log (newest first)

## 2026-11-05/06 — cryosuction stage 0b/0c, model-freeze-v2
- Session: 99675a49-73a2-42ba-ab59-ef7bd012f37d
- Done: analyzed stage 0b (v1): the do-nothing open-face phi BC is ill-posed
  (film 1.66 eps); crashes came from under-resolution plus the soft CFL and the
  shared +-0.05 guard. xi fixed at Moure & Fu values (no sweeps, no pre-solve).
  Ported the enceladus safeguards and tagged lunar_cryosuction/model-freeze-v2
  (physics band [-0.01, 1.01], solver guard [-0.01, 1.05], enforced
  interface-CFL rollback, dphimax 0.05, reservoir-face watch, 2 tau_sub
  start-up window for the resolution checker; gates in
  studies/cryosuction/verification/verify_safeguards.sh and verify_resolution.sh).
  All experiments use -phi_open_bc 0. Started the run ledger
  studies/cryosuction/RUNS.md. Stage 0c (v2, 5 runs) all completed clean:
  B1 harvests all ice in 2.3 days; C3/C4 change the harvest rate 1.25x/0.84x,
  matching rho_vs at the ice face (sublimation cooling 0.017 K matters);
  C2 shows shielding. Added studies/cryosuction/scripts/tgrad_envelope.py
  (temperature gradient vs suction: cold-trapping length l_T).
- Remaining (next session: the user presents the conference and manuscript
  simulation plan):
  1. Throat resolution: narrow-plug / open-side cap phases at Pi ~ 1 on the
     0.86 um meshes. Add an eps argument to
     preprocess/build_geometry_two_throat.py (EPS constant at line 118), then
     a 1-day eps-convergence pair at eps 0.328 um (comp_eps --Lx 8.5e-4
     --Ly 1.6e-4 --Rave 1.5e-5 --T0 -20 --alpha 1e-2 --vn 1.95e-9) against 0c-C1.
  2. B1 series law: compare with compare_velocity.py (rough fit: slope 1.15x
     1D diffusion, intercept ~0.5 beta0/sigma).
  3. Plan rules: every experiment carries -phi_open_bc 0 and keeps ice off
     the open face; -rhovfix_lo is a fraction of rho_vs(temp0), not of the
     local face T; no nucleation, so seed ice where deposition is expected;
     thermal BC must be stated per run; eps check at both T extremes once
     dT across the domain > ~1 K.
  4. Minor: -t_final 0 with output hangs; vapor_movie.py piles tracers at x = 0.
  5. Five local commits not yet pushed (ask first).
- Bottlenecks: refined throat mesh is ~7x the DOFs (~12,700 steps for 5 days).

## 2026-11-04 — project created from lunar_regolith_DSM
- Session: b7a604f2-0658-4583-934e-c4bd21b9422b (run from lunar_regolith_DSM)
- Done: froze the two-phase lunar model at tag `lunar_cryosuction/model-freeze-v1`
  (147a49f0); copied the solver side here (6ebaab34); renamed binary, main
  file, run scripts, venv and results paths (52e7ba16); pinned igakit to
  upstream commit 92ee0975, the code venv_lunar runs; added FREEZE.md and
  scripts/check_freeze.sh (edec892c); project CLAUDE.md, README, memory,
  settings and wrap-up skill.
- Remaining: push to origin (user's call), then pull, build and create
  venv_cryosuction on the HPC; plan the first paper-2 study.
- Bottlenecks: none.
