# LOG.md — running session log (newest first)

## 2026-11-13 — Slate slide template (AGU talk and defense)
- Session: ffbfac7f-55b4-489f-86d5-ae0ff9f85390
- Done: Designed a PowerPoint template for the AGU talk and later the defense.
  Concepts A-D, then four variations of B; the user chose B option 1 (burgundy
  accent, slate text, Franklin Gothic) with the frost background (#EEF2F5) from
  option 3. Built `Slate.potx` (15 layouts; official Caltech and JPL logos per
  identity.caltech.edu; autofit off), a 19-slide example deck from the
  2026-11-05 update, and `tools/update_outline.py` (uv script: draws the
  section-outline row from PowerPoint Sections, fills the Outline slide,
  numbers dividers; "_" sections excluded, "Full | Short" labels). Checked:
  validator passes; PowerPoint-measured text fit, 0 overflow in a 40-slide,
  408-box stress deck. Everything is in iCloud
  `Documents/Research/SlideTemplate/` (builder scripts in `source/`); the
  .potx is installed in PowerPoint's personal Templates folder. Nothing in
  this repository changed.
- Remaining: user to confirm Slate appears under File > New from Template >
  Personal, and to try a real figure in a figure layout. Then build the AGU
  deck on it.
- Bottlenecks: LibreOffice hangs here; renders and fit checks go through
  PowerPoint via AppleScript (scripts in SlideTemplate/source/).

## 2026-11-06/07 — simulation plan, analysis tools and tests, first production batches
- Session: be49e5a1-66f2-4f54-9dcf-073a3f0d6550
- Done: SimulationPlan_v2 (21 pp.) and AGU_SimulationPlan (4 pp., chi for
  curvature) in the iCloud SimulationProposalsSetup folder (build scripts in
  SimulationPlan_v2_build/). Built plan tools T1-T8 (ice_bodies, water_budget,
  section_flux, rhov_bc, check_open_bc in preflight, type_summary, vertical
  display, compare_velocity --at-ell/--at-r) plus figcheck.py and a spline
  snapshot loader in pplib; fixed pplib.read_opts precedence (solver.opts had
  been winning). Test suites: postprocess/tests 233, studies/cryosuction/tests
  71, preprocess/tests 28, all passing; every bug they found in the older tools
  was fixed (gt_balance temperature, meniscus_velocity split/mapped/restart,
  plot_mass weights, crossings on a node, preflight last line, ...). Closed
  plan issue 2 (0c-B1 slope 0.996-1.000x with D_v(-20 C), area-based face) and
  issue 4. Added STATUS.md (imported by CLAUDE.md, updated in wrap-up). Built
  the issue-1 fine mesh (eps 0.328 um) and all type 3-4 inputs
  (make_type34_inputs.py). Pushed. User submitted I1-fine (200k DoF/core),
  T3 (10 runs) and T4 (7 runs) at TARGET_DOFS_PER_CORE=50000. DoF/core:
  enceladus' 200k was set by its k_eff solve and 24 M-DoF multi-node meshes;
  on the 0.4-0.5 M-DoF lunar meshes, Sept runs (8-10 ranks, file timestamps,
  checked against SLURM to 2 s) ran 1.4-2.0 s/step against 4.2-7.0 s/step at
  2-3 ranks in Oct, at similar core-hours. The velocity-study README's
  "3.2 s/step at 8 ranks" is not reproduced by any record.
- Remaining:
  1. Download and post-process I1-fine, T3, T4 (RUNS.md rows): resolution.dat
     Pi, rollbacks (watch wedge om20, element aspect ~8 at its narrow end),
     s/step from the .o files (decide on a DoF/core scaling test after that;
     enceladus has submit_scaling_test.sh to adapt). I1: narrow-plug vanishing
     time vs 0c-C1 day 1 (< 5 % keeps 0.86 um). T3/T4: ice_bodies,
     water_budget, compare_velocity --at-ell / ice_bodies --at gap.
  2. Geometry builders: wall caps of a given theta (types 1, 2, first), single
     throat (5), dog-bone (6), 2 mm regolith pore with ice layout (7); then
     experiment files and the shakedown batch.
  3. User decisions open: D1-D4, D6 (assumed as the plan suggests), hypotheses.
  4. Offered, not done: fix the velocity-study README's 3.2 s/step passage
     (lunar_regolith_DSM); the full plan still uses kappa in the iCloud copy
     (the builder now writes chi; the user has edited SimulationPlan_v2.docx,
     so it was not overwritten).
- Bottlenecks: geometry builders for types 1, 2, 5, 6, 7; wall time of the
  180-day type-4 runs.

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
