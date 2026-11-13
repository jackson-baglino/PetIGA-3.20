# STATUS.md: where paper 2 and the AGU talk stand

Read this at the start of every session; update it whenever the state below
changes and always during wrap-up. Keep it to what is true now: move finished
items to "Done recently" (with a date), drop them from there after about two
weeks (LOG.md keeps the history), and never let it pass about 150 lines.
Run-level detail lives in `studies/cryosuction/RUNS.md`; session history in
`LOG.md`; the plan itself in the iCloud folder below.

**Last updated:** 2026-11-12 (session 6475a664)

## Goal and deadline

- AGU 2026 presentation, final by **Mon 30 Nov 2026**; then paper 2.
- Plan: `SimulationPlan_v2.docx` (full, 21 pp.) and `AGU_SimulationPlan.docx`
  (4 pp.) in
  `~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/lunar_regolith_DSM/Prompts/SimulationProposalsSetup/`.
  The build scripts are in `SimulationPlan_v2_build/` there.
- Eight simulation types: 1 droplet, 2 several droplets, 3 lens in a flat pore,
  4 wedge, 5 throat, 6 constriction, 7 regolith pore under mining conditions,
  8 lunar physics (paper 3, not this model).

## Phase

**First production batches (T3, T4, I1) finished and analyzed; next: build the
remaining geometries** (plan timeline: builders and experiment files 9-11 Nov;
shakedown Wed 11 Nov; production from 12-14 Nov).

## Decisions

| # | Decision | State |
|---|---|---|
| D1 | theta_ref = 60 deg | assumed (plan suggestion); confirm |
| D2 | type 7: 2 mm, h_open 0.99 | assumed; confirm |
| D3 | new thermal BC options after AGU | assumed |
| D4 | AGU types: 1, 2, 6, 7 certain; 3-5 as trends | assumed |
| D5 | accept B1 resolution | check done (issue 2 closed); confirm |
| D6 | type 8 in paper 3 | assumed |
| D7 | push local commits | approved and pushed 2026-11-07; ask again for later pushes |
| — | curvature symbol is chi, not kappa | decided 2026-11-06 |
| — | open phi boundary dropped: every experiment sets -phi_open_bc 0 | decided 2026-11-06 |

Hypotheses: the user writes them; none recorded yet.

## Open issues (plan section 5)

| # | Issue | State | Next step |
|---|---|---|---|
| 1 | Fast phases under-resolved on 0.86 um meshes (Pi ~ 1) | **closed** 2026-11-07 (user): I1-fine within 3 % (area 0.31 %, vanishing time 3.0 %; REPORT_issue1); production stays on eps 0.8584 um | check Pi in each new geometry's log |
| 2 | B1 series-law slope 1.15x | **closed** 2026-11-06: slope 0.996-1.000x with D_v(-20 C), area-based face (RUNS.md 0c-B1) | D5 can be accepted |
| 3 | Thermal BCs set the harvest rate | rule adopted | every run table has a thermal column |
| 4 | Default -phi_open_bc 1 | **closed** | check_open_bc.sh runs in preflight.sh (commit 4db30788) |
| 5 | Collector is h_open x rho_vs(T0) | rule adopted | rhov_bc.py prints exact values; methods text |
| 6 | No nucleation | rule adopted | seed ice where deposition is expected |
| 7 | -20 C analog | rule adopted | framing with advisor |

## Tools (plan section 8)

- Built and tested 2026-11-06: T1 ice_bodies, T2 water_budget, T3 section_flux,
  T4 rhov_bc, T5 check_open_bc, T6 type_summary, T7 vertical display (pplib
  helpers, used by T1/T6), T8 compare_velocity --at-ell/--at-r; figcheck.py.
- Tests: `postprocess/tests` (233), `studies/cryosuction/tests` (71),
  `preprocess/tests` (28), all passing; commands in postprocess/README.md. Every
  bug the tests found in the older tools was fixed on 2026-11-06.
- Types 3 and 4 inputs done (studies/cryosuction/scripts/make_type34_inputs.py).
- Type 6 inputs done 2026-11-08 (build_geometry_dogbone.py, make_type6_inputs.py;
  type6_tests.txt, preflight OK).
- Solver UNFROZEN by the user 2026-11-09 until production (FREEZE.md status;
  keep logging solver commits there as a change log).
- Type 7: build_ice_layout.py + -ic_type phi_file (044a2cf4) built and gated
  2026-11-09; layout A (2 mm rough pore, 3 lenses at 60 deg, 4 droplets and
  3 grains at a 40 deg water placeholder, 24.6 % ice, 2.41 M DOF) and the 7A
  experiment are in inputs; 7B-7D experiment files still to write.
- Types 1, 2, 5 inputs done 2026-11-09 (make_type125_inputs.py; 19 + 14 + 12
  runs; caps at 40 deg water placeholder, user; lenses at 60 deg). Per-run
  preview PDFs: preprocess/make_run_previews.py -> preprocess/previews/type<N>/.
- Still to write: experiment files for 7B-7D; S2, S5.
- 2026-11-13 (branch periodic-x, not merged): solver `-periodic_x 1` (44f8e40f) and periodic
  converging-diverging cells in build_ice_layout.py (`"periodic": true`, wall profile sinhyp/cos)
  for the redesigned type 5; verify_periodic_x.sh 56/56 PASS, G1 PASS (studies/periodic_x/README.md).

## HPC

- T6 ran 2026-11-09/11: 5 of 7 finished; wc 40 G 300 (theta 60, 90) will hit the
  48 h limit near day 10-12 of 30 (RUNS.md). Lens splits into wall films;
  only the 40 um neck delays transit. Overview of every type:
  `studies/cryosuction/RUN_PROGRESS.md`. User (11 Nov): G 100/300 too large
  (summer runs <= 50 K/m). Run plan 2026-11-11 proposes G 1-10 (+30), S7 scan
  first; T1/T2 gradient runs and T6/T7C-D on hold until approved.
- 2026-11-12: 6B-1 restart prepared in place (`resume_run.sh`, solver option
  -restart_append eae44490); 6B-2 (wc 80, theta 90, G 300) inputs built. User:
  G 1-3 K/m too slow; plan_estimates.py section 7 gives the G window (wide-gas
  G below ~18.5 K/m for 150 um faces; below ~44 K/m runs are dtmax-bound, wall
  time ~ 1/G). User chose G 10-20 K/m with a few runs at 50 (2026-11-12);
  RunPlan_2026-11-12 (7 pp.) uses them and adds 6D (theta-90 face stability,
  straight pore, flat face + seed); inputs to rebuild at the new G.
  6B-1 voided by the user (faces started 0.375 um out of flat); its resume is
  withdrawn. Flat-face reruns 6B-1/6B-2 (G 300, 401 snapshots) prepared
  (type6b_flat_tests.txt, solver 5815f290), to go --into the type-6 batch.
- 2026-11-13: ready to queue (preflight OK, pushed): type 1 (18), type 2
  (14), 6B flat (2), 6C (9), 6D (6). New opt-in solver options: event
  snapshots -outp_dphi (e5a593e7), flat-face seed -ice_flat_seed (fe8a7334).
  Type 7: six candidate layouts (A-F) for the user to choose on 13 Nov
  (review_scratch/2026-11-13_type7_layouts). No eps halving (user).
  I1-fine, T3, T4 finished 2026-11-07 (RUNS.md). T6 lenses start
  at the wall-BC angle (60 deg concave via the new opt-in IC option
  -bridge_concaveL/R, FREEZE 18270f6f; 6B-1 flat). Droplets/caps start at
  liquid water's angle on regolith: no source in Zotero (user); an effective
  value from the mineral mix is to be engineered (user, 2026-11-09). Lit
  search 2026-11-09: only lunar value is Zhang et al. 2025 (Chang'E-5,
  10.1016/j.apsusc.2025.163331), 108.7 deg apparent on rough particles
  (Cassie-Baxter); minerals ~15-45 deg; Cassie 1948 (10.1039/df9480300011)
  rule applied to mineral fractions in Arismendi Florez et al. 2024
  (10.1016/j.jmrt.2024.07.177, mercury). No pyroxene/glass data. Near-0 deg
  becomes a stated limitation (user). Realistic-pore ice: equal-curvature
  rule (draft shown 2026-11-09).
- DoF/core: decide on a scaling test after seeing these runs' s/step
  (user, 2026-11-07). Compare with Oct 2026 2-3-rank runs: 4.2-7.0 s/step.
  T3 result (job wall time / steps, includes setup): 1.5-2.2 s/step at
  47-50k DoF/core (11 ranks for 0.52 M DoF, 5 for 0.25 M), 17-24 core-s/step
  for 0.52 M DoF; in line with the Sept 2026 50k runs (1.5-2.0 s/step).
  T4: 1.3-1.7 s/step at 42-50k DoF/core (5-8 ranks, 0.22-0.40 M DoF).
- No shakedown batch (user, 2026-11-09): cheap runs are submitted, inspected and
  rerun; pre-check only jobs above ~30-40 cores (type 7, 2.5 M DOF).
- Pushed 2026-11-07 (D7 approved).
- Submit (from the project root on the HPC, after `git pull`):
  `./scripts/HPC/submit_batch.sh --tag <tag> --tests-file studies/cryosuction/<file> --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs -- --time=<t>`
  with I1: issue1_eps_pair, 2-00:00:00 (7.6 M DOF: default 200k DoF/core);
  T3: type3 and T4: type4, each `TARGET_DOFS_PER_CORE=50000` and 2-00:00:00.
  Evidence (2026-11-07, same 387k/495k-DoF meshes): Sept 2026 runs at 8-10
  ranks (~50k/core) 1.5-2.0 s/step, 12-20 core-s/step (wall time from file
  timestamps, a method that matches SLURM to 2 s); Oct 2026 runs at 2-3 ranks
  (~190k/core) 4.2-7.0 s/step, 8.5-21 core-s/step (SLURM). Same core-hours
  within ~20 %, ~3.4x shorter wall time at 50k. The velocity README's
  "3.2 s/step at 8 ranks" is not reproduced by any record. Not a controlled
  test (Sept used alpha_c 2e-3 and older code); confirm with a 2- vs 8-rank
  pair in the shakedown. The enceladus 200k optimum is for 24 M DoF meshes.

## Done recently

- 2026-11-07: T4 (type 4, 7 runs) analyzed: a taper moves ice to the narrow
  end in every wedge; series law with theta -/+ omega to 3-5 % (beta_eff 8.8e5
  in the flat pore, as in T3). Every theta 60 lens reached a vapor face (20 deg
  at day 8): future margins must be longer. I1-fine: 0.86 um mesh within 3 %
  of the fine mesh. Reports REPORT_type4, REPORT_issue1 in the batch folders.
  New: analyze_type4.py, compare_issue1.py.
- 2026-11-07: contact_angle.py fixed: keeps the remaining face when one face
  reaches the domain edge (was: all lost); reports contact-line speed and
  calls faces static (< 2 eps travel) or moving; extrapolates theta_inf only
  for static faces (was: theta_inf -6406 deg on a moving wedge). Moving faces
  show a dynamic angle below Young (3A-1: 25 vs 30 deg); static ones match to
  0.06 deg. Reports T3/T4 updated.

- 2026-11-07: T3 (type 3, 10 runs) finished and analyzed: series law holds to
  1 % with nothing fitted; report `REPORT_type3.md/.pdf` in the batch folder.
  Open: moving-face contact angle (fit) is off Young but the speed follows Young.
  New tools: run_health.py, md2pdf.py (batch reports), analyze_type3.py.
  water_budget fix 2db33eb8.
- 2026-11-07: I1-fine, T3, T4 submitted; types 3-4 inputs; DoF/core evidence (LOG.md).
- 2026-11-06: AGU_SimulationPlan (4 pp.) delivered; issue 2 closed; tool bugs fixed.
- 2026-11-06: SimulationPlan_v2 (21 pp.) delivered; tools T1-T8 and tests;
  pplib.read_opts precedence fix (runs had been read with solver.opts winning).
- 2026-11-05: stage 0c complete on model-freeze-v2 (RUNS.md).
