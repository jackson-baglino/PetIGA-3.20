# STATUS.md: where paper 2 and the AGU talk stand

Read this at the start of every session; update it whenever the state below
changes and always during wrap-up. Keep it to what is true now: move finished
items to "Done recently" (with a date), drop them from there after about two
weeks (LOG.md keeps the history), and never let it pass about 150 lines.
Run-level detail lives in `studies/cryosuction/RUNS.md`; session history in
`LOG.md`; the plan itself in the iCloud folder below.

**Last updated:** 2026-11-07 (session 6475a664)

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

**First production batches running; next: post-process them and build the
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
| 1 | Fast phases under-resolved on 0.86 um meshes (Pi ~ 1) | inputs ready (RUNS.md I1-fine) | user submits the 1-day eps 0.328 um run; compare with 0c-C1 day 1 |
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
- Still to build: geometry builders for the single throat (type 5), dog-bone
  (type 6), 2 mm regolith pore with its ice layout (type 7), wall caps of a
  given theta (types 1, 2); experiment files for types 1, 2, 5, 6, 7.

## HPC

- Running (submitted 2026-11-07): I1-fine (1 run), T4 (7 runs); T3 finished;
  RUNS.md rows. Nothing else is ready to submit.
- DoF/core: decide on a scaling test after seeing these runs' s/step
  (user, 2026-11-07). Compare with Oct 2026 2-3-rank runs: 4.2-7.0 s/step.
  T3 result (job wall time / steps, includes setup): 1.5-2.2 s/step at
  47-50k DoF/core (11 ranks for 0.52 M DoF, 5 for 0.25 M), 17-24 core-s/step
  for 0.52 M DoF; in line with the Sept 2026 50k runs (1.5-2.0 s/step).
- Next: shakedown (one short run per type) once builders and experiment files exist.
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
