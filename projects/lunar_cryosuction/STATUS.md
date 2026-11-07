# STATUS.md: where paper 2 and the AGU talk stand

Read this at the start of every session; update it whenever the state below
changes and always during wrap-up. Keep it to what is true now: move finished
items to "Done recently" (with a date), drop them from there after about two
weeks (LOG.md keeps the history), and never let it pass about 150 lines.
Run-level detail lives in `studies/cryosuction/RUNS.md`; session history in
`LOG.md`; the plan itself in the iCloud folder below.

**Last updated:** 2026-11-06 (session 99675a49 / be49e5a1)

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

**Closing the open issues, building inputs, starting HPC runs** (plan
timeline: issues closed by Mon 9 Nov; builders and experiment files 9-11 Nov;
shakedown Wed 11 Nov; production from 12-14 Nov).

## Decisions

| # | Decision | State |
|---|---|---|
| D1 | theta_ref = 60 deg | assumed (plan suggestion); confirm |
| D2 | type 7: 2 mm, h_open 0.99 | assumed; confirm |
| D3 | new thermal BC options after AGU | assumed |
| D4 | AGU types: 1, 2, 6, 7 certain; 3-5 as trends | assumed |
| D5 | accept B1 resolution once checked | pending the check (issue 2) |
| D6 | type 8 in paper 3 | assumed |
| D7 | push local commits before the shakedown | **ask the user before pushing** |
| — | curvature symbol is chi, not kappa | decided 2026-11-06 |
| — | open phi boundary dropped: every experiment sets -phi_open_bc 0 | decided 2026-11-06 |

Hypotheses: the user writes them; none recorded yet.

## Open issues (plan section 5)

| # | Issue | State | Next step |
|---|---|---|---|
| 1 | Fast phases under-resolved on 0.86 um meshes (Pi ~ 1) | open | `--eps` for build_geometry_two_throat.py, then a 1-day pair at eps 0.328 um vs 0c-C1 (user submits) |
| 2 | B1 series-law slope 1.15x | likely resolved (D_v(T): ratio 1.004) | confirm on the 0c-B1 output with ice_bodies + series fit |
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
  `preprocess/tests` (28); commands in postprocess/README.md. Bugs found in the
  older tools are being fixed (each was an expectedFailure test).
- Still to build: geometry builders for the single throat (type 5), dog-bone
  (type 6), 2 mm regolith pore with its ice layout (type 7), wall caps of a
  given theta (types 1, 2); experiment files for types 1-7.

## HPC

- Running: none. Queued: none.
- Next submission: issue-1 eps pair, then shakedown (one short run per type).
- 24+ local commits are not pushed (D7); the HPC needs them for new inputs.

## Done recently

- 2026-11-06: SimulationPlan_v2 (21 pp.) delivered; tools T1-T8 and tests;
  pplib.read_opts precedence fix (runs had been read with solver.opts winning).
- 2026-11-05: stage 0c complete on model-freeze-v2 (RUNS.md).
