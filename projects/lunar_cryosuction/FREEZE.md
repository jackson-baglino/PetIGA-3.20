# FREEZE.md — the frozen model for paper 2

## Status: UNFROZEN (user, 2026-11-09)

The user unfroze the solver on 2026-11-09: until production runs start, the
solver (equations, residual/Jacobian, ICs, BCs, defaults) may be changed as
the work needs. The freeze returns, with a new tag, when production starts;
the user decides when. Until then this file is a change log: add a ledger row
for every solver commit so each run can be tied to the code it used, and say
in the row whether the change alters results of runs that do not use it.
`scripts/check_freeze.sh` still reports unlogged solver commits; its OK is
not a gate while unfrozen. Runs made so far (T3, T4, I1) used
`model-freeze-v2` plus the logged additions.


## What is frozen

The paper-2 results must all come from one model, the two-phase lunar solver as
it stood when the thesis work was split (2026-11-04):

- **Tag:** `lunar_cryosuction/model-freeze-v1`, on commit `147a49f0` of
  `projects/lunar_regolith_DSM`. The last solver change before the tag is
  `a7af16ed` (2026-11-01, rename of `-t_out` to `-t_interv`).
- **Split commits:** `6ebaab34` (verbatim copy) and `52e7ba16` (names and
  paths). After them, `src/`, `include/` and `inputs/solver.opts` are
  byte-identical to the tag, except that `src/lunar_main.c` is renamed
  `src/cryosuction_main.c`.

The model keeps developing in `lunar_regolith_DSM` (paper 3). Changes made
there do not come here.

## Versions

| Tag | Date | What changed | Runs on it |
|---|---|---|---|
| `lunar_cryosuction/model-freeze-v1` | 2026-11-04 | the split (above) | stage 0 (cancelled) and stage 0b (`batch_2026-11-04__19.34.44_cryo_stage0b`), with the logged v1 additions 60bd33a5, 770a9bb5, 8c96ab09 |
| `lunar_cryosuction/model-freeze-v2` | 2026-11-05 | solver safeguards ported from enceladus_DSM, unfrozen by the user (ledger rows 118d69f2, a47d659d, d9c10138). Equations, constitutive laws and the residual/Jacobian are unchanged; what changed is which steps are accepted. | everything after stage 0b |

Why v2: in stage 0b, five of the eight coupled runs aborted. Their solver
logs and fields match the failure modes that enceladus_DSM had already fixed
(see the ledger rows). The review of 2026-11-05 found no regression in the
physics: four independent comparisons covered the BC code, enceladus runs,
lunar rectangular runs and lunar distorted-mesh runs.

## The rule: physics frozen, setup open

**Frozen.** The governing equations, the constitutive laws
(`material_properties.c`), the residual and Jacobian of the existing terms
(`assembly.c`), and every solver default. Defaults live in `inputs/solver.opts`
and in the default values of the options parsed in `src/cryosuction_main.c`.

**Allowed**, each logged below:
- **IC:** a new `-ic_type` or IC option.
- **geometry:** reading or building a new geometry.
- **BC option:** a new boundary condition or forcing, for example for
  cryosuction-like setups. It must be opt-in and off by default.
- **output:** monitoring, diagnostics, file output.
- **bug fix:** a fix that changes results is still a physics question. Ask the
  user before committing it, and record which runs it affects.

An opt-in addition must leave runs that do not use it unchanged. The check for
that is the regression gate in
`studies/contact_angle/verification/verify_wall_bc.sh` (G1). It builds a binary
from a base commit and compares `SSA_evo.dat` bit for bit on a short run
without the new option. Record the result in the ledger row.

Changing the physics or a default means unfreezing. Only the user decides
that, and it gets a new tag (`…-v2`), a ledger row, and a note on which runs
belong to which version.

## How to log a change

1. Commit the solver change by itself.
2. Add a row below with its 8-character hash (`git log -1 --format=%h --abbrev=8`),
   then commit the row.
3. `scripts/check_freeze.sh` must print `check_freeze: OK` before any run whose
   results are meant for the paper.

## Ledger

| Commit | Date | Category | What and why | Regression gate |
|---|---|---|---|---|
| 60bd33a5 | 2026-11-04 | BC option | `-rhovfix_faces`, `-Tfix_faces`: pin vapor or T only on named faces (open end plus dead pore). Unset = legacy logic. `-wall_faces` parser factored into `ParseFaceList`. | G1–G4 PASS (`verify_wall_bc.sh`, base 1443a0b0); F1 bitwise vs 1443a0b0 with legacy reservoirs, F2–F4 PASS (`studies/cryosuction/verification/verify_face_bc.sh`) |
| 770a9bb5 | 2026-11-04 | output | `bflux.dat`: physical (unscaled D_v, k) vapor and heat flux out through every Dirichlet face, one row per step; `-bflux_output 0` disables. Read-only. | G1–G4 PASS; H1 bitwise SSA vs 8484a311, H2–H3 analytic 1D flux within 1e-5 and independent of xi_v, H4 PASS (`studies/cryosuction/verification/verify_bflux.sh`) |
| 8c96ab09 | 2026-11-04 | output | Resolution checker: per-step K&P Eq. 45 number Pi = beta0 v_n eps/d0 from the measured front speed (boxed warning if Pi > 1), ice on open faces, start-up mesh/eps check, end summary, `resolution.dat`. Read-only; `-resolution_check 0` disables. | G1–G4 PASS; K1 bitwise SSA vs 1e608260, K2 v_n within +6.4 % of the area speed (conservative by design), K3–K5 warnings fire and stay quiet when they should (`studies/cryosuction/verification/verify_resolution.sh`) |
| 118d69f2 | 2026-11-05 | unfreeze v2 (default) | Separate solver guard from the phi physics band (port of enceladus c6dee625/60244164): `-phase_lo/hi` physics band [-0.01, 1.01] (rollback), `-snes_guard_lo/hi` [-0.01, 1.05] on phi only (domain error). Before: one band +-0.05 on phi and phi_a did both; stage-0b B1 stalled with phi_max = 1.0500 at the guard edge. | S1 bitwise vs 8d659590 with the legacy values on the command line (all SSA columns but the logged next-dt); S3 banner warning (`studies/cryosuction/verification/verify_safeguards.sh`) |
| a47d659d | 2026-11-05 | unfreeze v2 (behavior) | Interface-CFL limit enforced by rollback (port of enceladus caa209e0), monitor registered first, all other monitors skip a condemned step; BoundsRollbackPreStep re-anchors cfl_U_prev/res_U_prev (fix beyond the enceladus code). | S1 as above; S2 forced violation rolls back once, kept steps max \|dphi\| 0.045 <= 0.05, no duplicate step rows; S4 completes |
| d9c10138 | 2026-11-05 | unfreeze v2 (default) | `-dtCFL_dphimax` default 0.2 -> 0.05 (enceladus 60244164 ladder); cryosuction experiment files drop their explicit 0.2. | S1 unaffected (sets the value explicitly); S2 runs at 0.05 |
| 57d4e469 | 2026-11-05 | output (v2) | Resolution checker watches every vapor-reservoir face for ice (`res_watch_face`), not only faces open for phi, so the warning keeps working with `-phi_open_bc 0`; message renamed RESERVOIR-FACE WARNING. | K1 bitwise SSA vs HEAD^, K2–K5 PASS, new K6 (warns with -phi_open_bc 0) PASS (`verify_resolution.sh`) |
| 60460bd1 | 2026-11-05 | output (v2) | Resolution checker: Pi not judged inside a start-up window of `-resolution_skip_tau` x tau_sub (default 2); rows still written; the summary reports the window and its max Pi. | K1 bitwise SSA vs HEAD^, K2-K6 PASS (K3 now with -resolution_skip_tau 0: its run is shorter than its window), new K7 PASS: 20 vs 90 warnings with and without the window, first at t = 721 s >= 719 s (`verify_resolution.sh`) |
| 18270f6f | 2026-11-09 | IC | `-bridge_concaveL/R` (one 0/1 per bridge, default 0): a `-bridge_*` plug face becomes concave (ice outside the face disc, beyond its centre), so lenses can start at a wall-BC angle below 90 deg (user, 2026-11-09). Unset = the legacy convex IC. | G1–G4 PASS (`verify_wall_bc.sh`, base dc601583); C1 convex bridge bitwise vs dc601583 and explicit 0 bitwise vs unset, C2 IC area of a concave (convex) 60 deg lens within 0.02 % of the sharp shape, C3 guards PASS (`studies/cryosuction/verification/verify_concave_bridge.sh`) |
| 044a2cf4 | 2026-11-09 | IC (unfrozen) | `-ic_type phi_file -phi_file <path>`: ice field per control point from `preprocess/build_ice_layout.py` (PHIIC001 binary: x, y, phi), for type-7 layouts (equal-curvature lenses on rough walls, droplets, free grains). T and rho_v as every other IC. Runs without it unchanged. | G1 PASS (`verify_wall_bc.sh`, base f1b1738d); P1 IC area within 0.08 % of the sharp layout, P2 phi in [0, 1], P3 four guards (`studies/cryosuction/verification/verify_phi_file.sh`); layout A at 2.41 M DOF: area 6.0790e-8 vs 6.0783e-8 (0.01 %) |
| eae44490 | 2026-11-12 | restart (unfrozen) | `-restart_append` with `-initial_cond sol_N.dat`: step counter resumes at N, the per-step logs are appended (no new header, no row for the loaded step), sol_N.dat is not rewritten, %-change baseline from the first leg's step 0. Used by `scripts/HPC/resume_run.sh` (RESUME_DIR mode of `scripts/HPC/run_cryosuction.sh`, which trims the logs to step <= N after a backup). The generalized-alpha integrator restarts from the state alone, so the first steps differ slightly from an uninterrupted run. Runs without it unchanged. | G1 PASS (`verify_wall_bc.sh`, base d2acff29); local test on a copy of 6B-1 from sol_09165: rows continue at step 9166 under the original headers, no repeated row, sol_09170 written under the continued numbering; first step dt 150 s against 204 s in the original and 6 Newton iterations against 1; ice area at step 9169 within 3e-6 of the original interpolated to the same t |
| 5815f290 | 2026-11-12 | IC guard (unfrozen) | `-ic_type multi_grains` accepts `-ice_flat_x/-ice_flat_R/-ice_flat_height` as the only ice body (the guard required grains, wedge bands or bridges) and prints each flat body. Used for theta-90 lenses with exactly flat faces (build_geometry_dogbone.py, 'lens90flat'), replacing the crossed 7.5 mm discs whose faces sat 0.375 um out of flat (6B-1 voided by the user). Runs that do not use it unchanged. | G1 PASS (`verify_wall_bc.sh`, base 82296b53); local 2-step check on dogbone wc 40 lens90flat: phi = 1/2 at x = 875.000108 and 1024.999628 um on all 497 node rows (spread 0 to 1e-9 um), phi in [0, 1] |
