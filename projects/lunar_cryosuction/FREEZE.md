# FREEZE.md — the frozen model for paper 2

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
