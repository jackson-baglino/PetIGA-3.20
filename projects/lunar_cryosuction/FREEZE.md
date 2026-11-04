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
| — | — | — | (no changes since the freeze) | — |
