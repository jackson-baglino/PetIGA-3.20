---
paths:
  - "src/**"
  - "include/**"
---

# C solver code

- C99. Declare loop variables in the `for` statement. No variable-length
  arrays; use `PetscMalloc`/`PetscFree`.
- `CHKERRQ(ierr)` immediately after every PETSc/MPI call.
- Prefer `PetscMin`, `PetscMax`, `PetscSqrtReal`, `PetscTanhReal` over bare C
  math.
- Lint with `scripts/lint_c.sh <file>`. A Jacobian change is checked with
  `-snes_test_jacobian` passed through the run script.

# The model is frozen (FREEZE.md)

- Physics and defaults do not change: no edits to governing equations,
  constitutive laws, existing residual or Jacobian terms, or option defaults.
  New physics belongs in `lunar_regolith_DSM`.
- Allowed: new ICs, geometry, opt-in BC or forcing options (off by default),
  output, bug fixes (ask first if a fix changes results). Commit each change
  by itself, add its row to the `FREEZE.md` ledger, and run
  `scripts/check_freeze.sh` until it prints OK.
- An opt-in addition must leave runs without it bit-for-bit unchanged; check
  with the G1 gate of `studies/contact_angle/verification/verify_wall_bc.sh`.
