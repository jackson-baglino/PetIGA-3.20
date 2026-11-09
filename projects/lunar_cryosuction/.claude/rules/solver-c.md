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

# The solver is unfrozen until production (FREEZE.md, user 2026-11-09)

- Solver changes are allowed, including physics and defaults; discuss changes
  that alter existing results with the user first. Never cite the freeze as a
  reason not to do something.
- Commit each solver change by itself and add its row to the `FREEZE.md`
  ledger (change log), saying whether runs that do not use it are affected.
- A change meant to be inert for existing runs must leave them bit-for-bit
  unchanged; check with the G1 gate of
  `studies/contact_angle/verification/verify_wall_bc.sh`.
