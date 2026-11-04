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
