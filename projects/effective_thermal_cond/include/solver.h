#ifndef SOLVER_H
#define SOLVER_H

#include "app_ctx.h"

/*
 * ComputeInitialCondition — zero the solution vector before each solve.
 */
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user);

/*
 * CreateSolverObjects — allocate and configure the matrix A, RHS vector b,
 * and KSP solver ksp. Call once before the per-file loop.
 *
 * The KSP is configured with CG + GAMG (or HYPRE BoomerAMG if available).
 * It can be further tuned via PETSc command-line options.
 */
PetscErrorCode CreateSolverObjects(IGA iga, AppCtx *user,
                                   Mat *A, Vec *b, KSP *ksp);

/*
 * Solve — reassemble A and b for the current ice field, reset the KSP
 * operators, and solve.  Enforces zero-mean per component (parallel-safe).
 * Call once per file inside the loop.
 */
PetscErrorCode Solve(AppCtx *user, IGA iga, Mat A, Vec b, KSP ksp);

/*
 * DestroySolverObjects — free A, b, and ksp after the loop.
 */
PetscErrorCode DestroySolverObjects(Mat *A, Vec *b, KSP *ksp);

#endif /* SOLVER_H */

