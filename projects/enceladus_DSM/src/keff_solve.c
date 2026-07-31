#include "keff.h"

/* ---------------------------------------------------------------------------
 * keff_solve.c -- assemble and solve the cell problem for all dim correctors.
 *
 * Step 4 of the merge deliberately hardwires a DIRECT solve (PREONLY + LU,
 * SuperLU_DIST in parallel). It is slow and does not scale, but it removes the
 * linear solver as a variable while the assembly and the tensor are being
 * validated against the analytic layered result. Step 5 switches the default to
 * CG + GAMG and validates it against this one on an identical assembly, so any
 * disagreement then is unambiguously the solver's fault.
 *
 * The KSP carries the "keff_" option prefix from the outset. Without it the
 * corrector would inherit inputs/solver.opts's -ksp_type bcgs -pc_type asm
 * -sub_pc_type ilu -sub_pc_factor_levels 3, i.e. BiCGStab + ILU(3) on an SPD
 * Laplacian -- slow, non-symmetric, and a plausible origin of the "iterative
 * solvers are unreliable for this problem" claim in the standalone project's
 * README, whose own header file simultaneously documents CG + GAMG.
 * ------------------------------------------------------------------------- */

PetscErrorCode KeffCreateSolver(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;
  PC             pc;

  PetscFunctionBegin;

  ierr = IGACreateKSP(kc->iga, &kc->ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(kc->ksp, "keff_"); CHKERRQ(ierr);

  ierr = KSPSetType(kc->ksp, KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(kc->ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
  {
    PetscMPIInt size;
    ierr = MPI_Comm_size(PetscObjectComm((PetscObject)kc->iga), &size); CHKERRQ(ierr);
    if (size > 1) { ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST); CHKERRQ(ierr); }
  }

  /* Set once, never per sample: PETSc tracks the matrix's object state itself,
   * and re-calling this each time only obscures the reuse intent. */
  ierr = KSPSetOperators(kc->ksp, kc->A, kc->A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kc->ksp); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffSolveCell
 *
 * Assemble the operator once, then solve dim right-hand sides against it.
 *
 * The pure-Neumann/periodic cell problem has a constant null space per
 * component. We pin global dof 0 with MatZeroRowsColumns rather than attaching a
 * MatNullSpace: pinning leaves the system non-singular for ANY solver, direct
 * included, and keeps the matrix symmetric. The pinned value is 0, which fixes
 * the otherwise arbitrary additive constant in t_m; k_eff depends only on
 * grad t_m, so the choice is immaterial to the answer.
 *
 * MatZeroRowsColumns is collective, but each rank passes only the rows IT owns
 * -- rank 0 passes the single global row 0, everyone else passes none. (The
 * standalone version had every rank pass the same global indices, which PETSc
 * tolerates but is not the documented contract.)
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffSolveCell(AppCtx *app, PetscInt *its_out, KSPConvergedReason *reason_out)
{
  PetscErrorCode ierr;
  KeffCtx       *kc  = app->keff;
  const PetscInt dim = kc->dim;
  PetscInt       rlo, rhi, npin, row0 = 0, its_total = 0;
  KSPConvergedReason worst = KSP_CONVERGED_ITERATING;

  PetscFunctionBegin;

  ierr = IGASetFormMatrix(kc->iga, KeffFormMatrix, kc); CHKERRQ(ierr);
  ierr = IGAComputeMatrix(kc->iga, kc->A); CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(kc->A, &rlo, &rhi); CHKERRQ(ierr);
  npin = (rlo <= 0 && 0 < rhi) ? 1 : 0;
  ierr = MatZeroRowsColumns(kc->A, npin, &row0, 1.0, NULL, NULL); CHKERRQ(ierr);

  ierr = IGASetFormVector(kc->iga, KeffFormVector, kc); CHKERRQ(ierr);

  for (PetscInt m = 0; m < dim; m++) {
    PetscInt           its;
    KSPConvergedReason reason;

    kc->cur_dir = m;
    ierr = IGAComputeVector(kc->iga, kc->b[m]); CHKERRQ(ierr);

    /* Match the pinned row: the constraint is t_m(dof 0) = 0. */
    if (npin) { ierr = VecSetValue(kc->b[m], row0, 0.0, INSERT_VALUES); CHKERRQ(ierr); }
    ierr = VecAssemblyBegin(kc->b[m]); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(kc->b[m]); CHKERRQ(ierr);

    /* kc->T[m] still holds the previous sample's corrector, which is an
     * excellent initial guess once the iterative solver is switched on and is
     * simply ignored by the direct solver. */
    ierr = KSPSolve(kc->ksp, kc->b[m], kc->T[m]); CHKERRQ(ierr);

    ierr = KSPGetIterationNumber(kc->ksp, &its); CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(kc->ksp, &reason); CHKERRQ(ierr);
    its_total += its;
    if (reason < 0) worst = reason;
    else if (worst >= 0) worst = reason;
  }

  if (its_out)    *its_out    = its_total;
  if (reason_out) *reason_out = worst;

  PetscFunctionReturn(0);
}

PetscErrorCode KeffDestroySolver(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;

  PetscFunctionBegin;
  if (kc && kc->ksp) { ierr = KSPDestroy(&kc->ksp); CHKERRQ(ierr); }
  PetscFunctionReturn(0);
}
