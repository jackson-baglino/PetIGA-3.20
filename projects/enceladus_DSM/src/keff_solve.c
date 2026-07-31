#include "keff.h"

/* ---------------------------------------------------------------------------
 * keff_solve.c -- assemble and solve the cell problem for all dim correctors.
 *
 * DEFAULT: CG + GAMG. After pinning one dof the operator is symmetric positive
 * definite, and it is a plain scalar Laplacian with a piecewise-smooth
 * coefficient -- textbook algebraic multigrid territory. Direct LU stays one
 * flag away (-keff_ksp_type preonly -keff_pc_type lu) and is what the analytic
 * benchmark was first validated with; studies/snow_thermal/verification holds
 * the reference numbers.
 *
 * The standalone project's README asserted that iterative solvers are unreliable
 * for this problem class, while its own solver.h simultaneously documented
 * CG+GAMG and its code set PREONLY+LU. The likeliest explanation for the
 * original failure is the missing option prefix: without one the corrector KSP
 * inherits inputs/solver.opts's -ksp_type bcgs -pc_type asm -sub_pc_type ilu
 * -sub_pc_factor_levels 3, i.e. BiCGStab + ILU(3) on an SPD Laplacian. Hence
 * KSPSetOptionsPrefix below, before any KSPSetFromOptions.
 *
 * Why LU cannot stay the default: the study meshes reach Nx=Ny=2829, i.e. ~8M
 * scalar unknowns per corrector. A direct factorization of that, per sample, for
 * thousands of samples, is not viable in time or memory.
 *
 * AMG REUSE. Between consecutive samples the interface moves a fraction of a
 * cell, so the previous hierarchy is an excellent preconditioner for the next
 * solve. Three things make that pay off, and one of them is easy to miss:
 *   - MAT_KEEP_NONZERO_PATTERN on the matrix (set in KeffCreate). Without it
 *     MatZeroRowsColumns REMOVES the zeroed off-diagonals, which bumps
 *     A->nonzerostate; PCSetUp then sees DIFFERENT_NONZERO_PATTERN and redoes
 *     the full symbolic setup every sample, making reuse a no-op.
 *   - PCGAMGSetReuseInterpolation: keep aggregates and prolongators, rebuild
 *     only the Galerkin coarse operators from the current values. Cheap, and
 *     stays correct as phi evolves.
 *   - KSPSetInitialGuessNonzero: kc->T[m] still holds the previous sample's
 *     corrector, which for a slowly evolving microstructure is close.
 * A harder freeze (-keff_pc_freeze) additionally skips the coarse-operator
 * rebuild, refreshing only every -keff_pc_refresh samples, with an adaptive
 * backstop so a degrading hierarchy heals itself instead of silently poisoning
 * the tensor.
 * ------------------------------------------------------------------------- */

PetscErrorCode KeffCreateSolver(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;
  PC             pc;

  PetscFunctionBegin;

  ierr = IGACreateKSP(kc->iga, &kc->ksp); CHKERRQ(ierr);

  /* CRITICAL, and must precede KSPSetFromOptions: see the header comment. */
  ierr = KSPSetOptionsPrefix(kc->ksp, "keff_"); CHKERRQ(ierr);

  ierr = KSPSetType(kc->ksp, KSPCG); CHKERRQ(ierr);
  /* The true residual, not the preconditioned one: with AMG the preconditioned
   * norm can be optimistic by orders of magnitude, and this tolerance is what
   * the k_eff values inherit. */
  ierr = KSPSetNormType(kc->ksp, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);
  ierr = KSPSetTolerances(kc->ksp, 1.0e-10, PETSC_DEFAULT, PETSC_DEFAULT, 500); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(kc->ksp, PETSC_TRUE); CHKERRQ(ierr);

  ierr = KSPGetPC(kc->ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCGAMG); CHKERRQ(ierr);
  ierr = PCGAMGSetType(pc, PCGAMGAGG); CHKERRQ(ierr);
  ierr = PCGAMGSetReuseInterpolation(pc, PETSC_TRUE); CHKERRQ(ierr);

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

  /* -keff_pc_freeze: hold the whole preconditioner across samples, refreshing
   * only every -keff_pc_refresh. Off by default -- reuse-interpolation alone is
   * the safe win, since it still rebuilds the coarse operators from the current
   * coefficient. */
  if (kc->pc_freeze) {
    PC        pc;
    PetscBool refresh = (PetscBool)(kc->pc_refresh > 0 &&
                                    kc->nsamples % kc->pc_refresh == 0);
    ierr = KSPGetPC(kc->ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetReusePreconditioner(pc, (PetscBool)!refresh); CHKERRQ(ierr);
  }

  for (PetscInt m = 0; m < dim; m++) {
    PetscInt           its;
    KSPConvergedReason reason;

    kc->cur_dir = m;
    ierr = IGAComputeVector(kc->iga, kc->b[m]); CHKERRQ(ierr);

    /* Match the pinned row: the constraint is t_m(dof 0) = 0. */
    if (npin) { ierr = VecSetValue(kc->b[m], row0, 0.0, INSERT_VALUES); CHKERRQ(ierr); }
    ierr = VecAssemblyBegin(kc->b[m]); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(kc->b[m]); CHKERRQ(ierr);

    /* kc->T[m] still holds the previous sample's corrector, which for a slowly
     * evolving microstructure is a good initial guess (and is ignored by the
     * direct solver). */
    ierr = KSPSolve(kc->ksp, kc->b[m], kc->T[m]); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(kc->ksp, &its); CHKERRQ(ierr);
    ierr = KSPGetConvergedReason(kc->ksp, &reason); CHKERRQ(ierr);

    /* Adaptive backstop: a stale hierarchy shows up as a stalled or diverged
     * solve. Rebuild it and retry once, rather than writing a quietly wrong
     * tensor row. */
    if (reason < 0 || (kc->max_its > 0 && its >= kc->max_its)) {
      PC pc;
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "  [keff] corrector %d: %d its, reason %d -- rebuilding the preconditioner "
        "and retrying\n", (int)m, (int)its, (int)reason); CHKERRQ(ierr);
      ierr = KSPGetPC(kc->ksp, &pc); CHKERRQ(ierr);
      ierr = PCSetReusePreconditioner(pc, PETSC_FALSE); CHKERRQ(ierr);
      ierr = PCSetUp(pc); CHKERRQ(ierr);
      ierr = VecZeroEntries(kc->T[m]); CHKERRQ(ierr);
      ierr = KSPSolve(kc->ksp, kc->b[m], kc->T[m]); CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(kc->ksp, &its); CHKERRQ(ierr);
      ierr = KSPGetConvergedReason(kc->ksp, &reason); CHKERRQ(ierr);
      kc->nrebuilds++;
    }

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
