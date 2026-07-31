#include "keff.h"

/* ---------------------------------------------------------------------------
 * keff_sample.c -- the one function both the in-line and replay paths call.
 *
 * KeffSample is deliberately the whole of the shared core: project phi, solve
 * the cell problem, form the tensor, report. The in-line path gets (step, t, U)
 * from the TS monitor; the replay path gets it from the filesystem. Nothing
 * else differs, so the two cannot drift apart.
 * ------------------------------------------------------------------------- */

PetscErrorCode KeffSample(AppCtx *app, PetscInt step, PetscReal t, Vec U)
{
  PetscErrorCode     ierr;
  KeffCtx           *kc = app->keff;
  PetscReal          keff[9] = {0,0,0,0,0,0,0,0,0};
  PetscReal          phi_bar = 0.0, k_iso = 0.0;
  PetscInt           its = 0;
  KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  PetscLogDouble     t0, t1;

  PetscFunctionBegin;
  if (!kc || !kc->enabled) PetscFunctionReturn(0);

  ierr = PetscTime(&t0); CHKERRQ(ierr);

  ierr = KeffProjectIce(app, U); CHKERRQ(ierr);
  ierr = KeffSolveCell(app, &its, &reason); CHKERRQ(ierr);
  ierr = KeffComputeTensor(app, keff, &phi_bar); CHKERRQ(ierr);

  ierr = PetscTime(&t1); CHKERRQ(ierr);

  for (PetscInt d = 0; d < kc->dim; d++) k_iso += keff[d * kc->dim + d];
  k_iso /= (PetscReal)kc->dim;

  kc->nsamples++;

  if (reason < 0)
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  [keff] WARNING: corrector solve did not converge (KSPConvergedReason %d) "
      "at step %d\n", (int)reason, (int)step); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  [keff] step %-6d t = %.6e s   phi_bar = %.6f   k_iso = %.6e W/m/K"
    "   (%d its, %.2f s)\n",
    (int)step, (double)t, (double)phi_bar, (double)k_iso,
    (int)its, (double)(t1 - t0)); CHKERRQ(ierr);

  if (kc->dim == 2) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "         k = [ %.9e  %.9e ]\n"
      "             [ %.9e  %.9e ]\n",
      (double)keff[0], (double)keff[1],
      (double)keff[2], (double)keff[3]); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "         k = [ %.9e  %.9e  %.9e ]\n"
      "             [ %.9e  %.9e  %.9e ]\n"
      "             [ %.9e  %.9e  %.9e ]\n",
      (double)keff[0], (double)keff[1], (double)keff[2],
      (double)keff[3], (double)keff[4], (double)keff[5],
      (double)keff[6], (double)keff[7], (double)keff[8]); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}
