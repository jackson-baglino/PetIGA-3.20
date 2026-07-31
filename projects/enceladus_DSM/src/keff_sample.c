#include "keff.h"

/* ---------------------------------------------------------------------------
 * keff_sample.c -- the one function both the in-line and replay paths call.
 *
 * KeffSample is deliberately the whole of the shared core: project phi, solve
 * the cell problem, form the tensor, report. The in-line path gets (step, t, U)
 * from the TS monitor; the replay path gets it from the filesystem. Nothing
 * else differs, so the two cannot drift apart.
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * KeffCSVAppend
 *
 * One row per sample, through a persistent viewer opened on first use and
 * flushed each time. That is the user->ssa_view pattern from monitoring.c, which
 * exists because the standalone code's per-call fopen/fclose silently exhausted
 * file descriptors around step 15 (the calls had no error checking), truncating
 * the file while the run carried on.
 *
 * The schema carries `time`, which the old sol_index-keyed CSV lacked. Without
 * it k_eff(t) cannot be plotted, nor joined to SSA_evo.dat, without a separate
 * lookup table. phi_bar is free (one accumulator in the integrand) and is the
 * x-axis of every Hashin-Shtrikman/Wiener bound plot. ksp_its and ksp_reason
 * make a degrading frozen AMG hierarchy visible in the data file rather than
 * only in a log nobody keeps; wall_s is what sizes -keff_freq for the next run.
 * ------------------------------------------------------------------------- */
static PetscErrorCode KeffCSVAppend(AppCtx *app, PetscInt step, PetscReal t,
                                    const PetscReal keff[9], PetscReal phi_bar,
                                    PetscReal k_iso, PetscInt its,
                                    KSPConvergedReason reason, PetscReal wall)
{
  PetscErrorCode ierr;
  KeffCtx       *kc  = app->keff;
  const PetscInt dim = kc->dim;

  PetscFunctionBegin;

  if (!kc->csv) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, kc->csv_path, &kc->csv); CHKERRQ(ierr);
    /* Truncate, not append: one file per run, no "does it exist" branch. */
    ierr = PetscViewerFileSetMode(kc->csv, FILE_MODE_WRITE); CHKERRQ(ierr);

    ierr = PetscViewerASCIIPrintf(kc->csv, "step,time"); CHKERRQ(ierr);
    for (PetscInt i = 0; i < dim; i++)
      for (PetscInt j = 0; j < dim; j++) {
        ierr = PetscViewerASCIIPrintf(kc->csv, ",k_%d%d", (int)i, (int)j); CHKERRQ(ierr);
      }
    ierr = PetscViewerASCIIPrintf(kc->csv,
             ",phi_bar,k_iso,ksp_its,ksp_reason,wall_s\n"); CHKERRQ(ierr);
  }

  ierr = PetscViewerASCIIPrintf(kc->csv, "%d,%.12e", (int)step, (double)t); CHKERRQ(ierr);
  for (PetscInt q = 0; q < dim * dim; q++) {
    ierr = PetscViewerASCIIPrintf(kc->csv, ",%.12e", (double)keff[q]); CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPrintf(kc->csv, ",%.12e,%.12e,%d,%d,%.3f\n",
                                (double)phi_bar, (double)k_iso, (int)its,
                                (int)reason, (double)wall); CHKERRQ(ierr);
  ierr = PetscViewerFlush(kc->csv); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffDue -- is this step a sampling step?
 *
 * Mirrors OutputMonitor's cadence logic (monitoring.c:302-316) deliberately, so
 * the two behave the same way, INCLUDING the while-loop advance. A single
 * `t_next += t_interv` degenerates to every-step sampling as soon as dt exceeds
 * t_interv, because t_next then falls permanently behind t -- the same bug
 * OutputMonitor carries a comment about.
 * ------------------------------------------------------------------------- */
static PetscBool KeffDue(KeffCtx *kc, PetscInt step, PetscReal t)
{
  if (step == 0) return kc->at_step0;

  if (kc->t_interv > 0.0) {
    if (t < kc->t_next) return PETSC_FALSE;
    while (kc->t_next <= t) kc->t_next += kc->t_interv;
    return PETSC_TRUE;
  }
  return (PetscBool)(kc->freq > 0 && step % kc->freq == 0);
}

/* ---------------------------------------------------------------------------
 * KeffMonitor -- the TS hook.
 *
 * Registered AFTER Monitor (enceladus_main.c), and the order matters. Monitor
 * sets user->bounds_violated when phi leaves [phase_lo, phase_hi] and returns
 * early; BoundsRollbackPreStep then discards that step entirely. PETSc runs
 * monitors in registration order, so registering last lets this one see the flag
 * and skip. Otherwise a full corrector solve is spent on a state that is about
 * to be thrown away, AND a physically meaningless row lands in the CSV --
 * ThermalCond clamps phi to [0,1], so the number would look perfectly plausible
 * while corresponding to no state the run ever kept.
 *
 * The TSMonitor read-lock on ts->vec_sol is not a constraint here: this callback
 * only reads U, and everything it mutates (kc->ice, A, b, T) it privately owns.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx)
{
  PetscErrorCode ierr;
  AppCtx        *app = (AppCtx *)mctx;
  KeffCtx       *kc  = app->keff;

  PetscFunctionBegin;
  (void)ts;
  if (!kc || !kc->enabled) PetscFunctionReturn(0);
  if (app->bounds_violated) PetscFunctionReturn(0);
  if (!KeffDue(kc, step, t)) PetscFunctionReturn(0);

  ierr = KeffSample(app, step, t, U); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

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

  ierr = KeffCSVAppend(app, step, t, keff, phi_bar, k_iso, its, reason,
                       (PetscReal)(t1 - t0)); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
