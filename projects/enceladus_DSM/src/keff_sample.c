#include "keff.h"
#include <dirent.h>
#include <string.h>

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
static PetscBool KeffDueSteps(KeffCtx *kc, PetscInt step, PetscReal t)
{
  if (kc->t_interv > 0.0) {
    if (t < kc->t_next) return PETSC_FALSE;
    while (kc->t_next <= t) kc->t_next += kc->t_interv;
    return PETSC_TRUE;
  }
  return (PetscBool)(kc->freq > 0 && step % kc->freq == 0);
}

/* The SSA trigger (-keff_dlnssa). WHY: sampling on a step count spends samples
 * evenly in steps, but k_eff bends fastest early, so a fixed stride leaves
 * visible corners there and wastes samples on the nearly straight late curve.
 * Sampling on relative SSA change places them where the k-SSA curve -- the
 * figure -- changes. On batch 2's every-step data (studies/keff_sintering/
 * sampling_stride.py follow-up, 2026-09-29) a 0.3% trigger gave 0.06% max
 * interpolation error with 90 samples at -20 C, against 0.22% for every 2
 * steps with 164. And since temperature acts as a time rescaling, the same
 * option samples every temperature at the same states.
 *
 * Before dlnssa_t0_tau*tau_sub (the IC relaxation, shown greyed and never a
 * baseline) the step cadence applies. The reference ln SSA is the last
 * SAMPLE, so steps between samples accumulate the drop. A restart records its
 * first step as the reference without sampling it. */
static PetscBool KeffDue(KeffCtx *kc, PetscInt step, PetscReal t)
{
  AppCtx   *app = kc->app;
  PetscBool due;

  if (kc->dlnssa <= 0.0) {
    if (step == 0) return kc->at_step0;
    return KeffDueSteps(kc, step, t);
  }

  const PetscBool have_ssa = (PetscBool)(app->ssa_step == step && app->ssa_now > 0.0);
  if (!have_ssa) {
    /* Monitor did not run this step (-pf_monitor 0): no SSA to trigger on. */
    if (!kc->warned_no_ssa) {
      PetscPrintf(PETSC_COMM_WORLD, "  [keff] WARNING: -keff_dlnssa needs the per-step "
                  "SSA from Monitor (-pf_monitor 1); falling back to -keff_freq\n");
      kc->warned_no_ssa = PETSC_TRUE;
    }
    if (step == 0) return kc->at_step0;
    return KeffDueSteps(kc, step, t);
  }
  const PetscReal ln_ssa = PetscLogReal(app->ssa_now);
  const PetscReal tau    = app->tau_sub_run;

  if (step == 0) {
    due = kc->at_step0;
  } else if (!kc->have_last) {
    /* first step after a restart: this state is the reference, not a sample */
    kc->have_last = PETSC_TRUE; kc->ln_ssa_last = ln_ssa; kc->t_last = t;
    return PETSC_FALSE;
  } else if (tau <= 0.0 || t < kc->dlnssa_t0_tau * tau) {
    due = KeffDueSteps(kc, step, t);
  } else {
    due = (PetscBool)((kc->ln_ssa_last - ln_ssa) >= kc->dlnssa
                      || (kc->max_gap_tau > 0.0 && t - kc->t_last >= kc->max_gap_tau * tau));
  }
  if (due || step == 0) {
    kc->have_last = PETSC_TRUE; kc->ln_ssa_last = ln_ssa; kc->t_last = t;
  }
  return due;
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

  /* -keff_write_corrector: the corrector fields t_m of this sample, one file
   * per macroscopic direction m, beside the CSV:
   *     igakeff.dat             the scalar corrector IGA (written once)
   *     t_vec_<step>_<m>.dat    corrector for direction m (0 = x, 1 = y, ...)
   * Read with igakit: PetIGA().read("igakeff.dat"), then read_vec(). The
   * local temperature for a unit macroscopic gradient e_m is x_m + t_m, and
   * the flux is -K (e_m + grad t_m). For figures and checks; ~64 MB per
   * direction at the campaign mesh, so it is off by default. */
  if (kc->write_corrector) {
    char  dir[PETSC_MAX_PATH_LEN], fn[PETSC_MAX_PATH_LEN];
    char *slash;
    ierr = PetscStrncpy(dir, kc->csv_path, sizeof(dir)); CHKERRQ(ierr);
    slash = strrchr(dir, '/');
    if (slash) *slash = '\0';
    else { ierr = PetscStrncpy(dir, ".", sizeof(dir)); CHKERRQ(ierr); }
    if (kc->nsamples == 0) {
      ierr = PetscSNPrintf(fn, sizeof(fn), "%s/igakeff.dat", dir); CHKERRQ(ierr);
      ierr = IGAWrite(kc->iga, fn); CHKERRQ(ierr);
    }
    for (PetscInt m = 0; m < kc->dim; m++) {
      ierr = PetscSNPrintf(fn, sizeof(fn), "%s/t_vec_%05d_%d.dat", dir, (int)step, (int)m); CHKERRQ(ierr);
      ierr = IGAWriteVec(kc->iga, kc->T[m], fn); CHKERRQ(ierr);
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD, "  [keff] wrote correctors t_vec_%05d_{0..%d}.dat in %s\n",
                       (int)step, (int)kc->dim - 1, dir); CHKERRQ(ierr);
  }

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

/* ---------------------------------------------------------------------------
 * REPLAY MODE  (-keff_replay <run_dir>)
 *
 * Computes k_eff for every sol_*.dat in a finished run directory, replacing the
 * standalone binary that used to do this.
 *
 * The decisive simplification over that binary: because the replaying process
 * builds its IGA from the SAME options files as the producing run, user->iga IS
 * the IGA that wrote the snapshots. So a snapshot is read with a plain
 * IGAReadVec onto the live 3-dof vector, and handed to the identical KeffSample
 * the in-line path uses.
 *
 * The old code instead reconstructed an input IGA from re-parsed
 * -Nx/-p/-C/-periodic and compared vector lengths to catch a mismatch. Its own
 * comment (field_init.c:56-58) admitted that check could not detect a -periodic
 * mismatch -- which is precisely the mismatch that silently yields wrong
 * answers, since a periodic and a non-periodic mesh of the same Nx have
 * different node counts only for some degrees. Reusing the producer's IGA makes
 * the entire failure class unrepresentable.
 * ------------------------------------------------------------------------- */

/* Rank 0 scans <dir> for sol_NNNNN.dat and returns the sorted step indices.
 * Two fixes over the original: the count is taken first and the array heap
 * allocated (the original used a fixed PetscInt[10000] stack array, 80 KB, that
 * silently truncated), and the parse is "sol_%d.dat" rather than "sol_%05d.dat"
 * so it keeps working past 99999 steps. */
static PetscErrorCode KeffScanSnapshots(const char *dir, PetscInt **steps_out,
                                        PetscInt *n_out)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       n = 0, *steps = NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    DIR           *d;
    struct dirent *ent;

    d = opendir(dir);
    if (!d) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                    "-keff_replay: cannot open directory '%s'", dir);

    while ((ent = readdir(d)) != NULL) {
      int  idx; char tail[16];
      if (sscanf(ent->d_name, "sol_%d.dat%15s", &idx, tail) == 1 && idx >= 0) n++;
    }
    rewinddir(d);

    if (n > 0) {
      PetscInt k = 0;
      ierr = PetscMalloc1(n, &steps); CHKERRQ(ierr);
      while ((ent = readdir(d)) != NULL) {
        int idx; char tail[16];
        if (sscanf(ent->d_name, "sol_%d.dat%15s", &idx, tail) == 1 && idx >= 0 && k < n)
          steps[k++] = (PetscInt)idx;
      }
      n = k;
      ierr = PetscSortInt(n, steps); CHKERRQ(ierr);
    }
    closedir(d);
  }

  /* MPIU_INT, not MPI_INT: PetscInt is 64-bit under --with-64-bit-indices, and
   * the original's MPI_INT would have silently corrupted the list there. */
  ierr = MPI_Bcast(&n, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  if (rank != 0 && n > 0) { ierr = PetscMalloc1(n, &steps); CHKERRQ(ierr); }
  if (n > 0) { ierr = MPI_Bcast(steps, n, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr); }

  *steps_out = steps;
  *n_out     = n;
  PetscFunctionReturn(0);
}

/* Recover simulated time per step from the producing run's SSA_evo.dat, whose
 * columns are: sub_interf/eps, tot_ice, t, step, dt, tot_air, tot_rhov,
 * tot_mass (monitoring.c). Without this the CSV's `time` column -- the whole
 * reason the schema exists -- would be unavailable on the replay path. */
static PetscErrorCode KeffLoadTimeMap(const char *path, PetscInt **steps_out,
                                      PetscReal **times_out, PetscInt *n_out)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscInt       n = 0, *steps = NULL;
  PetscReal     *times = NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

  if (rank == 0) {
    FILE *fp = fopen(path, "r");
    if (fp) {
      char line[512];
      while (fgets(line, sizeof(line), fp)) {
        double c[8];
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
                   &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7]) == 8) n++;
      }
      if (n > 0) {
        PetscInt k = 0;
        rewind(fp);
        ierr = PetscMalloc1(n, &steps); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &times); CHKERRQ(ierr);
        while (fgets(line, sizeof(line), fp) && k < n) {
          double c[8];
          if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
                     &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7]) == 8) {
            steps[k] = (PetscInt)(c[3] + 0.5);
            times[k] = (PetscReal)c[2];
            k++;
          }
        }
        n = k;
      }
      fclose(fp);
    }
  }

  ierr = MPI_Bcast(&n, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  if (n > 0) {
    if (rank != 0) {
      ierr = PetscMalloc1(n, &steps); CHKERRQ(ierr);
      ierr = PetscMalloc1(n, &times); CHKERRQ(ierr);
    }
    ierr = MPI_Bcast(steps, n, MPIU_INT,  0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(times, n, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
  }

  *steps_out = steps;
  *times_out = times;
  *n_out     = n;
  PetscFunctionReturn(0);
}

PetscErrorCode KeffReplay(AppCtx *app, Vec U)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;
  PetscInt      *steps = NULL, nsteps = 0;
  PetscInt      *tsteps = NULL, ntimes = 0;
  PetscReal     *times = NULL;
  PetscBool      warned_no_times = PETSC_FALSE;

  PetscFunctionBegin;

  ierr = KeffScanSnapshots(kc->replay_dir, &steps, &nsteps); CHKERRQ(ierr);
  if (nsteps == 0)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "-keff_replay: no sol_*.dat snapshots found in '%s'", kc->replay_dir);

  ierr = KeffLoadTimeMap(kc->replay_times, &tsteps, &times, &ntimes); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  [keff] replay: %d snapshot(s) in %s, stride %d%s\n",
    (int)nsteps, kc->replay_dir, (int)kc->replay_stride,
    (ntimes > 0) ? "" : "   (no step->time map; `time` will be -1)"); CHKERRQ(ierr);

  for (PetscInt i = 0; i < nsteps; i += kc->replay_stride) {
    char      path[PETSC_MAX_PATH_LEN];
    PetscReal t = -1.0;

    ierr = PetscSNPrintf(path, sizeof(path), "%s/sol_%05d.dat",
                         kc->replay_dir, (int)steps[i]); CHKERRQ(ierr);

    for (PetscInt j = 0; j < ntimes; j++)
      if (tsteps[j] == steps[i]) { t = times[j]; break; }

    if (t < 0.0 && ntimes == 0 && !warned_no_times) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "  [keff] WARNING: %s not readable; the CSV's time column will be -1. "
        "Supply one with -keff_replay_times.\n", kc->replay_times); CHKERRQ(ierr);
      warned_no_times = PETSC_TRUE;
    }

    /* The producer's own IGA, so no length check or mesh reconstruction. */
    ierr = IGAReadVec(app->iga, U, path); CHKERRQ(ierr);
    ierr = KeffSample(app, steps[i], t, U); CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  [keff] replay complete: %d sample(s) written to %s\n",
    (int)kc->nsamples, kc->csv_path); CHKERRQ(ierr);

  ierr = PetscFree(steps);  CHKERRQ(ierr);
  ierr = PetscFree(tsteps); CHKERRQ(ierr);
  ierr = PetscFree(times);  CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
