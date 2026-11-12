#include "monitoring.h"
#include "assembly.h"
#include "material_properties.h"
#include <petscstring.h>
#include <petsc/private/tsimpl.h>   /* for direct access to ts->dtmin */

/* Open a per-step log file. Normally the file is created (truncated) and the
 * caller writes its header. Under -restart_append the file is opened for
 * append and *append is set, so the caller skips the header: the first leg's
 * header and rows (trimmed to step <= N by the job script) stay in place. */
static PetscErrorCode OpenLogViewer(AppCtx *user, const char *fn, PetscViewer *v,
                                    PetscBool *append)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  *append = user->restart_append;
  if (!*append) {
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, fn, v); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = PetscViewerCreate(PETSC_COMM_WORLD, v); CHKERRQ(ierr);
  ierr = PetscViewerSetType(*v, PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(*v, FILE_MODE_APPEND); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(*v, fn); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* True on the monitor call for the state loaded by -restart_append: that
 * step's rows are already in the logs and its sol_N.dat is on disk. */
static PetscBool RestartLoadedStep(const AppCtx *user, PetscInt step)
{
  return (PetscBool)(user->restart_append && step == user->restart_step);
}

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  AppCtx *user = (AppCtx *)mctx;

  /* InterfaceCFLMonitor runs first and may already have condemned this step;
   * never record a state BoundsRollbackPreStep is about to discard. */
  if (user->bounds_violated) PetscFunctionReturn(0);

  //-------- compute beta_sub
  Vec localU;
  const PetscScalar *arrayU;
  IGAElement element;
  IGAPoint point;
  PetscScalar *UU;
  PetscScalar rhovs, arg_kin, sigm0, sigm_surf, v_kin, alp;
  PetscInt indd=0;
  PetscReal a1=5.0, a2=0.1581, bet_max=0.0, bet_min=1.0e30;
  PetscReal bet0, d0, rho_rhovs, d0_sub,  beta_sub, lambda_sub, tau_sub;

  if (user->flag_Tdep) {
    ierr = IGAGetLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
    ierr = IGABeginElement(user->iga,&element);CHKERRQ(ierr);
    while (IGANextElement(user->iga,element)) {
        ierr = IGAElementGetValues(element,arrayU,&UU);CHKERRQ(ierr);
        ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
        while (IGAElementNextPoint(element,point)) {
            PetscScalar solS[3];
            ierr = IGAPointFormValue(point,UU,&solS[0]);CHKERRQ(ierr);
            RhoVS_I(user,solS[1],&rhovs,NULL);
            sigm_surf=fabs(solS[2]-rhovs)/rhovs;
            rho_rhovs = user->rho_ice/rhovs;

            arg_kin = 1.38e-23*(solS[1]+273.15)/(2.0*3.14159*3.0e-26);
            v_kin = pow(arg_kin,0.5)/rho_rhovs;

            Sigma0(solS[1],&sigm0);
            if(sigm0<=0.0) PetscPrintf(PETSC_COMM_SELF,"ERROR: Negative Sigma0 value %e\n",sigm0);
            if(sigm_surf < sigm0/69.0775) alp = 1.0e-30;
            else alp = exp(-sigm0/sigm_surf);

            if(alp*v_kin<1.0e-30) bet0 = 1.0e30;
            else bet0 = 1.0/(alp*v_kin);
            /* Capillary length: the single value computed at start-up from K&P
             * Eq. 13 at temp0 (lunar_main.c), used for the whole run. */
            d0 = user->d0_sub0;

            if(bet0>bet_max) bet_max = bet0;
            if(bet0<bet_min) bet_min = bet0;
            d0_sub   = d0/rho_rhovs;  
            beta_sub = bet0/rho_rhovs;
            lambda_sub = a1*user->eps/d0_sub;
            tau_sub    = user->eps*lambda_sub*(beta_sub/a1 + a2*user->eps/user->diff_sub + a2*user->eps/user->dif_vap);

            user->mob[indd]  = user->eps/3.0/tau_sub;
            user->alph[indd] = lambda_sub/tau_sub;

            indd ++;
        }
        ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga,&element);CHKERRQ(ierr);
    ierr = IGARestoreLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
    PetscReal B_min, B_max;
    ierr = MPI_Allreduce(&bet_max,&B_max,1,MPI_DOUBLE,MPI_MAX,PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&bet_min,&B_min,1,MPI_DOUBLE,MPI_MIN,PETSC_COMM_WORLD);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD," b_min %.2e b_max %.2e\n",B_min,B_max);

    // After computing beta_sub, we set the flag to 0...
    user->flag_Tdep = PETSC_FALSE;

    // Print the new mobility and alpha values for verification
    PetscPrintf(PETSC_COMM_WORLD, "M0_sub new: %.6e\n", user->mob_sub);
    PetscPrintf(PETSC_COMM_WORLD, "alpha_sub new: %.6e\n", user->alph_sub);
  }


  //-------- domain integrals
  #define N_INTEGRALS 5
  PetscScalar stats[N_INTEGRALS] = {0.0,0.0,0.0,0.0,0.0};
  ierr = IGAComputeScalar(user->iga,U,N_INTEGRALS,&stats[0],Integration,mctx);CHKERRQ(ierr);
  PetscReal tot_ice    = PetscRealPart(stats[0]);
  PetscReal sub_interf = PetscRealPart(stats[1]);
  PetscReal tot_air    = PetscRealPart(stats[2]);
  PetscReal tot_temp   = PetscRealPart(stats[3]);
  PetscReal tot_rhov   = PetscRealPart(stats[4]);

  /* Total system mass: only counts vapor in the air phase (TOT_RHOV is already
   * weighted by phi_a). This is the conservation quantity that should stay flat. */
  PetscReal tot_mass = user->rho_ice * tot_ice + tot_rhov;

  /* Store initial integrals at step 0 for later percentage reporting. */
  if (step == 0) {
    user->tot_ice_0  = tot_ice;
    user->tot_air_0  = tot_air;
    user->tot_rhov_0 = tot_rhov;
    user->tot_mass_0 = tot_mass;
  }

  //-------- phase-field min/max bounds (printed every step for out-of-bounds detection)
  {
    PetscReal phi_ice_min =  1.0e30, phi_ice_max = -1.0e30;
    PetscReal phi_air_min =  1.0e30, phi_air_max = -1.0e30;

    Vec localUb; const PetscScalar *arrayUb;
    IGAElement elementb; IGAPoint pointb; PetscScalar *UUb;
    ierr = IGAGetLocalVecArray(user->iga, U, &localUb, &arrayUb); CHKERRQ(ierr);
    ierr = IGABeginElement(user->iga, &elementb); CHKERRQ(ierr);
    while (IGANextElement(user->iga, elementb)) {
      ierr = IGAElementGetValues(elementb, arrayUb, &UUb); CHKERRQ(ierr);
      ierr = IGAElementBeginPoint(elementb, &pointb); CHKERRQ(ierr);
      while (IGAElementNextPoint(elementb, pointb)) {
        PetscScalar solb[3];
        ierr = IGAPointFormValue(pointb, UUb, &solb[0]); CHKERRQ(ierr);
        PetscReal fi = PetscRealPart(solb[0]);
        PetscReal fa = 1.0 - fi;
        if (fi < phi_ice_min) phi_ice_min = fi;
        if (fi > phi_ice_max) phi_ice_max = fi;
        if (fa < phi_air_min) phi_air_min = fa;
        if (fa > phi_air_max) phi_air_max = fa;
      }
      ierr = IGAElementEndPoint(elementb, &pointb); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &elementb); CHKERRQ(ierr);
    ierr = IGARestoreLocalVecArray(user->iga, U, &localUb, &arrayUb); CHKERRQ(ierr);

    PetscReal Gfi_min, Gfi_max, Gfa_min, Gfa_max;
    ierr = MPI_Allreduce(&phi_ice_min, &Gfi_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_ice_max, &Gfi_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_air_min, &Gfa_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_air_max, &Gfa_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

    /* Scientific format, not %.4f. This line is the per-step record of how far
     * the phase field strays outside [0,1], and postprocess/plot_phi_bounds.py
     * reads it; at fixed precision a -4e-5 excursion printed as "-0.0000",
     * indistinguishable from a healthy -1e-23. The magnitude IS the
     * diagnostic. */
    PetscPrintf(PETSC_COMM_WORLD,
        "  BOUNDS: phi_ice [%.4e, %.4e]  phi_air [%.4e, %.4e]\n",
        Gfi_min, Gfi_max, Gfa_min, Gfa_max);

    /* If any phase field has left [phase_lo, phase_hi], roll the just-finished
     * step back, halve dt, and let TS retry. Abort only if dt would fall below
     * dtmin, or if the violation is at step 0 (i.e., a bad initial condition). */
    PetscBool oob = (Gfi_min < user->phase_lo || Gfi_max > user->phase_hi ||
                     Gfa_min < user->phase_lo || Gfa_max > user->phase_hi);
    if (oob) {
      PetscReal cur_dt;
      ierr = TSGetTimeStep(ts, &cur_dt); CHKERRQ(ierr);
      PetscReal new_dt = cur_dt / 2.0;
      PetscBool can_retry = (step > 0) && (new_dt >= ts->dtmin);

      if (!can_retry) {
        PetscPrintf(PETSC_COMM_WORLD,
            "\033[31m[ABORT] Phase field out of bounds [%.2f, %.2f] at step %d\n"
            "  phi_ice [%.4f, %.4f]  phi_air [%.4f, %.4f]\n"
            "  Cannot retry: %s\033[0m\n",
            user->phase_lo, user->phase_hi, step,
            Gfi_min, Gfi_max, Gfa_min, Gfa_max,
            (step == 0) ? "bad initial condition (step 0)"
                        : "dt already at dtmin");
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_NOT_CONVERGED,
                "Phase field out of bounds at step %" PetscInt_FMT
                " — phi_ice [%.4g, %.4g]  phi_air [%.4g, %.4g]"
                "  (bounds [%.2g, %.2g])",
                step,
                Gfi_min, Gfi_max, Gfa_min, Gfa_max,
                user->phase_lo, user->phase_hi);
      }

      PetscPrintf(PETSC_COMM_WORLD,
          "\033[33m[WARN] Phase field out of bounds at step %d — "
          "deferring rollback to next pre-step, dt %.3e -> %.3e\n"
          "  phi_ice [%.4f, %.4f]  phi_air [%.4f, %.4f]\033[0m\n",
          step, cur_dt, new_dt,
          Gfi_min, Gfi_max, Gfa_min, Gfa_max);

      /* Can't call TSRollBack here — ts->vec_sol is read-locked inside
       * TSMonitor. Defer to the BoundsRollbackPreStep callback, which
       * fires before the next TSStep (when the vector is writable). */
      user->bounds_violated = PETSC_TRUE;
      user->bounds_new_dt   = new_dt;

      /* Skip the rest of the monitor for this step (it's being undone). */
      PetscFunctionReturn(0);
    }
  }

  /* ---- Stall detector -------------------------------------------------
   * A run that cannot advance is otherwise INVISIBLE: TS keeps stepping, the
   * clock reaches t_final, and the job reports "Solution completed" with a
   * solution that stopped changing thousands of steps earlier. Both eps =
   * 0.75 um runs of batch 2026-09-12 died this way; the only symptom was that
   * 52 of 61 sol_*.dat files were byte-identical, which nobody checks by
   * default.
   *
   * ||U|| not moving AT ALL -- bit-identical, not merely small -- for
   * stall_limit consecutive steps is not a physical steady state in a
   * relaxation problem; it means the solve is not updating the solution. Abort
   * and say so rather than burn the rest of the wall clock.
   *
   * Disable with -stall_limit 0. */
  if (user->stall_limit > 0) {
    PetscReal unorm;
    ierr = VecNorm(U, NORM_2, &unorm); CHKERRQ(ierr);
    if (step > 0 && unorm == user->stall_norm) {
      user->stall_count++;
      if (user->stall_count >= user->stall_limit)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_NOT_CONVERGED,
                "Solution has not changed for %" PetscInt_FMT " consecutive steps "
                "(||U|| = %.17g, bit-identical) at step %" PetscInt_FMT ", t = %g s.\n"
                "  The time stepper is advancing but the solve is not updating the "
                "solution, so every remaining step is wasted and the run would have\n"
                "  reported SUCCESS at t_final with a frozen field. Usual cause: the "
                "phase field left [phase_lo, phase_hi], after which the residual's\n"
                "  domain guard returns zero and SNES reads ||F|| = 0 as converged. "
                "Check phi_min/phi_max, and -dtmax against tau_sub.\n"
                "  Disable this check with -stall_limit 0 if a frozen solution is "
                "genuinely expected.",
                user->stall_count, (double)unorm, step, (double)t);
    } else {
      user->stall_count = 0;
    }
    user->stall_norm = unorm;
  }

  //-------------
  PetscReal dt;
  TSGetTimeStep(ts,&dt);

  //------ printf information (robust table header + aligned columns)
  {
    const PetscInt headerEvery = 1;  // print header every step (can change to 10 later)
    PetscErrorCode ierr2;

    if ((step % headerEvery) == 0) {
      // Build header line using the SAME fixed-width fields as the data row
      char header[512];
      ierr2 = PetscSNPrintf(header, sizeof(header),
                            " %5s | %12s | %9s | %10s | %10s | %9s | %9s | %10s | %10s",
                            "STEP",
                            "TIME [s]",
                            "DT [s]",
                            "TOT_ICE",
                            "TOT_AIR",
                            "TEMP",
                            "TOT_RHOV",
                            "I-A INTERF",
                            "TOTAL_MASS");
      CHKERRQ(ierr2);

      // Separator line: exactly matches header length
      size_t hlen = 0;
      ierr2 = PetscStrlen(header, &hlen);
      CHKERRQ(ierr2);

      char sep[512];
      size_t i;
      for (i = 0; i < hlen && i < sizeof(sep) - 1; i++) sep[i] = '-';
      sep[i] = '\0';

      PetscPrintf(PETSC_COMM_WORLD, "  ===============================================================================\n");
      PetscPrintf(PETSC_COMM_WORLD, "  >>> DOMAIN INTEGRALS (PER TIME STEP)\n");
      PetscPrintf(PETSC_COMM_WORLD, "  ===============================================================================\n");
      PetscPrintf(PETSC_COMM_WORLD, "  %s\n", sep);
      PetscPrintf(PETSC_COMM_WORLD, "  %s\n", header);
      PetscPrintf(PETSC_COMM_WORLD, "  %s\n", sep);
    }

    // Data row: uses matching widths so it lines up under the header
    PetscPrintf(PETSC_COMM_WORLD,
                " %5d | %12.5e | %9.3e | %10.3e | %10.3e | %9.3e | %9.3e | %10.3e | %10.3e\n",
                step, t, dt,
                tot_ice, tot_air, tot_temp, tot_rhov,
                sub_interf, tot_mass);

    /* Percentage-change row (relative to initial values at step 0). */
    if (step > 0 && user->tot_ice_0 > 0.0) {
      PetscReal pct_ice  = (tot_ice  - user->tot_ice_0)  / user->tot_ice_0  * 100.0;
      PetscReal pct_air  = (tot_air  - user->tot_air_0)  / user->tot_air_0  * 100.0;
      PetscReal pct_rhov = (user->tot_rhov_0 > 0.0)
                           ? (tot_rhov - user->tot_rhov_0) / user->tot_rhov_0 * 100.0 : 0.0;
      PetscReal pct_mass = (user->tot_mass_0 > 0.0)
                           ? (tot_mass - user->tot_mass_0) / user->tot_mass_0 * 100.0 : 0.0;
      PetscPrintf(PETSC_COMM_WORLD,
                  "       |              |           | %+9.3f%% | %+9.3f%% |           | %+8.3f%% |            | %+9.3f%%\n",
                  pct_ice, pct_air, pct_rhov, pct_mass);
    }

    // Optional: add a blank line every N rows for readability (set to 0 to disable)
    // if (step > 0 && (step % 25) == 0) PetscPrintf(PETSC_COMM_WORLD, "\n");
  }

  /* SSA_evo.dat gets a row at every time step, whatever -outp and -t_interv say:
   * the row is cheap and we want the full record (the restart path also reads
   * it). -outp / -t_interv govern only the sol_*.dat snapshots (OutputMonitor). */
  PetscInt print = 1;

  /* solver_evo.dat: Newton and Krylov iterations spent on the step just taken.
   * TS keeps running totals (rejected attempts included), so the per-step cost
   * is the difference since the last call. Three integer reads and one short
   * line per step -- no extra solves, norms or reductions.
   *   step  t  dt  newton_its  krylov_its  krylov_per_newton  rejections
   * krylov_per_newton is the conditioning diagnostic: how hard each linear
   * solve was. rejections counts step attempts thrown away and retried. */
  {
    PetscInt snes_tot = 0, ksp_tot = 0, rej_tot = 0;
    ierr = TSGetSNESIterations(ts, &snes_tot);CHKERRQ(ierr);
    ierr = TSGetKSPIterations(ts, &ksp_tot);CHKERRQ(ierr);
    ierr = TSGetStepRejections(ts, &rej_tot);CHKERRQ(ierr);
    PetscInt d_snes = snes_tot - user->solver_snes_prev;
    PetscInt d_ksp  = ksp_tot  - user->solver_ksp_prev;
    PetscInt d_rej  = rej_tot  - user->solver_rej_prev;
    user->solver_snes_prev = snes_tot;
    user->solver_ksp_prev  = ksp_tot;
    user->solver_rej_prev  = rej_tot;

    if (user->solver_view == NULL) {
      char fname[256];
      const char *sdir = getenv("folder");
      if (sdir) { sprintf(fname,"%s/solver_evo.dat",sdir); }
      else      { sprintf(fname,"solver_evo.dat"); }
      PetscBool app;
      ierr = OpenLogViewer(user, fname, &user->solver_view, &app);CHKERRQ(ierr);
      if (!app) {
        ierr = PetscViewerASCIIPrintf(user->solver_view,
                 "# step t dt newton_its krylov_its krylov_per_newton rejections\n");CHKERRQ(ierr);
      }
    }
    if (!RestartLoadedStep(user, step)) {
    ierr = PetscViewerASCIIPrintf(user->solver_view,"%d %e %e %d %d %.2f %d\n",
                                  (int)step, (double)t, (double)dt, (int)d_snes, (int)d_ksp,
                                  d_snes > 0 ? (double)d_ksp/(double)d_snes : 0.0,
                                  (int)d_rej);CHKERRQ(ierr);
    ierr = PetscViewerFlush(user->solver_view);CHKERRQ(ierr);
    }
  }

  if(print==1) {
    /* Open the SSA_evo.dat viewer ONCE (first call), then reuse it and flush
     * every step. Re-opening a viewer per step previously exhausted file
     * descriptors around step 15 and, with no CHKERRQ on any of the calls,
     * failed silently -- truncating SSA_evo.dat while outp.txt kept growing.
     * Flushing each step keeps the on-disk file current for mid-run reads. */
    if (user->ssa_view == NULL) {
      char filedata[256];
      const char *dir = getenv("folder");
      if (dir) { sprintf(filedata,"%s/SSA_evo.dat",dir); }
      else     { sprintf(filedata,"SSA_evo.dat"); }
      PetscBool app;   /* SSA_evo.dat has no header either way */
      ierr = OpenLogViewer(user, filedata, &user->ssa_view, &app);CHKERRQ(ierr);
    }
    if (!RestartLoadedStep(user, step)) {
      ierr = PetscViewerASCIIPrintf(user->ssa_view,"%e %e %e %d %e %e %e %e\n",
                                    sub_interf/user->eps, tot_ice, t, step, dt,
                                    tot_air, tot_rhov, tot_mass);CHKERRQ(ierr);
      ierr = PetscViewerFlush(user->ssa_view);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OutputMonitor(TS ts, PetscInt step, PetscReal t, Vec U, 
                              void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  /* InterfaceCFLMonitor runs first and may already have condemned this step;
   * never record a state BoundsRollbackPreStep is about to discard. */
  if (user->bounds_violated) PetscFunctionReturn(0);

  // Check if it's the first step
  if (step == 0) {
    const char *env = "folder";
    char *dir;
    dir = getenv(env);

    char fileiga[256];
    sprintf(fileiga, "%s/igasol.dat", dir);

    ierr = IGAWrite(user->iga, fileiga);CHKERRQ(ierr);
  }

  /* Always write step 0 (the raw initial condition) and step 1 (the first
   * computed step). Step 0 already falls out of both schedules below, but
   * step 1 generally does not: under time-uniform output the next scheduled
   * snapshot is a whole t_interv away, so the first *solved* state is
   * otherwise never on disk. That state is what movies and diagnostics
   * should start from -- the IC's vapor field is a uniform hum0*rho_vs, so
   * its supersaturation is a constant that carries no dynamics and skews
   * any percentile range computed across the series. Forcing step 1 gives
   * an "initial frame" that is a real solution, so the IC can be skipped. */
  PetscInt print = (step <= 1) ? 1 : 0;

  /* -restart_append: sol_N.dat is the file just loaded. Do not rewrite it;
   * only move the next output time past t, on the same t_interv grid as the
   * first leg (t_out starts at 0 and steps by t_interv either way). */
  if (RestartLoadedStep(user, step)) {
    if (user->outp <= 0) { while (user->t_out <= t) user->t_out += user->t_interv; }
    PetscFunctionReturn(0);
  }

  // Check if it's time to print output
  if (!print) {
    if (user->outp > 0) {   // Print output every user->outp steps
      if (step % user->outp == 0) print = 1;
    }
    else {                  // Print output every user->t_interv seconds
      if (t >= user->t_out) print = 1;
    }
  }

  // If it's time to print output, do the following
  if (print == 1) {
    PetscPrintf(PETSC_COMM_WORLD, "OUTPUT print!\n");
    /* Advance t_out past the current time: a single += degenerates to
     * every-step output whenever dt > t_interv (t_out falls permanently
     * behind t), defeating the 1000-snapshot output cap. */
    while (user->t_out <= t) user->t_out += user->t_interv;

    // Get the directory path from the environment variable
    const char *env = "folder";
    char *dir;
    dir = getenv(env);

    // Create the filename for the output file
    char filename[256];
    sprintf(filename, "%s/sol_%05d.dat", dir, step);

    // Write the vector U to the output file
    ierr = IGAWriteVec(user->iga, U, filename);
    CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


/* =========================================================================
 * InterfaceCFLMonitor
 *
 * Interface-CFL timestep limiter. After every ACCEPTED step, measure the
 * fastest pointwise phase change rate of that step,
 *
 *     rate = ||phi^n - phi^{n-1}||_inf / (t^n - t^{n-1}),
 *
 * and clamp the controller's proposed next dt so that no DOF can change by
 * more than user->cfl_dphimax (default 0.2) in one step:
 *
 *     dt_next <= cfl_dphimax / rate.
 *
 * Rationale (2026-07-09 two-grain campaign): during fast interface events —
 * a grain's terminal collapse, initial neck formation — the front velocity
 * diverges (v = d0/(beta*r) as r -> 0) and a quadratic B-spline front moved
 * many cells per step develops Gibbs undershoot ripples that pass both the
 * phase guard (field dips between quadrature points) and the atol test
 * (O(0.1) phi errors have ~1e-12 residuals in this dimensional scaling).
 * A static dtmax cannot track a diverging velocity: the ladder showed 4e2
 * clean but 8e2 rippling at grain-vanish. This limiter is the dynamic
 * version: NRmin/NRmax and dtmax still propose dt; the measured phase rate
 * disposes. Quiet phases are untouched (rate ~ 0 -> cap huge); events
 * throttle automatically and dt regrows the moment they end.
 *
 * ENFORCEMENT (model-freeze-v2, 2026-11-05; port of enceladus_DSM
 * caa209e0). Clamping only the NEXT dt is advisory: the step that overshot
 * has already happened and stays in the record (stage-0b B1 kept 264 steps
 * with |dphi| > 0.2, the largest 2.04). Now a step whose MEASURED |dphi|
 * exceeds cfl_dphimax is rolled back through the deferred mechanism the
 * bounds guard uses, and retried at
 *
 *     dt_retry = CFL_RETRY_SAFETY * cfl_dphimax * dt_used / dphi_max,
 *
 * so max |dphi| per step <= cfl_dphimax on every step the run keeps. The
 * forward clamp is kept (at CFL_FORWARD_MARGIN) because it prevents most
 * rollbacks. Registered FIRST among the monitors, and every other monitor
 * returns early when bounds_violated is set, so a rejected step leaves no
 * row and no snapshot. BoundsRollbackPreStep resets cfl_U_prev to the state
 * it rolls back to, so the retried step is measured against the right
 * reference whichever guard asked for the rollback.
 *
 * Runs as a TS monitor (after a converged step).
 * ========================================================================= */

/* Fractions of the dt that would hit cfl_dphimax exactly.
 * CFL_RETRY_SAFETY   after a rollback: well inside, so the retry does not
 *                    violate again.
 * CFL_FORWARD_MARGIN on the predicted next dt: the prediction uses the last
 *                    step's rate, so an accelerating front overshoots it;
 *                    targeting the limit exactly turns every acceleration
 *                    into a rollback (enceladus smoke test: 1 step in 3). */
#define CFL_RETRY_SAFETY   0.8
#define CFL_FORWARD_MARGIN 0.9
PetscErrorCode InterfaceCFLMonitor(TS ts, PetscInt step, PetscReal t, Vec U,
                                   void *mctx)
{
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  PetscFunctionBegin;

  /* Snap-to-dtmax: the NRmin growth heuristic refuses to grow dt when
   * dt*factor would exceed dtmax, stranding dt just below the cap forever
   * (e.g. 790 s against dtmax = 1000: "Max dt! can not Increase"). When
   * growth is blocked ONLY by the cap, take dtmax itself. Runs before the
   * CFL cap below, which can still pull the result back down. */
  if (ts->adap && ts->dtmax > 0.0) {
    PetscReal dt_now;
    PetscErrorCode ierr2 = TSGetTimeStep(ts, &dt_now);CHKERRQ(ierr2);
    if (dt_now < ts->dtmax && dt_now * ts->factor > ts->dtmax) {
      ierr2 = TSSetTimeStep(ts, ts->dtmax);CHKERRQ(ierr2);
    }
  }

  if (!user->flag_dtCFL) PetscFunctionReturn(0);

  if (!user->cfl_U_prev) {                    /* first call: just record */
    ierr = VecDuplicate(U, &user->cfl_U_prev); CHKERRQ(ierr);
    ierr = VecDuplicate(U, &user->cfl_diff);   CHKERRQ(ierr);
    ierr = VecCopy(U, user->cfl_U_prev);       CHKERRQ(ierr);
    user->cfl_t_prev = t;
    PetscFunctionReturn(0);
  }

  PetscReal dt_used = t - user->cfl_t_prev;
  if (dt_used > 0.0) {
    /* ||dphi||_inf of the last step, into a SCRATCH vector: the rollback
     * path returns early, so cfl_U_prev must never be used as the buffer. */
    PetscReal dphi_max;
    ierr = VecWAXPY(user->cfl_diff, -1.0, user->cfl_U_prev, U); CHKERRQ(ierr);
    ierr = VecStrideNorm(user->cfl_diff, 0, NORM_INFINITY, &dphi_max); CHKERRQ(ierr);

    if (dphi_max > 0.0) {
      /* --- hard enforcement: undo a step that already violated it ------- */
      if (dphi_max > user->cfl_dphimax) {
        PetscReal dt_retry = CFL_RETRY_SAFETY * user->cfl_dphimax * dt_used / dphi_max;
        PetscBool can_retry = (step > 0) && (dt_retry >= ts->dtmin);
        if (can_retry) {
          user->cfl_n_reject++;
          PetscPrintf(PETSC_COMM_WORLD,
              "\033[33m[WARN] Interface-CFL violated at step %d: max |dphi| "
              "%.3e > %.3e, rolling back, dt %.3e -> %.3e%s\033[0m\n",
              (int)step, (double)dphi_max, (double)user->cfl_dphimax,
              (double)dt_used, (double)dt_retry,
              (user->cfl_n_reject > 3) ? "  [3+ consecutive: check -dtCFL_dphimax]" : "");
          user->bounds_violated = PETSC_TRUE;
          user->bounds_new_dt   = dt_retry;
          /* Return WITHOUT recording U: cfl_U_prev must keep the last
           * ACCEPTED state. */
          PetscFunctionReturn(0);
        }
        PetscPrintf(PETSC_COMM_WORLD,
            "\033[33m[WARN] Interface-CFL violated at step %d (max |dphi| "
            "%.3e > %.3e) but cannot retry: %s. Accepting the step.\033[0m\n",
            (int)step, (double)dphi_max, (double)user->cfl_dphimax,
            (step == 0) ? "step 0" : "dt would fall below dtmin");
      }
      user->cfl_n_reject = 0;

      /* --- forward clamp: keep the NEXT step inside the criterion too --- */
      PetscReal dt_cap = CFL_FORWARD_MARGIN * user->cfl_dphimax * dt_used / dphi_max;
      PetscReal dt_next;
      ierr = TSGetTimeStep(ts, &dt_next); CHKERRQ(ierr);
      if (dt_cap < dt_next) {
        ierr = TSSetTimeStep(ts, dt_cap); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,
                    "Interface-CFL cap: dt %.3e -> %.3e (max |dphi|/step %.3e)\n",
                    (double)dt_next, (double)dt_cap, (double)dphi_max);
      }
    }
  }

  ierr = VecCopy(U, user->cfl_U_prev); CHKERRQ(ierr);
  user->cfl_t_prev = t;
  PetscFunctionReturn(0);
}


/* =========================================================================
 * BoundsRollbackPreStep
 *
 * TS pre-step callback. If Monitor() flagged a bounds violation on the last
 * step, undo that step and shrink dt here — where ts->vec_sol is writable
 * (it's read-locked during TSMonitor, so TSRollBack cannot be called there).
 *
 * Registered via TSSetPreStep() in src/<project>_main.c. The AppCtx pointer is
 * retrieved via TSGetApplicationContext, which is set in main() with
 * TSSetApplicationContext.
 * ========================================================================= */
PetscErrorCode BoundsRollbackPreStep(TS ts)
{
  PetscErrorCode ierr;
  AppCtx *user;

  PetscFunctionBegin;
  ierr = TSGetApplicationContext(ts, &user); CHKERRQ(ierr);

  if (user && user->bounds_violated) {
    ierr = TSRollBack(ts); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, user->bounds_new_dt); CHKERRQ(ierr);
    /* Re-anchor the step-to-step diagnostics on the state just restored.
     * Whichever guard asked for the rollback (bounds in Monitor, or the CFL
     * check), InterfaceCFLMonitor or ResolutionMonitor may already have
     * recorded the rejected state as "previous"; the retried step must be
     * measured against the state it actually starts from. */
    {
      Vec       X;
      PetscReal t0;
      ierr = TSGetSolution(ts, &X); CHKERRQ(ierr);
      ierr = TSGetTime(ts, &t0); CHKERRQ(ierr);
      if (user->cfl_U_prev) { ierr = VecCopy(X, user->cfl_U_prev); CHKERRQ(ierr); user->cfl_t_prev = t0; }
      if (user->res_U_prev) { ierr = VecCopy(X, user->res_U_prev); CHKERRQ(ierr); user->res_t_prev = t0; }
    }
    user->bounds_violated = PETSC_FALSE;
    user->bounds_new_dt   = 0.0;
  }

  PetscFunctionReturn(0);
}

/* Physical flux out through the Dirichlet faces, one row per step in bflux.dat:
 *   vapor:  J = -int D_v(T) phi_aef grad(rho_v).n dGamma   [kg/s]
 *   heat:   Q = -int k(phi) grad(T).n dGamma               [W]
 * per unit depth in 2D, and the full 2*pi*r surface in axisymmetric mode.
 * Positive is out of the domain. D_v and k are the UNSCALED coefficients:
 * the residual multiplies them by xi_v and xi_T to accelerate the run, but
 * the physical flux uses the true values (Kaempfer and Plapp 2009, Sec. II D).
 * phi_aef = max(1 - phi, air_lim) and the clamped phi in k(phi) match
 * Residual_A1. Read-only: it does not change the solution. */
PetscErrorCode BoundaryFluxMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  AppCtx  *user = (AppCtx *)mctx;
  PetscInt dim  = user->dim;
  /* InterfaceCFLMonitor runs first and may already have condemned this step;
   * never record a state BoundsRollbackPreStep is about to discard. */
  if (user->bounds_violated) PetscFunctionReturn(0);
  PetscReal flux[3][2][2];   /* [axis][side][0 = vapor, 1 = heat] */
  PetscReal gflux[3][2][2];
  ierr = PetscMemzero(flux, sizeof(flux)); CHKERRQ(ierr);

  Vec localU; const PetscScalar *arrayU;
  IGAElement element; IGAPoint point; PetscScalar *UU;
  ierr = IGAGetLocalVecArray(user->iga, U, &localU, &arrayU); CHKERRQ(ierr);
  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
    ierr = IGAElementGetValues(element, arrayU, &UU); CHKERRQ(ierr);
    while (IGAElementNextForm(element, user->flux_face)) {
      if (!element->atboundary) continue;
      PetscInt axis = element->boundary_id / 2, side = element->boundary_id % 2;
      ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
      while (IGAElementNextPoint(element, point)) {
        /* FormGrad packs [dof][dim] with stride dim, so index g[dof*dim + l]. */
        PetscScalar sol[3], g[9];
        ierr = IGAPointFormValue(point, UU, &sol[0]); CHKERRQ(ierr);
        ierr = IGAPointFormGrad (point, UU, &g[0]); CHKERRQ(ierr);
        PetscReal phi_c = PetscMin(PetscMax(PetscRealPart(sol[0]), 0.0), 1.0);
        PetscReal phi_aef = PetscMax(1.0 - phi_c, user->air_lim);
        PetscScalar thcond, dif_vap;
        ThermalCond(user, phi_c, &thcond, NULL);
        VaporDiffus(user, sol[1], &dif_vap, NULL);
        PetscReal dT_dn = 0.0, drhov_dn = 0.0;
        for (PetscInt l = 0; l < dim; l++) {
          dT_dn    += PetscRealPart(g[1 * dim + l]) * point->normal[l];
          drhov_dn += PetscRealPart(g[2 * dim + l]) * point->normal[l];
        }
        PetscReal JW = point->detJac[0] * point->weight[0];
        if (user->axisym) {
          PetscReal xphys[3] = {0.0, 0.0, 0.0};
          IGAPointFormPoint(point, xphys);
          JW *= 2.0 * PETSC_PI * xphys[1];
        }
        flux[axis][side][0] += -PetscRealPart(dif_vap) * phi_aef * drhov_dn * JW;
        flux[axis][side][1] += -PetscRealPart(thcond) * dT_dn * JW;
      }
      ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(user->iga, U, &localU, &arrayU); CHKERRQ(ierr);
  ierr = MPI_Allreduce(flux, gflux, 12, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

  static const char axes[3] = {'x', 'y', 'z'};
  if (user->bflux_view == NULL) {
    char filedata[256];
    const char *dir = getenv("folder");
    if (dir) { snprintf(filedata, sizeof(filedata), "%s/bflux.dat", dir); }
    else     { snprintf(filedata, sizeof(filedata), "bflux.dat"); }
    PetscBool app;
    ierr = OpenLogViewer(user, filedata, &user->bflux_view, &app); CHKERRQ(ierr);
    if (!app) {
    ierr = PetscViewerASCIIPrintf(user->bflux_view,
             "# outward physical flux (unscaled D_v, k); vapor kg/s, heat W; %s\n# t step",
             user->axisym ? "axisymmetric, full 2*pi*r surface" : "per unit depth in 2D");
    CHKERRQ(ierr);
    for (PetscInt l = 0; l < dim; l++)
      for (PetscInt m = 0; m < 2; m++)
        if (user->flux_face[l][m]) {
          ierr = PetscViewerASCIIPrintf(user->bflux_view, " J_vap_%c%d Q_heat_%c%d",
                                        axes[l], (int)m, axes[l], (int)m); CHKERRQ(ierr);
        }
    ierr = PetscViewerASCIIPrintf(user->bflux_view, "\n"); CHKERRQ(ierr);
    }
  }
  if (RestartLoadedStep(user, step)) PetscFunctionReturn(0);
  ierr = PetscViewerASCIIPrintf(user->bflux_view, "%.10e %d", (double)t, (int)step); CHKERRQ(ierr);
  for (PetscInt l = 0; l < dim; l++)
    for (PetscInt m = 0; m < 2; m++)
      if (user->flux_face[l][m]) {
        ierr = PetscViewerASCIIPrintf(user->bflux_view, " %.10e %.10e",
                                      (double)gflux[l][m][0], (double)gflux[l][m][1]); CHKERRQ(ierr);
      }
  ierr = PetscViewerASCIIPrintf(user->bflux_view, "\n"); CHKERRQ(ierr);
  ierr = PetscViewerFlush(user->bflux_view); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


/* =========================================================================
 * Resolution checker
 *
 * Kaempfer and Plapp (2009, Eq. 45): the phase field relaxes on a time scale
 * that must be shorter than the time the interface needs to cross its own
 * width, which bounds eps < d0/(beta0 v_n). preprocess/comp_eps.py sizes eps
 * at HALF that bound (safety 0.5) for an expected v_n. Here v_n is MEASURED
 * every step, v_n = |dphi/dt| / |grad phi| at quadrature points inside the
 * interface band (0.1 < phi < 0.9), and the number
 *
 *     Pi = beta_sub0 * v_n * eps / d0_sub0
 *
 * is reported: Pi <= 0.5 is inside comp_eps's design margin, 0.5 < Pi <= 1
 * is outside the margin but inside the bound, Pi > 1 violates Eq. 45 (the
 * mesh/eps is too coarse for the speed the run actually reached). The same
 * pass records the largest phi on any vapor-reservoir face (res_watch_face,
 * whatever its phi condition): ice should not reach the reservoir, and with
 * the do-nothing phi condition (-phi_open_bc 1) it spreads into a film.
 * dphi/dt is a backward difference over the last accepted step. Read-only.
 * ========================================================================= */
PetscErrorCode ResolutionMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  AppCtx  *user = (AppCtx *)mctx;
  PetscInt dim  = user->dim;
  /* InterfaceCFLMonitor runs first and may already have condemned this step;
   * never record a state BoundsRollbackPreStep is about to discard. */
  if (user->bounds_violated) PetscFunctionReturn(0);

  if (!user->res_U_prev) {
    ierr = VecDuplicate(U, &user->res_U_prev); CHKERRQ(ierr);
    ierr = VecCopy(U, user->res_U_prev); CHKERRQ(ierr);
    user->res_t_prev = t;
    user->res_t0     = t;
    PetscFunctionReturn(0);
  }
  PetscReal dt = t - user->res_t_prev;
  if (dt <= 0.0) PetscFunctionReturn(0);

  PetscReal vmax = 0.0, open_max = 0.0;
  Vec lU, lP; const PetscScalar *aU, *aP;
  IGAElement el; IGAPoint pt; PetscScalar *UU, *PP;
  ierr = IGAGetLocalVecArray(user->iga, U, &lU, &aU); CHKERRQ(ierr);
  ierr = IGAGetLocalVecArray(user->iga, user->res_U_prev, &lP, &aP); CHKERRQ(ierr);
  ierr = IGABeginElement(user->iga, &el); CHKERRQ(ierr);
  while (IGANextElement(user->iga, el)) {
    ierr = IGAElementGetValues(el, aU, &UU); CHKERRQ(ierr);
    ierr = IGAElementGetValues(el, aP, &PP); CHKERRQ(ierr);
    while (IGAElementNextForm(el, user->res_watch_face)) {
      ierr = IGAElementBeginPoint(el, &pt); CHKERRQ(ierr);
      while (IGAElementNextPoint(el, pt)) {
        PetscScalar s[3];
        ierr = IGAPointFormValue(pt, UU, &s[0]); CHKERRQ(ierr);
        PetscReal phi = PetscRealPart(s[0]);
        if (el->atboundary) {                     /* reservoir-face point */
          open_max = PetscMax(open_max, phi);
          continue;
        }
        if (phi <= 0.1 || phi >= 0.9) continue;   /* interface band only */
        PetscScalar sp[3], g[9];
        ierr = IGAPointFormValue(pt, PP, &sp[0]); CHKERRQ(ierr);
        ierr = IGAPointFormGrad(pt, UU, &g[0]); CHKERRQ(ierr);
        PetscReal gn = 0.0;
        for (PetscInt l = 0; l < dim; l++) gn += PetscSqr(PetscRealPart(g[l]));
        gn = PetscSqrtReal(gn);
        if (gn <= 0.0) continue;
        PetscReal v = PetscAbsReal(phi - PetscRealPart(sp[0])) / dt / gn;
        vmax = PetscMax(vmax, v);
      }
      ierr = IGAElementEndPoint(el, &pt); CHKERRQ(ierr);
    }
  }
  ierr = IGAEndElement(user->iga, &el); CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(user->iga, U, &lU, &aU); CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(user->iga, user->res_U_prev, &lP, &aP); CHKERRQ(ierr);
  PetscReal loc[2] = {vmax, open_max}, glob[2];
  ierr = MPI_Allreduce(loc, glob, 2, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
  vmax = glob[0]; open_max = glob[1];

  PetscReal Pi = user->beta_sub0 * vmax * user->eps / user->d0_sub0;
  PetscBool startup = (PetscBool)(t - user->res_t0 < user->res_t_skip);
  user->res_open_phi_max = PetscMax(user->res_open_phi_max, open_max);
  if (open_max > 0.01) user->res_n_open++;
  if (startup) {
    user->res_pi_max_startup = PetscMax(user->res_pi_max_startup, Pi);
  } else {
    user->res_pi_max = PetscMax(user->res_pi_max, Pi);
    if (Pi > 0.5) user->res_n_warn++;
    if (Pi > 1.0) user->res_n_bad++;
  }

  if (startup) {
    /* not judged: see res_t_skip in NASA_types.h */
  } else if (Pi > 1.0) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      "!!! RESOLUTION WARNING step %d: interface too fast for this eps/mesh (K&P Eq. 45)\n"
      "!!!   v_n = %.3e m/s  ->  Pi = beta0 v_n eps/d0 = %.3f > 1\n"
      "!!!   eps must be below d0/(beta0 v_n) = %.3e m; rerun comp_eps.py --vn %.3e\n"
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n",
      (int)step, (double)vmax, (double)Pi,
      (double)(user->d0_sub0 / (user->beta_sub0 * vmax)), (double)vmax);
  } else if (Pi > 0.5) {
    PetscPrintf(PETSC_COMM_WORLD,
      "RESOLUTION NOTE step %d: Pi = %.3f (outside comp_eps margin 0.5, inside the Eq. 45 bound), v_n = %.3e m/s\n",
      (int)step, (double)Pi, (double)vmax);
  }
  if (open_max > 0.01) {
    PetscPrintf(PETSC_COMM_WORLD,
      "!!! RESERVOIR-FACE WARNING step %d: phi = %.3f on a vapor-reservoir face (ice has "
      "reached the reservoir; the setup assumes it never does)\n", (int)step, (double)open_max);
  }

  if (user->res_view == NULL) {
    char fn[256];
    const char *dir = getenv("folder");
    if (dir) { snprintf(fn, sizeof(fn), "%s/resolution.dat", dir); }
    else     { snprintf(fn, sizeof(fn), "resolution.dat"); }
    PetscBool app;
    ierr = OpenLogViewer(user, fn, &user->res_view, &app); CHKERRQ(ierr);
    if (!app) {
    ierr = PetscViewerASCIIPrintf(user->res_view,
             "# K&P Eq. 45 check: Pi = beta0 v_n eps/d0 (<= 0.5 design, > 1 violated); "
             "eps %.6e m, beta0 %.6e s/m, d0 %.6e m\n# t step dt vn_max Pi phi_max_reservoir_face\n",
             (double)user->eps, (double)user->beta_sub0, (double)user->d0_sub0); CHKERRQ(ierr);
    }
  }
  ierr = PetscViewerASCIIPrintf(user->res_view, "%.10e %d %.6e %.6e %.6e %.6e\n",
           (double)t, (int)step, (double)dt, (double)vmax, (double)Pi, (double)open_max); CHKERRQ(ierr);
  ierr = PetscViewerFlush(user->res_view); CHKERRQ(ierr);

  ierr = VecCopy(U, user->res_U_prev); CHKERRQ(ierr);
  user->res_t_prev = t;
  PetscFunctionReturn(0);
}


PetscErrorCode ResolutionStaticCheck(AppCtx *user)
{
  PetscFunctionBegin;
  const PetscReal hmax = user->eps / PetscSqrtReal(2.0);
  PetscReal L[3] = {user->Lx, user->Ly, user->Lz};
  PetscInt  N[3] = {user->Nx, user->Ny, user->Nz};
  const char ax[3] = {'x', 'y', 'z'};
  PetscBool bad = PETSC_FALSE;
  PetscPrintf(PETSC_COMM_WORLD, "\nRESOLUTION CHECK (start-up)\n");
  for (PetscInt l = 0; l < user->dim; l++) {
    if (N[l] <= 0 || L[l] <= 0.0) continue;
    PetscReal h = L[l] / (PetscReal)N[l];
    PetscPrintf(PETSC_COMM_WORLD, "   h_%c = L/N = %.4e m   eps/sqrt(2) = %.4e m   ratio %.3f\n",
                ax[l], (double)h, (double)hmax, (double)(h / hmax));
    if (h > 1.05 * hmax) bad = PETSC_TRUE;
  }
  PetscReal rmin = PETSC_MAX_REAL;
  for (PetscInt k = 0; k < user->n_act; k++)
    if (user->radius && user->radius[k] > 0.0) rmin = PetscMin(rmin, user->radius[k]);
  if (rmin < PETSC_MAX_REAL) {
    PetscPrintf(PETSC_COMM_WORLD, "   eps / smallest IC grain radius = %.4f (comp_eps limit 0.05)\n",
                (double)(user->eps / rmin));
    if (user->eps / rmin > 0.05) bad = PETSC_TRUE;
  }
  PetscPrintf(PETSC_COMM_WORLD,
    "   eps is only valid for front speeds v_n < d0/(2 beta0 eps) = %.3e m/s (comp_eps margin)\n",
    (double)(user->d0_sub0 / (2.0 * user->beta_sub0 * user->eps)));
  if (bad) {
    PetscPrintf(PETSC_COMM_WORLD,
      "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      "!!! RESOLUTION WARNING: mesh coarser than eps/sqrt(2) by > 5 %% or eps > 5 %% of the\n"
      "!!! smallest grain radius. Size eps and Nx/Ny/Nz with preprocess/comp_eps.py.\n"
      "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");
  }
  PetscFunctionReturn(0);
}


PetscErrorCode ResolutionSummary(AppCtx *user)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,
    "\nRESOLUTION SUMMARY (after the %.3g s start-up window; max Pi inside it %.3f, not judged): "
    "max Pi = %.3f; steps with Pi > 0.5: %d, Pi > 1: %d; "
    "max phi on a reservoir face %.3e (%d steps > 0.01)%s\n",
    (double)user->res_t_skip, (double)user->res_pi_max_startup,
    (double)user->res_pi_max, (int)user->res_n_warn, (int)user->res_n_bad,
    (double)user->res_open_phi_max, (int)user->res_n_open,
    (user->res_n_bad > 0 || user->res_n_open > 0) ? "   <-- CHECK THIS RUN" : "");
  PetscFunctionReturn(0);
}
