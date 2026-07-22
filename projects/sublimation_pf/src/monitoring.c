#include "monitoring.h"
#include "assembly.h"
#include "material_properties.h"
#include <petscstring.h>
#include <petsc/private/tsimpl.h>   /* for direct access to ts->dtmin */

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  AppCtx *user = (AppCtx *)mctx;

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
            PetscScalar solS[4];
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
            /* Capillary length from the model's own surface energy (K&P
             * Eq. 13): d0 = gamma_iv*V_m/(R*T). gamma_iv = user->Sigma_i (the
             * -Sigma_i opts parameter, default 0.109 J/m^2), V_m = 1.963e-5
             * m^3/mol, R = 8.314 J/mol/K — the same constants comp_eps.py
             * uses, so preprocessing and runtime agree. Replaces the stale
             * hardcoded 2.548e-7/T_K (~1% low, and blind to -Sigma_i). */
            d0 = user->Sigma_i * 1.963e-5 / 8.314 / (solS[1]+273.15);

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
        PetscScalar solb[4];
        ierr = IGAPointFormValue(pointb, UUb, &solb[0]); CHKERRQ(ierr);
        PetscReal fi = PetscRealPart(solb[0]);
        PetscReal fs = (user->dof == 4) ? PetscRealPart(solb[3]) : 0.0;
        PetscReal fa = 1.0 - fi - fs;  /* 3-phase air fraction */
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

    PetscPrintf(PETSC_COMM_WORLD,
        "  BOUNDS: phi_ice [%.4f, %.4f]  phi_air [%.4f, %.4f]\n",
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

  PetscInt print=0;
  if(user->outp > 0) {
    if(step % user->outp == 0) print=1;
  } else {
    if (t>= user->t_out) print=1;
  }

  print = 1;

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
      ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filedata,&user->ssa_view);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(user->ssa_view,"%e %e %e %d %e %e %e %e\n",
                                  sub_interf/user->eps, tot_ice, t, step, dt,
                                  tot_air, tot_rhov, tot_mass);CHKERRQ(ierr);
    ierr = PetscViewerFlush(user->ssa_view);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}


PetscErrorCode OutputMonitor(TS ts, PetscInt step, PetscReal t, Vec U, 
                              void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;  

  // Check if it's the first step
  if (step == 0) {
    const char *env = "folder";
    char *dir;
    dir = getenv(env);

    char fileiga[256];
    sprintf(fileiga, "%s/igasol.dat", dir);

    ierr = IGAWrite(user->iga, fileiga);CHKERRQ(ierr);
  }

  // Check if it's time to print output
  PetscInt print = 0;
  if (user->outp > 0) {   // Print output every user->outp steps
    if (step % user->outp == 0) print = 1;
  } 
  else {                  // Print output every user->t_interv seconds
    if (t >= user->t_out) print = 1;
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
 * Runs as a TS monitor (accepted steps only), so the rate is always
 * measured on a solution the solver already accepted; the cap is applied
 * with TSSetTimeStep after the controller's own adjustment.
 * ========================================================================= */
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
    ierr = VecCopy(U, user->cfl_U_prev);       CHKERRQ(ierr);
    user->cfl_t_prev = t;
    PetscFunctionReturn(0);
  }

  PetscReal dt_used = t - user->cfl_t_prev;
  if (dt_used > 0.0) {
    /* ||dphi||_inf of the last step: reuse cfl_U_prev as the diff buffer */
    PetscReal dphi_max;
    ierr = VecAYPX(user->cfl_U_prev, -1.0, U); CHKERRQ(ierr);  /* prev = U - prev */
    ierr = VecStrideNorm(user->cfl_U_prev, 0, NORM_INFINITY, &dphi_max); CHKERRQ(ierr);

    if (dphi_max > 0.0) {
      PetscReal dt_cap = user->cfl_dphimax * dt_used / dphi_max;
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
 * Registered via TSSetPreStep() in permafrost2.c. The AppCtx pointer is
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
    user->bounds_violated = PETSC_FALSE;
    user->bounds_new_dt   = 0.0;
  }

  PetscFunctionReturn(0);
}