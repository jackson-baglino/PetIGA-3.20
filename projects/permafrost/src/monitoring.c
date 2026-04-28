#include "monitoring.h"
#include "assembly.h"
#include "material_properties.h"
#include <petscstring.h>

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

  if(user->flag_Tdep==1){
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
            d0 = 2.548e-7/(solS[1]+273.15);

            if(bet0>bet_max) bet_max = bet0;
            if(bet0<bet_min) bet_min = bet0;
            d0_sub   = d0/rho_rhovs;  
            beta_sub = bet0/rho_rhovs;
            lambda_sub = a1*user->eps/d0_sub;
            tau_sub    = user->eps*lambda_sub*(beta_sub/a1 + a2*user->eps/user->diff_sub + a2*user->eps/user->dif_vap);

            user->mob[indd]         = user->eps/3.0/tau_sub;
            user->mob_sed_arr[indd] = user->mob_sed;  /* placeholder: update when sediment RHS is defined */
            user->alph[indd]        = lambda_sub/tau_sub;

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
    user->flag_Tdep = 0;

    // Print the new mobility and alpha values for verification
    PetscPrintf(PETSC_COMM_WORLD, "M0_sub new: %.6e\n", user->mob_sub);
    PetscPrintf(PETSC_COMM_WORLD, "alpha_sub new: %.6e\n", user->alph_sub);
  }


  //-------- domain integrals
  PetscScalar stats[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  ierr = IGAComputeScalar(user->iga,U,9,&stats[0],Integration,mctx);CHKERRQ(ierr);
  PetscReal tot_ice        = PetscRealPart(stats[0]);
  PetscReal tot_trip       = PetscRealPart(stats[1]);
  PetscReal tot_air        = PetscRealPart(stats[2]);
  PetscReal tot_temp       = PetscRealPart(stats[3]);
  PetscReal tot_rhov       = PetscRealPart(stats[4]);
  PetscReal sub_interf     = PetscRealPart(stats[5]);
  PetscReal tot_sed        = PetscRealPart(stats[6]);
  PetscReal sed_air_interf = PetscRealPart(stats[7]);
  PetscReal ice_sed_interf = PetscRealPart(stats[8]);
 
  //-------- phase-field min/max bounds (printed every step for out-of-bounds detection)
  {
    PetscReal phi_ice_min =  1.0e30, phi_ice_max = -1.0e30;
    PetscReal phi_sed_min =  1.0e30, phi_sed_max = -1.0e30;
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
        PetscReal fs = PetscRealPart(solb[3]);
        PetscReal fa = 1.0 - fi - fs;
        if (fi < phi_ice_min) phi_ice_min = fi;
        if (fi > phi_ice_max) phi_ice_max = fi;
        if (fs < phi_sed_min) phi_sed_min = fs;
        if (fs > phi_sed_max) phi_sed_max = fs;
        if (fa < phi_air_min) phi_air_min = fa;
        if (fa > phi_air_max) phi_air_max = fa;
      }
      ierr = IGAElementEndPoint(elementb, &pointb); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &elementb); CHKERRQ(ierr);
    ierr = IGARestoreLocalVecArray(user->iga, U, &localUb, &arrayUb); CHKERRQ(ierr);

    PetscReal Gfi_min, Gfi_max, Gfs_min, Gfs_max, Gfa_min, Gfa_max;
    ierr = MPI_Allreduce(&phi_ice_min, &Gfi_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_ice_max, &Gfi_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_sed_min, &Gfs_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_sed_max, &Gfs_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_air_min, &Gfa_min, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&phi_air_max, &Gfa_max, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD,
        "  BOUNDS: phi_ice [%.4f, %.4f]  phi_sed [%.4f, %.4f]  phi_air [%.4f, %.4f]\n",
        Gfi_min, Gfi_max, Gfs_min, Gfs_max, Gfa_min, Gfa_max);

    /* Abort if any phase field exceeds the configured bounds */
    PetscBool oob = (Gfi_min < user->phase_lo || Gfi_max > user->phase_hi ||
                     Gfs_min < user->phase_lo || Gfs_max > user->phase_hi ||
                     Gfa_min < user->phase_lo || Gfa_max > user->phase_hi);
    if (oob) {
      PetscPrintf(PETSC_COMM_WORLD,
          "\033[31m[ABORT] Phase field out of bounds [%.2f, %.2f] at step %d — stopping TS.\033[0m\n",
          user->phase_lo, user->phase_hi, step);
      ierr = TSSetConvergedReason(ts, TS_DIVERGED_NONLINEAR_SOLVE); CHKERRQ(ierr);
    }
  }

  //-------------
  PetscReal dt;
  TSGetTimeStep(ts,&dt);
  if(step==1) user->flag_it0 = 0;

  //------------- initial condition
  if(user->flag_tIC==1) if(step==user->nsteps_IC) {
    user->flag_tIC = 0; user->t_IC = t; //user->flag_rtol = 1;
    PetscPrintf(PETSC_COMM_WORLD,"INITIAL_CONDITION!!! \n");
  }

  //------------- sediment freeze (mode 1: triggered by interface stability OR nsteps_sed fallback)
  if (user->flag_sed_mode == 1 && !user->flag_sed_frozen) {
    /* Compute relative per-step change in ice-sed interface integral.
     * When this drops below sed_freeze_tol the interface has relaxed to the
     * 3-phase equilibrium and the sediment penalty can safely be activated. */
    PetscReal delta_icesd = PetscAbsReal(ice_sed_interf - user->prev_ice_sed_interf);
    PetscReal rel_change  = (user->prev_ice_sed_interf > 0.0)
                             ? delta_icesd / user->prev_ice_sed_interf
                             : 1.0;   /* not yet initialized — skip until step > 1 */

    PetscBool interface_stable = (step > 1 && rel_change < user->sed_freeze_tol);
    PetscBool fallback_reached = (t >= user->t_sed_freeze);

    if (interface_stable || fallback_reached) {
      user->flag_sed_frozen = 1;
      const char *reason = interface_stable ? "interface stabilized" : "t_sed_freeze reached";

      /* Snapshot the relaxed sediment field as the penalty reference.
       * Phi_sed0[] was populated at init from the tanh IC (poorly-defined
       * interface).  Overwrite it now with the current 3-phase equilibrium
       * profile so the penalty k_sed*(sed - sed_ref) drives sediment back
       * toward its properly-relaxed state, not the original IC. */
      {
        Vec localU_snap; const PetscScalar *arrayU_snap;
        IGAElement element_snap; IGAPoint point_snap; PetscScalar *UU_snap;
        PetscInt indd_snap = 0;
        ierr = IGAGetLocalVecArray(user->iga, U, &localU_snap, &arrayU_snap); CHKERRQ(ierr);
        ierr = IGABeginElement(user->iga, &element_snap); CHKERRQ(ierr);
        while (IGANextElement(user->iga, element_snap)) {
          ierr = IGAElementGetValues(element_snap, arrayU_snap, &UU_snap); CHKERRQ(ierr);
          ierr = IGAElementBeginPoint(element_snap, &point_snap); CHKERRQ(ierr);
          while (IGAElementNextPoint(element_snap, point_snap)) {
            PetscScalar sol_snap[4];
            ierr = IGAPointFormValue(point_snap, UU_snap, &sol_snap[0]); CHKERRQ(ierr);
            user->Phi_sed0[indd_snap] = PetscRealPart(sol_snap[3]);
            indd_snap++;
          }
          ierr = IGAElementEndPoint(element_snap, &point_snap); CHKERRQ(ierr);
        }
        ierr = IGAEndElement(user->iga, &element_snap); CHKERRQ(ierr);
        ierr = IGARestoreLocalVecArray(user->iga, U, &localU_snap, &arrayU_snap); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "  sed0 reference updated to relaxed profile (%d quadrature points)\n", indd_snap);
      }

      TSSetTimeStep(ts, 1.0e-2);
      PetscPrintf(PETSC_COMM_WORLD, "\033[34m"
          "╔════════════════════════════════════════════════════════════╗\n"
          "║              SEDIMENT PENALTY ACTIVATED                    ║\n"
          "║  3-phase ice eq. unchanged; sediment pinned to sed0        ║\n"
          "╠════════════════════════════════════════════════════════════╣\n"
          "║  Step: %-5d  Reason: %-30s   ║\n"
          "║  Rel. change in ice-sed interface: %.2e                    ║\n"
          "║  New time step: %.2e s                                      ║\n"
          "╚════════════════════════════════════════════════════════════╝\n"
          "\033[0m\n",
          step, reason, rel_change, dt);
    }
  }
  /* Always update the previous interface value for next step's comparison */
  user->prev_ice_sed_interf = ice_sed_interf;

  //------ printf information (robust table header + aligned columns)
  {
    const PetscInt headerEvery = 1;  // print header every step (can change to 10 later)
    PetscErrorCode ierr2;

    if ((step % headerEvery) == 0) {
      // Build header line using the SAME fixed-width fields as the data row
      char header[512];
      ierr2 = PetscSNPrintf(header, sizeof(header),
                            " %5s | %12s | %9s | %10s | %10s | %10s | %9s | %9s | %10s | %10s",
                            "STEP",
                            "TIME [s]",
                            "DT [s]",
                            "TOT_ICE",
                            "TOT_AIR",
                            "TOT_SED",
                            "TEMP",
                            "TOT_RHOV",
                            "I-A INTERF",
                            "TRIPL_JUNC");
      CHKERRQ(ierr2);

      // Separator line: exactly matches header length
      PetscInt hlen = 0;
      ierr2 = PetscStrlen(header, &hlen);
      CHKERRQ(ierr2);

      char sep[512];
      PetscInt i;
      for (i = 0; i < hlen && i < (PetscInt)sizeof(sep) - 1; i++) sep[i] = '-';
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
                " %5d | %12.5e | %9.3e | %10.3e | %10.3e | %10.3e | %9.3e | %9.3e | %10.3e | %10.3e\n",
                step, t, dt,
                tot_ice, tot_air, tot_sed, tot_temp, tot_rhov,
                sub_interf, tot_trip);

    // Optional: add a blank line every N rows for readability (set to 0 to disable)
    // if (step > 0 && (step % 25) == 0) PetscPrintf(PETSC_COMM_WORLD, "\n");
  }

  //-------- write per-step interface metrics for relaxation analysis
  {
    char relax_file[256];
    const char *renv = "folder"; char *rdir; rdir = getenv(renv);
    sprintf(relax_file, "%s/relax_monitor.dat", rdir);

    PetscViewer rv;
    PetscViewerCreate(PETSC_COMM_WORLD, &rv);
    PetscViewerSetType(rv, PETSCVIEWERASCII);
    if (step == 0) {
      PetscViewerFileSetMode(rv, FILE_MODE_WRITE);
      PetscViewerFileSetName(rv, relax_file);
      PetscViewerASCIIPrintf(rv,
          "# step  t  tot_ice  tot_sed  ice_air_interf  sed_air_interf  ice_sed_interf  tot_trip\n");
    } else {
      PetscViewerFileSetMode(rv, FILE_MODE_APPEND);
      PetscViewerFileSetName(rv, relax_file);
    }
    PetscViewerASCIIPrintf(rv, "%d %e %e %e %e %e %e %e\n",
        step, t,
        tot_ice, tot_sed,
        sub_interf, sed_air_interf, ice_sed_interf,
        tot_trip);
    PetscViewerDestroy(&rv);
  }

  PetscInt print=0;
  if(user->outp > 0) {
    if(step % user->outp == 0) print=1;
  } else {
    if (t>= user->t_out) print=1;
  }

  print = 1;

  if(print==1) {
    char filedata[256];
    const char *env = "folder"; char *dir; dir = getenv(env);

    sprintf(filedata,"%s/SSA_evo.dat",dir);
    PetscViewer       view;
    PetscViewerCreate(PETSC_COMM_WORLD,&view);
    PetscViewerSetType(view,PETSCVIEWERASCII);

    if (step==0){
      PetscViewerFileSetMode(view,FILE_MODE_WRITE);
    } else {
      PetscViewerFileSetMode(view,FILE_MODE_APPEND);
    }

    PetscViewerFileSetName(view,filedata);
    PetscViewerASCIIPrintf(view,"%e %e %e %d\n",sub_interf/user->eps, tot_ice, t, step);

    PetscViewerDestroy(&view);
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
    user->t_out += user->t_interv;

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