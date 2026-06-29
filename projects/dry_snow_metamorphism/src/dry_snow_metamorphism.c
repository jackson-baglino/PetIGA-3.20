#include "NASA_main.h"

int main(int argc, char *argv[]) {

  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,NULL,NULL);CHKERRQ(ierr);

  PetscLogDouble itim;
  ierr = PetscTime(&itim); CHKERRQ(ierr);

  AppCtx user;
  PetscInt flag_sedgrav, flag_BC_Tfix, flag_BC_rhovfix;

  // ---- Default physical parameters (overridable via .opts / CLI) ----
  user.xi_v       = 1.0e-5;
  user.xi_T       = 1.0e-4;
  user.flag_xiT   = 1;
  user.Lambd      = 1.0;
  user.air_lim    = 1.0e-6;
  user.nsteps_IC  = 10;
  user.lat_sub    = 2.83e6;
  user.thcond_ice = 2.29;
  user.thcond_met = 36.0;
  user.thcond_air = 0.02;
  user.cp_ice     = 1.96e3;
  user.cp_met     = 4.86e2;
  user.cp_air     = 1.044e3;
  user.rho_ice    = 919.0;
  user.rho_met    = 7753.0;
  user.rho_air    = 1.341;
  user.dif_vap    = 2.178e-5;
  user.T_melt     = 0.0;
  user.flag_it0   = 1;
  user.flag_tIC   = 0;
  user.flag_Tdep  = 1;
  user.d0_sub0    = 1.0e-9;
  user.beta_sub0  = 1.4e5;
  // Default initial conditions
  user.temp0         = -10.0;
  user.hum0          = 1.0;
  user.eps           = 1.0e-6;
  user.readFlag      = 1;
  user.grad_temp0[0] = 0.0;
  user.grad_temp0[1] = 0.0;
  user.grad_temp0[2] = 0.0;
  user.output_dir[0] = '\0';
  user.grains_file[0] = '\0';
  user.initial_cond[0] = '\0';
  user.PFgeom[0] = '\0';
  // Default I/O paths
  PetscStrcpy(user.output_dir, ".");

  // ---- Parse all parameters from .opts files / CLI ----
  PetscInt  Nx=100, Ny=100, Nz=1;
  PetscReal Lx=1.0e-3, Ly=1.0e-3, Lz=1.0e-4;
  PetscReal delt_t=1.0e-4, t_final=1.0e6;
  PetscInt  n_out=100, dim=2, p=1, C=0;
  PetscBool output=PETSC_TRUE, monitor=PETSC_TRUE;

  ierr = GetOptions(&user, &Nx, &Ny, &Nz, &Lx, &Ly, &Lz,
                    &delt_t, &t_final, &n_out,
                    &dim, &p, &C, &output, &monitor); CHKERRQ(ierr);

  user.p = p; user.C = C;

  // ---- Gibbs-Thomson / kinetic parameters ----
  PetscReal gamma_im = 0.033, gamma_iv = 0.109, gamma_mv = 0.056;
  PetscReal rho_rhovs, rhoI_vs, d_rhovs;
  RhoVS_I(&user, user.temp0, &rhoI_vs, &d_rhovs);
  rho_rhovs = user.rho_ice / rhoI_vs;

  ierr = PetscPrintf(PETSC_COMM_WORLD, "rho_rhovs = %e \n", rho_rhovs);CHKERRQ(ierr);

  // ---- Grain configuration ----
  if (user.readFlag == 0) {
    user.NCsed = 20;
  } else {
    user.NCsed = 0;
  }
  flag_sedgrav    = 0;
  user.RCsed      = 0.8e-5;
  user.RCsed_dev  = 0.4;
  user.NCice      = 140;
  user.RCice      = 0.9e-4;
  user.RCice_dev  = 0.55;

  // ---- Boundary conditions ----
  user.periodic   = 0;
  flag_BC_Tfix    = 1;
  flag_BC_rhovfix = 0;
  if(user.periodic==1 && flag_BC_Tfix==1) flag_BC_Tfix=0;
  if(user.periodic==1 && flag_BC_rhovfix==1) flag_BC_rhovfix=0;

  // ---- Output timing ----
  user.outp     = 0;
  user.t_out    = 0;
  user.t_interv = t_final / (n_out - 1);

  PetscInt adap = 1;
  PetscInt NRmin = 2, NRmax = 5;
  PetscReal factor = pow(10.0, 1.0/8.0);
  PetscReal dtmin = 0.01*delt_t, dtmax = 0.5*user.t_interv;
  PetscInt max_rej = 10;
  if(adap==1) PetscPrintf(PETSC_COMM_WORLD,"Adaptive time stepping: NR_iter %d-%d  factor %.3f  dt0 %.2e  dt_range %.2e-%.2e  \n\n",NRmin,NRmax,factor,delt_t,dtmin,dtmax);

  PetscInt size;
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  PetscPrintf(PETSC_COMM_WORLD, "Running on %d processes.\n\n\n", size);

  // ---- Gibbs-Thomson kinetic parameters ----
  user.diff_sub = 0.5*(user.thcond_air/user.rho_air/user.cp_air + user.thcond_ice/user.rho_ice/user.cp_ice);
  user.Etai = gamma_iv + gamma_im - gamma_mv;
  user.Etam = gamma_mv + gamma_im - gamma_iv;
  user.Etaa = gamma_iv + gamma_mv - gamma_im;
  PetscReal a1=5.0, a2=0.1581;
  PetscReal d0_sub, beta_sub, lambda_sub, tau_sub;
  d0_sub     = user.d0_sub0 / rho_rhovs;
  beta_sub   = user.beta_sub0 / rho_rhovs;
  lambda_sub = a1 * user.eps / d0_sub;
  tau_sub    = user.eps * lambda_sub * (beta_sub/a1 + a2*user.eps/user.diff_sub + a2*user.eps/user.dif_vap);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Beta_sub = %e \n", beta_sub);CHKERRQ(ierr);

  user.mob_sub  = 1.0 * user.eps / 3.0 / tau_sub;
  user.alph_sub = 1.0 * lambda_sub / tau_sub;

  if(user.flag_Tdep==0) PetscPrintf(PETSC_COMM_WORLD,"FIXED PARAMETERS: tau %.4e  lambda %.4e  M0 %.4e  alpha %.4e \n\n",tau_sub,lambda_sub,user.mob_sub,user.alph_sub);
  else PetscPrintf(PETSC_COMM_WORLD,"TEMPERATURE DEPENDENT G-T PARAMETERS \n\n");

  // ---- IGA setup ----
  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,dim);CHKERRQ(ierr);
  ierr = IGASetDof(iga,3);CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,0,"phaseice"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,1,"temperature"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,2,"vap_density"); CHKERRQ(ierr);

  PetscInt l, m;
  IGAAxis axis0, axis1, axis2;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis0,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,Nx,0.0,Lx,C);CHKERRQ(ierr);
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis1,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,Ny,0.0,Ly,C);CHKERRQ(ierr);
  if(dim==3){
    ierr = IGAGetAxis(iga,2,&axis2);CHKERRQ(ierr);
    if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis2,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axis2,p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis2,Nz,0.0,Lz,C);CHKERRQ(ierr);
  }

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;

  PetscInt nmb = iga->elem_width[0]*iga->elem_width[1]*SQ(p+1);
  if(dim==3) nmb = iga->elem_width[0]*iga->elem_width[1]*iga->elem_width[2]*CU(p+1);
  ierr = PetscMalloc(sizeof(PetscReal)*(nmb),&user.Phi_sed);CHKERRQ(ierr);
  ierr = PetscMemzero(user.Phi_sed,sizeof(PetscReal)*(nmb));CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*(nmb),&user.alph);CHKERRQ(ierr);
  ierr = PetscMemzero(user.alph,sizeof(PetscReal)*(nmb));CHKERRQ(ierr);
  ierr = PetscMalloc(sizeof(PetscReal)*(nmb),&user.mob);CHKERRQ(ierr);
  ierr = PetscMemzero(user.mob,sizeof(PetscReal)*(nmb));CHKERRQ(ierr);

  ierr = IGASetFormIFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIJacobian(iga,Jacobian,&user);CHKERRQ(ierr);

  // ---- Boundary conditions ----
  if(flag_BC_rhovfix==1){
    PetscReal rho0_vs;
    RhoVS_I(&user,user.temp0,&rho0_vs,NULL);
    for(l=0;l<dim;l++) for(m=0;m<2;m++) ierr = IGASetBoundaryValue(iga,l,m,2,user.hum0*rho0_vs);CHKERRQ(ierr);
  }
  if(flag_BC_Tfix==1){
    PetscReal T_BC[dim][2], LL[dim];
    LL[0] = Lx; LL[1]=Ly; LL[2]=Lz;
    for(l=0;l<dim;l++) for(m=0;m<2;m++) T_BC[l][m] = user.temp0 + (2.0*m-1)*user.grad_temp0[l]*0.5*LL[l];
    for(l=0;l<dim;l++) for(m=0;m<2;m++) ierr = IGASetBoundaryValue(iga,l,m,1,T_BC[l][m]);CHKERRQ(ierr);
  }

  // ---- Time stepping ----
  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts,t_final);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,delt_t);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSALPHA);CHKERRQ(ierr);
  ierr = TSAlphaSetRadius(ts,0.5);CHKERRQ(ierr);
  if (monitor) {ierr = TSMonitorSet(ts,Monitor,&user,NULL);CHKERRQ(ierr);}
  if (output)  {ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);}
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ts->adap = adap;
  ts->NRmin = NRmin;
  ts->NRmax = NRmax;
  ts->factor = factor;
  ts->dtmax = dtmax;
  ts->dtmin = dtmin;
  ts->max_reject = max_rej;
  ts->max_snes_failures = -1;

  SNES nonlin;
  ierr = TSGetSNES(ts,&nonlin);CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(nonlin,SNESDOFConvergence,&user,NULL);CHKERRQ(ierr);

  // ---- Sediment grains (if any) ----
  if(user.NCsed>0){
    if(flag_sedgrav==1) {
      if(dim==2) {ierr = InitialSedGrainsGravity(iga,&user);CHKERRQ(ierr);}
      else {
        PetscPrintf(PETSC_COMM_WORLD,"Pluviation Script not prepared for 3D. Run no-gravity \n");
        ierr = InitialSedGrains(iga,&user);CHKERRQ(ierr);
      }
    } else {ierr = InitialSedGrains(iga,&user);CHKERRQ(ierr);}

    IGA igaS;   IGAAxis axis0S, axis1S, axis2S;
    ierr = IGACreate(PETSC_COMM_WORLD,&igaS);CHKERRQ(ierr);
    ierr = IGASetDim(igaS,dim);CHKERRQ(ierr);
    ierr = IGASetDof(igaS,1);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaS,0,&axis0S);CHKERRQ(ierr);
    if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis0S,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axis0S,p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis0S,Nx,0.0,Lx,C);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaS,1,&axis1S);CHKERRQ(ierr);
    if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis1S,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axis1S,p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis1S,Ny,0.0,Ly,C);CHKERRQ(ierr);
    if(dim==3){
      ierr = IGAGetAxis(igaS,2,&axis2S);CHKERRQ(ierr);
      if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis2S,PETSC_TRUE);CHKERRQ(ierr);}
      ierr = IGAAxisSetDegree(axis2S,p);CHKERRQ(ierr);
      ierr = IGAAxisInitUniform(axis2S,Nz,0.0,Lz,C);CHKERRQ(ierr);
    }
    ierr = IGASetFromOptions(igaS);CHKERRQ(ierr);
    ierr = IGASetUp(igaS);CHKERRQ(ierr);

    Vec S;
    ierr = IGACreateVec(igaS,&S);CHKERRQ(ierr);
    if(dim==2) {ierr = FormInitialSoil2D(igaS,S,&user);CHKERRQ(ierr);}
    else {ierr = FormInitialSoil3D(igaS,S,&user);CHKERRQ(ierr);}

    char filename[256], filevect[256];
    sprintf(filename, "%s/igasoil.dat", user.output_dir);
    ierr=IGAWrite(igaS,filename);CHKERRQ(ierr);
    sprintf(filevect, "%s/soil.dat", user.output_dir);
    ierr=IGAWriteVec(igaS,S,filevect);CHKERRQ(ierr);

    ierr = VecDestroy(&S);CHKERRQ(ierr);
    ierr = IGADestroy(&igaS);CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"No sed grains\n\n");
    user.n_actsed = 0;
  }

  ierr = InitialIceGrains(iga,&user);CHKERRQ(ierr);

  // ---- Print resolved parameters ----
  PetscPrintf(PETSC_COMM_WORLD, "Nx: %f\n",           user.Nx);
  PetscPrintf(PETSC_COMM_WORLD, "Ny: %f\n",           user.Ny);
  PetscPrintf(PETSC_COMM_WORLD, "Nz: %f\n",           user.Nz);
  PetscPrintf(PETSC_COMM_WORLD, "Lx: %f\n",           user.Lx);
  PetscPrintf(PETSC_COMM_WORLD, "Ly: %f\n",           user.Ly);
  PetscPrintf(PETSC_COMM_WORLD, "Lz: %f\n",           user.Lz);
  PetscPrintf(PETSC_COMM_WORLD, "delt_t: %f\n",       delt_t);
  PetscPrintf(PETSC_COMM_WORLD, "t_final: %f\n",      t_final);
  PetscPrintf(PETSC_COMM_WORLD, "humidity: %f\n",     user.hum0);
  PetscPrintf(PETSC_COMM_WORLD, "temp: %f\n",         user.temp0);
  PetscPrintf(PETSC_COMM_WORLD, "grad_temp0X: %f\n",  user.grad_temp0[0]);
  PetscPrintf(PETSC_COMM_WORLD, "grad_temp0Y: %f\n",  user.grad_temp0[1]);
  PetscPrintf(PETSC_COMM_WORLD, "grad_temp0Z: %f\n",  user.grad_temp0[2]);
  PetscPrintf(PETSC_COMM_WORLD, "dim: %d\n",          user.dim);
  PetscPrintf(PETSC_COMM_WORLD, "eps: %e\n",          user.eps);
  PetscPrintf(PETSC_COMM_WORLD, "output_dir: %s\n",   user.output_dir);
  PetscPrintf(PETSC_COMM_WORLD, "grains_file: %s\n",  user.grains_file);

  // ---- Initial condition ----
  PetscReal t=0; Vec U;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = VecZeroEntries(U);CHKERRQ(ierr);
  if(dim==2) {
    ierr = FormLayeredInitialCondition2D(iga,t,U,&user,user.initial_cond,user.PFgeom);CHKERRQ(ierr);
  }
  else {ierr = FormInitialCondition3D(iga,t,U,&user,user.initial_cond,user.PFgeom);CHKERRQ(ierr);}

  ierr = TSSolve(ts,U);CHKERRQ(ierr);

  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);

  ierr = PetscFree(user.Phi_sed);CHKERRQ(ierr);
  ierr = PetscFree(user.alph);CHKERRQ(ierr);
  ierr = PetscFree(user.mob);CHKERRQ(ierr);

  PetscLogDouble ltim, tim;
  ierr = PetscTime(&ltim); CHKERRQ(ierr);
  tim = ltim - itim;

  int days    = (int)(tim / 86400);
  int hours   = (int)((tim - days * 86400) / 3600);
  int minutes = (int)((tim - days * 86400 - hours * 3600) / 60);
  double seconds = tim - days * 86400 - hours * 3600 - minutes * 60;

  PetscPrintf(PETSC_COMM_WORLD, "\033[1mcomp time: ");
  if (days > 0)
    PetscPrintf(PETSC_COMM_WORLD, "%d day%s ", days, days == 1 ? "" : "s");
  if (hours > 0)
    PetscPrintf(PETSC_COMM_WORLD, "%d hour%s ", hours, hours == 1 ? "" : "s");
  if (minutes > 0)
    PetscPrintf(PETSC_COMM_WORLD, "%d min%s ", minutes, minutes == 1 ? "" : "s");
  if (seconds > 0 || (days == 0 && hours == 0 && minutes == 0))
    PetscPrintf(PETSC_COMM_WORLD, "%.2f sec", seconds);
  PetscPrintf(PETSC_COMM_WORLD, "\033[0m\n\n");

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
