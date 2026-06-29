#include "options_helper.h"

PetscErrorCode GetOptions(AppCtx *user,
                          PetscInt *Nx, PetscInt *Ny, PetscInt *Nz,
                          PetscReal *Lx, PetscReal *Ly, PetscReal *Lz,
                          PetscReal *delt_t, PetscReal *t_final, PetscInt *n_out,
                          PetscInt *dim, PetscInt *p, PetscInt *C,
                          PetscBool *output, PetscBool *monitor)
{
  PetscErrorCode ierr;
  PetscInt ngrad = 3;
  PetscBool flg;

  PetscOptionsBegin(PETSC_COMM_WORLD, "", "DSM Options", "IGA");

  // Mesh & domain
  ierr = PetscOptionsInt("-Nx",  "elements in x",         __FILE__, *Nx,  Nx,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Ny",  "elements in y",         __FILE__, *Ny,  Ny,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Nz",  "elements in z",         __FILE__, *Nz,  Nz,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Lx", "domain length x [m]",  __FILE__, *Lx,  Lx,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Ly", "domain length y [m]",  __FILE__, *Ly,  Ly,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Lz", "domain length z [m]",  __FILE__, *Lz,  Lz,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-dim", "spatial dimension",     __FILE__, *dim, dim, NULL);CHKERRQ(ierr);

  // IGA basis
  ierr = PetscOptionsInt("-p", "B-spline polynomial order",  __FILE__, *p, p, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-C", "global continuity order",    __FILE__, *C, C, NULL);CHKERRQ(ierr);

  // Time stepping
  ierr = PetscOptionsReal("-delt_t",  "initial time step [s]",      __FILE__, *delt_t,  delt_t,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-t_final", "final simulation time [s]",  __FILE__, *t_final, t_final, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-n_out",    "number of output frames",     __FILE__, *n_out,   n_out,   NULL);CHKERRQ(ierr);

  // Physics
  ierr = PetscOptionsReal("-temp",     "initial temperature [C]",      __FILE__, user->temp0,       &user->temp0,       NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-humidity", "relative humidity [0-1]",      __FILE__, user->hum0,        &user->hum0,        NULL);CHKERRQ(ierr);
  ngrad = 3;
  ierr = PetscOptionsRealArray("-grad_temp0", "temperature gradient x,y,z [K/m]", __FILE__, user->grad_temp0, &ngrad, &flg);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-eps",      "phase-field interface width [m]", __FILE__, user->eps,       &user->eps,         NULL);CHKERRQ(ierr);

  // Grain input
  ierr = PetscOptionsInt("-readFlag",    "1=read grains from file, 0=generate",  __FILE__, user->readFlag,    &user->readFlag,    NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-grains_file", "path to grains.dat",                __FILE__, user->grains_file, user->grains_file,  sizeof(user->grains_file), NULL);CHKERRQ(ierr);

  // Output
  ierr = PetscOptionsString("-output_dir", "directory for simulation output", __FILE__, user->output_dir, user->output_dir, sizeof(user->output_dir), NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-NASA_output",  "enable sol*.dat output files",    __FILE__, *output,  output,  NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-NASA_monitor", "enable monitor/SSA output",       __FILE__, *monitor, monitor, NULL);CHKERRQ(ierr);

  // Restart / geometry loading
  ierr = PetscOptionsString("-initial_cond",   "load initial solution from file",      __FILE__, user->initial_cond, user->initial_cond, sizeof(user->initial_cond), NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-initial_PFgeom", "load initial ice geometry from file",  __FILE__, user->PFgeom,       user->PFgeom,       sizeof(user->PFgeom),       NULL);CHKERRQ(ierr);

  // Material properties (overridable, defaults set in main before this call)
  ierr = PetscOptionsReal("-thcond_ice", "ice thermal conductivity [W/m/K]",    __FILE__, user->thcond_ice, &user->thcond_ice, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-thcond_air", "air thermal conductivity [W/m/K]",    __FILE__, user->thcond_air, &user->thcond_air, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-thcond_met", "metal thermal conductivity [W/m/K]",  __FILE__, user->thcond_met, &user->thcond_met, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-cp_ice",     "ice heat capacity [J/kg/K]",          __FILE__, user->cp_ice,     &user->cp_ice,     NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-rho_ice",    "ice density [kg/m3]",                 __FILE__, user->rho_ice,    &user->rho_ice,    NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-lat_sub",    "latent heat of sublimation [J/kg]",   __FILE__, user->lat_sub,    &user->lat_sub,    NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-dif_vap",    "vapor diffusivity [m2/s]",            __FILE__, user->dif_vap,    &user->dif_vap,    NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-xi_v",       "vapor time scaling",                  __FILE__, user->xi_v,       &user->xi_v,       NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-xi_T",       "temperature time scaling",            __FILE__, user->xi_T,       &user->xi_T,       NULL);CHKERRQ(ierr);

  PetscOptionsEnd();

  // Propagate parsed scalars into AppCtx fields used by IGA setup
  user->Nx  = (PetscReal)(*Nx);
  user->Ny  = (PetscReal)(*Ny);
  user->Nz  = (PetscReal)(*Nz);
  user->Lx  = *Lx;
  user->Ly  = *Ly;
  user->Lz  = *Lz;
  user->dim = *dim;

  PetscFunctionReturn(0);
}
