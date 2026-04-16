#include "setup_thermal.h"
#include "material.h"

/*-----------------------------------------------------------------------------
  InitializeUserContext
  Set default material properties and discretisation order in AppCtx.
  Must be called before SetupIGA.
-----------------------------------------------------------------------------*/
void InitializeUserContext(AppCtx *user)
{
  user->thcond_ice = 2.29;     /* W/m·K  */
  user->thcond_air = 0.02;     /* W/m·K  */
  user->cp_ice     = 1.96e3;   /* J/kg·K */
  user->cp_air     = 1.044e3;  /* J/kg·K */
  user->rho_ice    = 919.0;    /* kg/m³  */
  user->rho_air    = 1.341;    /* kg/m³  */

  user->p = 2;  /* Quadratic B-splines                    */
  user->C = 1;  /* C^1 inter-element continuity (p-1)    */
}

/*-----------------------------------------------------------------------------
  SetupIGA
  Create and configure the IGA object for the periodic homogenization problem:
    - dim DOFs (one correction-temperature component per spatial direction)
    - Fully periodic in all axes
    - Uniform mesh on [0,Lx] × [0,Ly] [× [0,Lz]]
  Writes igaice.dat to user->output_dir.
-----------------------------------------------------------------------------*/
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga)
{
  PetscErrorCode ierr;
  IGAAxis        axisX, axisY, axisZ;

  PetscFunctionBegin;

  ierr = IGACreate(PETSC_COMM_WORLD, iga); CHKERRQ(ierr);
  ierr = IGASetDim(*iga, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(*iga, user->dim); CHKERRQ(ierr);
  ierr = IGASetFieldName(*iga, 0, "temperature_gradients"); CHKERRQ(ierr);

  /* X axis */
  ierr = IGAGetAxis(*iga, 0, &axisX); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisX, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axisX, PETSC_TRUE); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axisX, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);

  /* Y axis */
  ierr = IGAGetAxis(*iga, 1, &axisY); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisY, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axisY, PETSC_TRUE); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axisY, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);

  /* Z axis (3-D only) */
  if (user->dim == 3) {
    ierr = IGAGetAxis(*iga, 2, &axisZ); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisZ, user->p); CHKERRQ(ierr);
    ierr = IGAAxisSetPeriodic(axisZ, PETSC_TRUE); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisZ, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
  }

  /* Force AIJ storage so AMG (GAMG/HYPRE) can build its graph.
     GAMG/HYPRE do not support creategraph on BAIJ in PETSc 3.20. */
  ierr = IGASetMatType(*iga, MATAIJ); CHKERRQ(ierr);

  ierr = IGASetFromOptions(*iga); CHKERRQ(ierr);
  ierr = IGASetUp(*iga); CHKERRQ(ierr);

  /* Save the IGA descriptor alongside the other output files */
  char filename[PETSC_MAX_PATH_LEN];
  snprintf(filename, sizeof(filename), "%s/igaice.dat", user->output_dir);
  ierr = IGAWrite(*iga, filename); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
