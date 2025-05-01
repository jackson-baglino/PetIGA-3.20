#include "assembly.h"
#include "material_properties.h"

// Function declarations
static PetscErrorCode GetIceAtGaussPoint(IGAPoint pnt, AppCtx *user, PetscScalar *ice);

/*------------------------------------------------------------------------------
  Function: AssembleStiffnessMatrix
  Purpose: Assemble FEM stiffness matrix and apply bottom flux BC.
------------------------------------------------------------------------------*/
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Retrieve user context and defiene variables
  AppCtx *user = (AppCtx *)ctx;
  PetscInt a, b, l;
  PetscInt nen = pnt->nen, dim = user->dim;
  const PetscScalar *N0;
  const PetscReal (*N1)[dim];

  // Retrieve basis functions
  ierr = IGAPointGetShapeFuns(pnt, 0, &N0); CHKERRQ(ierr);
  ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

  // Get ice phase at Gauss point
  PetscScalar ice;
  ierr = GetIceAtGaussPoint(pnt, user, &ice); CHKERRQ(ierr);

  // Evaluate thermal conductivity
  PetscReal thcond;
  ThermalCond(user, ice, &thcond, NULL);

  // Get the weight of the Gauss point
  PetscReal W = *pnt->weight;

  // Initialize stiffness matrix and force vector
  for (a = 0; a < nen; a++) {
    for (b = 0; b < nen; b++) {
      for (l = 0; l < dim; l++) {
        // Assemble stiffness matrix
        K[a * nen + b] += thcond * N1[a][l] * N1[b][l] * W;
      }
    }
    // Apply Neumann boundary condition for bottom flux
    if (pnt->boundary_id == 2) {
      F[a] += N0[a] * user->q_bottom * W;
    }
  }
    
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: ApplyBoundaryConditions
  Purpose: Set top temperature, bottom flux, and side zero-flux BCs.
------------------------------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;
  
    PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");
  
    /* Fixed temperature at the top */
    ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n\n", user->T_top);

    // /* Fixed temperature at the top */
    // ierr = IGASetBoundaryValue(iga, 1, 0, 0, user->T_top-3.0); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n\n", user->T_top);

    /* Prescribed flux at the bottom */
    ierr = IGASetBoundaryForm(iga, 1, 0, PETSC_TRUE); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Prescribed flux at bottom: q = %g W/m²\n", user->q_bottom);

    // /* Insulated side boundaries */
    // ierr = IGASetBoundaryForm(iga, 0, 0, PETSC_TRUE); CHKERRQ(ierr);
    // ierr = IGASetBoundaryForm(iga, 0, 1, PETSC_TRUE); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "  - Insulated side boundaries applied\n");
  
    PetscFunctionReturn(0);
  }

/*------------------------------------------------------------------------------
  Function: InitialCondition
  Purpose: Set initial temperature using linear profile from BCs.
------------------------------------------------------------------------------*/
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscScalar *T_array;
  ierr = VecGetArray(T, &T_array); CHKERRQ(ierr);

  PetscInt Nx = user->Nx, Ny = user->Ny;
  PetscReal Ly = user->Ly;
  PetscReal T_top = user->T_top;
  PetscReal q_bottom = user->q_bottom;
  PetscReal k_eff = 0.5 * (user->thcond_ice + user->thcond_air);

  for (PetscInt j = 0; j <= Ny; j++) {
    PetscReal y = (PetscReal)j * Ly / (PetscReal)Ny;
    PetscReal T_init = T_top - (q_bottom / k_eff) * (Ly - y);
    for (PetscInt i = 0; i <= Nx; i++) {
      PetscInt idx = j * (Nx + 1) + i;
      T_array[idx] = T_init;
    }
  }

  ierr = VecRestoreArray(T, &T_array); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: GetIceAtGaussPoint -- CONSIDERING MOVING THIS TO A SCRIPT OF HELPER FUNCTIONS
  Helper: retrieve ice field value at Gauss point index.
------------------------------------------------------------------------------*/
static PetscErrorCode GetIceAtGaussPoint(IGAPoint pnt, AppCtx *user, PetscScalar *ice) {\

  PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
  *ice = user->ice[indGP];

  return 0;
}

/*------------------------------------------------------------------------------
  Function: SetupIGA
  Purpose:  Sets up the IGA object by defining the problem domain, degrees of 
            freedom, field names, and the uniform grid.
------------------------------------------------------------------------------*/
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = IGACreate(PETSC_COMM_WORLD, iga); CHKERRQ(ierr);
  ierr = IGASetDim(*iga, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(*iga, 1); CHKERRQ(ierr);
  ierr = IGASetFieldName(*iga, 0, "temperature"); CHKERRQ(ierr);

  /* Set up uniform grid along each axis */
  IGAAxis axisX, axisY, axisZ;
  ierr = IGAGetAxis(*iga, 0, &axisX); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisX, user->p); CHKERRQ(ierr);
  ierr = IGAGetAxis(*iga, 1, &axisY); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisY, user->p); CHKERRQ(ierr);
  if (user->dim == 3) {
      ierr = IGAGetAxis(*iga, 2, &axisZ); CHKERRQ(ierr);
  }
  
  /* Initialize each axis */
  ierr = IGAAxisInitUniform(axisX, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axisY, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);
  if (user->dim == 3) {
      ierr = IGAAxisInitUniform(axisZ, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
  }

  /* Finalize IGA setup */
  ierr = IGASetFromOptions(*iga); CHKERRQ(ierr);
  ierr = IGASetUp(*iga); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "IGA setup complete.\n\n");
  PetscFunctionReturn(0);
}