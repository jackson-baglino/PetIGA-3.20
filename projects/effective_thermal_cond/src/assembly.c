#include "assembly.h"
#include "material_properties.h"

// Function declarations
static PetscErrorCode GetIceAtGaussPoint(IGAPoint pnt, AppCtx *user, PetscScalar *ice);

/*------------------------------------------------------------------------------
  Function: AssembleStiffnessMatrix
Assemble FEM stiffness matrix and apply bottom flux BC.
------------------------------------------------------------------------------*/
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

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

  // Apply bottom Neumann flux BC (boundary_id 2 = y=0)
  if (pnt->atboundary && pnt->boundary_id == 2) {
    for (a = 0; a < nen; a++) {
      F[a] -= N0[a] * user->q_bottom / user->Lx;
    }
  } else {
    for (a = 0; a < nen; a++) {
      for (b = 0; b < nen; b++) {
          for (l = 0; l < dim; l++) {
            K[a * nen + b] += thcond * N1[a][l] * N1[b][l];
          }
      }
    }
  }

  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
Function: ComputeInitialCondition
// Set initial temperature using linear profile from BCs.
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
// Function: GetIceAtGaussPoint -- CONSIDERING MOVING THIS TO A SCRIPT OF HELPER FUNCTIONS
// Helper: retrieve ice field value at Gauss point index.
------------------------------------------------------------------------------*/
static PetscErrorCode GetIceAtGaussPoint(IGAPoint pnt, AppCtx *user, PetscScalar *ice) {
  PetscErrorCode ierr;
  const PetscScalar *iceArray;
  ierr = VecGetArrayRead(user->ice, &iceArray); CHKERRQ(ierr);
  PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
  *ice = iceArray[indGP];
  ierr = VecRestoreArrayRead(user->ice, &iceArray); CHKERRQ(ierr);
  return 0;
}

/*-----------------------------------------------------------
  Function: ApplyBoundaryConditions
  Set top temperature, bottom flux, and side zero-flux BCs.
-----------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");

  /* Fixed temperature at the top */
  ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n", user->T_top);

  /* Prescribed flux at the bottom */
  ierr = IGASetBoundaryForm(iga, 1, 0, PETSC_TRUE); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "  - Prescribed flux at bottom: q = %g W/m²\n", user->q_bottom);

  /* Zero-flux (Neumann) conditions on the z-axis (if 3D) */
  if (user->dim == 3) {
      ierr = IGASetBoundaryForm(iga, 2, 0, PETSC_TRUE); CHKERRQ(ierr);
      ierr = IGASetBoundaryForm(iga, 2, 1, PETSC_TRUE); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}