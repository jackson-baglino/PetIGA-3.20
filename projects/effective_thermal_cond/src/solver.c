#include "solver.h"
#include "material_properties.h"

/**
 * @brief Assembles the stiffness matrix for the finite element method.
 * 
 * This function computes the element stiffness matrix based on the ice phase field and
 * thermal conductivity. It also applies a Neumann boundary condition (flux) at the bottom boundary.
 *
 * @param pnt        IGAPoint representing the Gauss point in the element.
 * @param K          Pointer to the stiffness matrix to be filled.
 * @param F          Pointer to the force vector to be filled.
 * @param ctx        Context containing user-defined parameters.
 * @return PetscErrorCode Error code for PETSc error handling.
 */
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // Retrieve user context
    AppCtx *user = (AppCtx *)ctx;

    PetscInt a, b, l;
    PetscInt nen = pnt->nen, dim = user->dim;
    const PetscScalar *N0, (*N1)[dim];

    // Get the ice phase value from the PETSc Vec
    PetscScalar ice;
    {
      const PetscScalar *iceArray;
      ierr = VecGetArrayRead(user->ice, &iceArray); CHKERRQ(ierr);
      /* 
         Compute an index for the Gauss point. Note: If the ice field is stored on the nodal grid,
         then a proper projection or interpolation might be necessary. Here we assume that the nodal 
         value corresponding to the Gauss point index is acceptable.
      */
      PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
      ice = iceArray[indGP];
      ierr = VecRestoreArrayRead(user->ice, &iceArray); CHKERRQ(ierr);
    }

    // Compute thermal conductivity based on the ice phase
    PetscReal thcond;
    ThermalCond(user, ice, &thcond, NULL);

    // Get basis functions
    ierr = IGAPointGetShapeFuns(pnt, 0, &N0); CHKERRQ(ierr);
    ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

    /* ===================================================
       APPLY FLUX BC AT THE BOTTOM BOUNDARY (y=0)
       =================================================== */
    if (pnt->atboundary && pnt->boundary_id == 2) {  // Check if the Gauss point is on a boundary
        // Boundary in y-direction (dir = 1) and bottom side (side = 0)
        for (a = 0; a < nen; a++) {
            F[a] -= N0[a] * user->q_bottom / user->Lx; // Apply flux as a Neumann condition
        }
    } else {
        // Loop over test functions (a) and trial functions (b)
        for (a = 0; a < nen; a++) {
            for (b = 0; b < nen; b++) {
                for (l = 0; l < dim; l++) {
                    K[a * nen + b] += thcond * N1[a][l] * N1[b][l]; // Diffusion term
                }
            }
        }
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Computes the initial temperature distribution.
 * 
 * This function initializes the temperature field with a linear profile based on a prescribed
 * top temperature and a heat flux condition at the bottom. It assumes an effective thermal conductivity
 * computed as the average of ice and air conductivity.
 *
 * @param T          PETSc vector for storing the temperature field.
 * @param user       Pointer to the user-defined application context containing parameters.
 * @return PetscErrorCode Error code for PETSc error handling.
 */
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // Get array to modify Vec directly
    PetscScalar *T_array;
    ierr = VecGetArray(T, &T_array); CHKERRQ(ierr);

    // Get mesh dimensions
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    PetscReal Ly = user->Ly;

    // Get boundary conditions
    PetscReal T_top = user->T_top;
    PetscReal q_bottom = user->q_bottom;  // Prescribed heat flux at the bottom

    // Approximate an effective thermal conductivity (use average of min and max)
    PetscReal k_eff = 0.5 * (user->thcond_ice + user->thcond_air);  // Can be refined

    // Loop over the mesh points and set initial temperature
    for (PetscInt j = 0; j <= Ny; j++) {
        PetscReal y = (PetscReal)j * Ly / (PetscReal)Ny;  // Physical y-position
        PetscReal T_init = T_top - (q_bottom / k_eff) * (Ly - y); // Linear profile

        for (PetscInt i = 0; i <= Nx; i++) {
            PetscInt idx = j * (Nx + 1) + i;  // Flattened 2D index
            T_array[idx] = T_init;  // Assign temperature
        }
    }

    // Restore array back to PETSc Vec
    ierr = VecRestoreArray(T, &T_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}