#include "thermal_solver.h"

int main (int argc, char *argv[]) {
    PetscErrorCode ierr;
    AppCtx user;

    // /* Initialize PETSc */
    // ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

    // /* Read simulation parameters */
    // ierr = GetEnvironment(&user); CHKERRQ(ierr);

    // /* Setup simulation */
    // ierr = SetupIGA(&user); CHKERRQ(ierr);
    // ierr = ApplyBoundaryConditions(user.iga, &user); CHKERRQ(ierr);

    // /* Solve the system */
    // ierr = SolveThermalDiffusion(&user); CHKERRQ(ierr);

    // /* Finalize PETSc */
    // ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}