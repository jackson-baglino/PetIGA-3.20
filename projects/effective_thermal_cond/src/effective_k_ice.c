#include "user_context.h"
#include "setup_thermal.h"
#include "assembly.h"
#include "io_thermal.h"
#include "material_properties.h"
#include "env_config.h"

int main (int argc, char *argv[]) {
    AppCtx              user;
    PetscErrorCode      ierr;

    /* ------------------ Initialize PETSc ------------------ */
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

    /* ------------------ Set options ------------------ */
    // Look into all of these options. See if there are other options that can be set.
    // These options are set in the command line when running the program.
    PetscBool print_error = PETSC_TRUE;
    PetscBool check_error = PETSC_FALSE;
    PetscBool save = PETSC_FALSE;
    PetscBool draw = PETSC_FALSE;
    PetscOptionsBegin(PETSC_COMM_WORLD,"","Laplace Options","IGA");CHKERRQ(ierr);
    ierr = PetscOptionsBool("-print_error","Prints the error of the solution",__FILE__,print_error,&print_error,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-check_error","Checks the error of the solution",__FILE__,check_error,&check_error,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-save","Save the solution to file",__FILE__,save,&save,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsBool("-draw","If dim <= 2, then draw the solution to the screen",__FILE__,draw,&draw,NULL);CHKERRQ(ierr);
    ierr = PetscOptionsString("-init_mode", "Set initial ice field mode (circle, layered, file)", __FILE__, "circle", user.init_mode, sizeof(user.init_mode), NULL); CHKERRQ(ierr);
    PetscOptionsEnd();CHKERRQ(ierr);

    /* ------------------ Read Environment Variables ------------------ */
    ierr = GetEnvironment(&user); CHKERRQ(ierr);

    /* ------------------ Define user context ------------------ */
    InitializeUserContext(&user);

    /* ------------------ Initialize IGA ------------------ */
    IGA iga;
    ierr = SetupIGA(&user, &iga); CHKERRQ(ierr); // Create and set up the IGA object
    user.iga = iga;

    /* ------------------ Initialize Field Variables ------------------ */
    ierr = InitializeFields(&user, iga); CHKERRQ(ierr); // Initialize the ice field

    /* ------------------ Define Boundary Conditions ------------------ */
    // Apply the boundary conditions
    ierr = ApplyBoundaryConditions(iga, &user); CHKERRQ(ierr);

    /* ------------------ Set Up KSP Solver ------------------ */
    // Creat KSP solver
    ierr = SetupAndSolve(&user, iga); CHKERRQ(ierr);

    /* ------------------ Write Output ------------------ */
    ierr = WriteOutput(&user, user.T_sol, "temperature.bin"); CHKERRQ(ierr); // Write the solution to file
    ierr = WriteIceFieldToFile("ice_field.dat", &user); CHKERRQ(ierr); // Write the ice field to a .dat file

    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.ice); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}