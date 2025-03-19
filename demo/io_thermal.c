#include "thermal_solver.h"

PetscErrorCode WriteBinaryFile(Vec field, const char *filename) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscBool fileExists;
    
    PetscFunctionBegin;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

    // Debug: Check if the vector is valid before writing
    PetscScalar norm;
    ierr = VecNorm(field, NORM_2, &norm);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Norm of vector %s before writing: %g\n", filename, (double)norm);

    ierr = VecView(field, viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    ierr = PetscTestFile(filename, 'r', &fileExists);CHKERRQ(ierr);
    if (!fileExists) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: Binary output file %s was not created!\n", filename);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Binary output successfully written to %s\n", filename);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K) {
    PetscInt i, N;
    PetscScalar ice, cond, dcond_ice;
    PetscErrorCode ierr;

    N = user->Nx * user->Ny; // Total number of grid points

    for (i = 0; i < N; i++) {
        ice = user->ice[i];
        ThermalCond(user, ice, &cond, &dcond_ice);

        // Store thermal conductivity value in K
        ierr = VecSetValue(K, i, cond, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(K); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(K); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Thermal conductivity field computed and stored.\n");

    PetscFunctionReturn(0);
}

PetscErrorCode WriteOutput(AppCtx *user, Vec x) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (user->outputBinary) {
        PetscPrintf(PETSC_COMM_WORLD, "Writing binary...\n");

        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
        ierr = WriteBinaryFile(x, "temperature.bin"); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
