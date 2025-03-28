#include "thermal_solver.h"

PetscErrorCode LoadIceField(AppCtx *user, const char *iga_filename, const char *vec_filename) {
    PetscErrorCode ierr;
    IGA iga_input;  // IGA object from the saved file (NOT the simulation's IGA)
    Vec U;
    PetscScalar *array;

    PetscFunctionBegin;

    // Read the IGA structure from the saved file
    ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
    ierr = IGARead(iga_input, iga_filename); CHKERRQ(ierr);

    // Read the solution vector associated with this IGA
    ierr = IGACreateVec(iga_input, &U); CHKERRQ(ierr);
    ierr = IGAReadVec(iga_input, U, vec_filename); CHKERRQ(ierr);

    // Access array from PETSc Vec
    ierr = VecGetArray(U, &array); CHKERRQ(ierr);

    // Extract metadata from the loaded IGA
    PetscInt N, dof;
    ierr = VecGetSize(U, &N); CHKERRQ(ierr);
    ierr = IGAGetDof(iga_input, &dof); CHKERRQ(ierr);

    if (N % dof != 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Mismatch between vector size and DOF count");
    }

    PetscInt num_points = N / dof;

    // Assign to user->ice (assuming ice phase is stored in DOF index 0)
    for (PetscInt i = 0; i < num_points; i++) {
        user->ice[i] = array[i * dof]; // Extracting first DOF (phase ice field)
    }

    // Restore and clean up
    ierr = VecRestoreArray(U, &array); CHKERRQ(ierr);
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


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

PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    if (user->outputBinary) {
        PetscPrintf(PETSC_COMM_WORLD, "Writing binary...\n");

        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
        ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#include <petscsys.h>
#include <petscviewer.h>

// Function to write the ice field to a .dat file
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user) {
    FILE *file;
    PetscInt num_points = (user->Nx) * (user->Ny) * (user->dim == 3 ?user->Nz : 1);
    
    PetscFunctionBegin; // PETSc standard error handling start
    
    // Open file for writing
    file = fopen(filename, "w");
    if (!file) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening file for writing: %s", filename);
    }

    // Write each value to the file, one per line
    for (PetscInt i = 0; i < num_points; i++) {
        if (fprintf(file, "%.16e\n", user->ice[i]) < 0) {
            fclose(file);
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_WRITE, "Error writing to ice field file.");
        }
    }

    // Close the file
    fclose(file);

    PetscPrintf(PETSC_COMM_WORLD, "✅ Ice field successfully written to file: %s\n", filename);
    
    PetscFunctionReturn(0); // PETSc standard return
}
