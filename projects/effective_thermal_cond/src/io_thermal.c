#include "io_thermal.h"

/* 
Function: ExtractIceField
   Extracts the ice field (first DOF) from a PETSc vector and creates a new Vec.
*/
static PetscErrorCode ExtractIceField(Vec U, PetscInt dof, Vec *ice_out) {
  PetscErrorCode ierr;
  PetscInt N;
  ierr = VecGetSize(U, &N); CHKERRQ(ierr);
  if (N % dof != 0) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Mismatch between vector size and DOF count");
  }
  PetscInt num_points = N / dof;
  Vec ice;
  ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
  ierr = VecSetSizes(ice, num_points, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = VecSetFromOptions(ice); CHKERRQ(ierr);
  
  PetscScalar *ice_array;
  ierr = VecGetArray(ice, &ice_array); CHKERRQ(ierr);
  const PetscScalar *U_array;
  ierr = VecGetArrayRead(U, &U_array); CHKERRQ(ierr);
  for (PetscInt i = 0; i < num_points; i++) {
    ice_array[i] = U_array[i * dof];
  }
  ierr = VecRestoreArrayRead(U, &U_array); CHKERRQ(ierr);
  ierr = VecRestoreArray(ice, &ice_array); CHKERRQ(ierr);
  
  *ice_out = ice;
  return 0;
}

/* Function: LoadIceField
   Loads a saved IGA and its corresponding solution vector, then extracts the ice field.
*/
PetscErrorCode LoadIceField(AppCtx *user, const char *iga_filename, const char *vec_filename) {
  PetscErrorCode ierr;
  IGA iga_input;  // IGA object read from file (external to simulation)
  Vec U;

  PetscFunctionBegin;

  ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
  ierr = IGARead(iga_input, iga_filename); CHKERRQ(ierr);

  ierr = IGACreateVec(iga_input, &U); CHKERRQ(ierr);
  ierr = IGAReadVec(iga_input, U, vec_filename); CHKERRQ(ierr);

  // Get DOF from the loaded IGA and extract the ice field
  PetscInt dof;
  ierr = IGAGetDof(iga_input, &dof); CHKERRQ(ierr);
  ierr = ExtractIceField(U, dof, &user->ice); CHKERRQ(ierr);

  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* Function: ComputeAndStoreThermalConductivity
   Computes thermal conductivity based on the ice field and stores it in Vec K.
*/
PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K) {
  PetscErrorCode ierr;
  PetscInt i, N = user->Nx * user->Ny;
  PetscScalar ice, cond, dcond_ice;

  PetscFunctionBegin;
  const PetscScalar *ice_array;
  ierr = VecGetArrayRead(user->ice, &ice_array); CHKERRQ(ierr);

  for (i = 0; i < N; i++) {
    ice = ice_array[i];
    ThermalCond(user, ice, &cond, &dcond_ice);
    ierr = VecSetValue(K, i, cond, INSERT_VALUES); CHKERRQ(ierr);
  }

  ierr = VecRestoreArrayRead(user->ice, &ice_array); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(K); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(K); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "Thermal conductivity field computed and stored.\n");
  PetscFunctionReturn(0);
}

/* Function: WriteOutput
   Writes the given PETSc vector to a file (binary if enabled).
*/
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (user->outputBinary) {
    PetscPrintf(PETSC_COMM_WORLD, "Writing binary output...\n");
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* Function: WriteIceFieldToFile
   Writes the ice field from user->ice to a text (.dat) file.
*/
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user) {
    FILE *file;
    PetscInt num_points = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
    PetscErrorCode ierr;

    PetscFunctionBegin;
    file = fopen(filename, "w");
    if (!file) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening file for writing: %s", filename);
    }

    const PetscScalar *ice_array;
    ierr = VecGetArrayRead(user->ice, &ice_array); CHKERRQ(ierr);

    for (PetscInt i = 0; i < num_points; i++) {
      if (fprintf(file, "%.16e\n", ice_array[i]) < 0) {
        fclose(file);
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_WRITE, "Error writing to ice field file.");
      }
    }

    ierr = VecRestoreArrayRead(user->ice, &ice_array); CHKERRQ(ierr);
    fclose(file);
    PetscPrintf(PETSC_COMM_WORLD, "Ice field successfully written to file: %s\n", filename);
    PetscFunctionReturn(0);
}

/* Function: WriteBinaryFile
   Writes a PETSc vector field to a binary file using a PetscViewer.
*/
PetscErrorCode WriteBinaryFile(Vec field, const char *filename) {
  PetscErrorCode ierr;
  PetscViewer viewer;
  PetscBool fileExists;

  PetscFunctionBegin;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

  PetscScalar norm;
  ierr = VecNorm(field, NORM_2, &norm); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Norm of vector before writing: %g\n", (double)norm);

  ierr = VecView(field, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscTestFile(filename, 'r', &fileExists); CHKERRQ(ierr);
  if (!fileExists) {
    PetscPrintf(PETSC_COMM_WORLD, "Warning: Binary output file %s was not created!\n", filename);
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "Binary output successfully written to %s\n", filename);
  }
  PetscFunctionReturn(0);
}