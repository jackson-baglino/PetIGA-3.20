#include "io_thermal.h"
#include "utils.h"

/*--------------------------------------------------------------------------------------------------
  Function: EvaluateFieldAtGaussPoints
  Purpose: Evaluates the scalar field from a PETSc vector at Gauss points and stores it in 
           user->ice array.
  Inputs: 
    - user: Application context containing simulation parameters and output field pointer.
    - iga:  IGA object describing geometry and discretization.
    - vec_phase: Vector containing the solution to be evaluated.
--------------------------------------------------------------------------------------------------*/
static PetscErrorCode EvaluateFieldAtGaussPoints(AppCtx *user, IGA iga, Vec vec_phase) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGAElement            element;
  IGAPoint              point;

  PetscInt              idx = 0;
  Vec                   localU;
  const PetscScalar     *arrayU;
  PetscScalar           *U;

  ierr = IGAGetLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);

  ierr = IGABeginElement(iga, &element); CHKERRQ(ierr);
  while (IGANextElement(iga, element)) {
    ierr = IGAElementGetValues(element, arrayU, &U); CHKERRQ(ierr);
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);

    while (IGAElementNextPoint(element, point)) {
      PetscScalar sol[3];
      IGAPointFormValue(point, U, &sol[0]);
      PetscScalar phi = sol[0];
      user->ice[idx++] = phi;
    }

    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  
  PetscInt nmb_expected = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (idx != nmb_expected) {
    PetscPrintf(PETSC_COMM_WORLD, "⚠️ Warning: Assigned %d Gauss point values, expected %d\n", idx, nmb_expected);
  }

  ierr = IGAEndElement(iga, &element); CHKERRQ(ierr);

  ierr = IGARestoreLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "✅ Gauss point ice field computed with %d values.\n", idx);

  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: ReadSolutionVec
  Purpose: Reads an IGA object and solution vector from file, evaluates the scalar field 
           at Gauss points, and stores it in user->ice.
  Inputs:
    - iga_file: Path to the IGA file (.iga).
    - vec_file: Path to the vector file (.vec).
    - user: Application context (to store loaded field and IGA).
  Outputs:
    - iga_out: Pointer to loaded IGA object.
    - vec_out: Pointer to loaded solution vector (before destruction).
--------------------------------------------------------------------------------------------------*/
PetscErrorCode ReadSolutionVec(const char *iga_file, const char *vec_file, 
                               IGA *iga_out, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGA iga_input;  // IGA object read from file (external to simulation)

  // Create the IGA object (NOTE: WE CAN PROBABLY REPLACE ALL THIS WITH A PREDEFINED FUNCTION)
  ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
  ierr = IGASetDim(iga_input, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(iga_input, 3); CHKERRQ(ierr);

  // Set up the axes
  IGAAxis axis0, axis1, axis2;
  ierr = IGAGetAxis(iga_input, 0, &axis0); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0, user->p); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);

  ierr = IGAGetAxis(iga_input, 1, &axis1); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1, user->p); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);

  if (user->dim == 3) {
    ierr = IGAGetAxis(iga_input, 2, &axis2); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axis2, user->p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis2, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
  }

  // Set the IGA from options and set up
  ierr = IGASetFromOptions(iga_input); CHKERRQ(ierr);
  ierr = IGASetUp(iga_input); CHKERRQ(ierr);

  // Set the input IGA object
  user->iga_input = iga_input;

  // Allocate memory fro user->ice
  ierr = AllocateAppCtxFields(iga_input, user, &user->ice); CHKERRQ(ierr);

  // Create the solution vector
  Vec ice_phase;
  ierr = IGACreateVec(iga_input, &ice_phase); CHKERRQ(ierr);
  ierr = IGAReadVec(iga_input, ice_phase, vec_file); CHKERRQ(ierr);

  // Assign Gauss points via a function (to be created)
  ierr = EvaluateFieldAtGaussPoints(user, iga_input, ice_phase); CHKERRQ(ierr);

  // Clean up
  ierr = VecDestroy(&ice_phase); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}







/*--------------------------------------------------------------------------------------------------
------------------------------------------- SCRATCH BELOW ------------------------------------------
--------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------
  Function: ExtractIceField
  Purpose: Extracts the first DOF field from a multi-component PETSc vector.
  Inputs:
    - U: Global vector with multiple DOFs per point.
    - dof: Number of degrees of freedom per point.
  Outputs:
    - ice_out: Pointer to new Vec containing only the ice component.
--------------------------------------------------------------------------------------------------*/
static PetscErrorCode ExtractIceField(Vec U, PetscInt dof, PetscScalar *ice_out) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt N;
  ierr = VecGetSize(U, &N); CHKERRQ(ierr);
  if (N % dof != 0) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ, "Mismatch between vector size and DOF count");
  }

  PetscInt num_points = N / dof;

  const PetscScalar *U_array;
  ierr = VecGetArrayRead(U, &U_array); CHKERRQ(ierr);
  for (PetscInt i = 0; i < num_points; i++) {
    ice_out[i] = U_array[i * dof];
  }
  ierr = VecRestoreArrayRead(U, &U_array); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: LoadIceField
  Purpose: Loads an IGA and associated solution vector, extracts the ice field (first DOF),
           and stores it in user->ice.
  Inputs:
    - user: Application context containing sim configuration.
    - iga_filename: Path to .iga file.
    - vec_filename: Path to solution vector file.
--------------------------------------------------------------------------------------------------*/
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
  PetscInt dof = 3;
  // ierr = IGAGetDof(iga_input, &dof); CHKERRQ(ierr);
  ierr = ExtractIceField(U, dof, user->ice);

  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// PROBABLY WANT TO DELETE THIS FUNCTION
// /*--------------------------------------------------------------------------------------------------
//   Function: ComputeAndStoreThermalConductivity
//   Purpose: Computes the thermal conductivity field based on the ice field and stores it in Vec K.
//   Inputs:
//     - user: Application context with ice field.
//   Outputs:
//     - K: PETSc Vec to store the thermal conductivity values.
// --------------------------------------------------------------------------------------------------*/
// PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K) {
//   PetscErrorCode ierr;
//   PetscInt i, N = user->Nx * user->Ny;
//   PetscScalar ice, cond, dcond_ice;

//   PetscFunctionBegin;
//   const PetscScalar *ice_array = user->ice;

//   for (i = 0; i < N; i++) {
//     ice = ice_array[i];
//     ThermalCond(user, ice, &cond, &dcond_ice);
//     ierr = VecSetValue(K, i, cond, INSERT_VALUES); CHKERRQ(ierr);
//   }

//   ierr = VecRestoreArrayRead(user->ice, &ice_array); CHKERRQ(ierr);
//   ierr = VecAssemblyBegin(K); CHKERRQ(ierr);
//   ierr = VecAssemblyEnd(K); CHKERRQ(ierr);

//   PetscPrintf(PETSC_COMM_WORLD, "Thermal conductivity field computed and stored.\n");
//   PetscFunctionReturn(0);
// }

/*--------------------------------------------------------------------------------------------------
  Function: WriteOutput
  Purpose: Writes the given PETSc vector to file in binary format if user->outputBinary is set.
  Inputs:
    - user: Application context containing the binary flag.
    - x: Vector to be written.
    - filename: Name of output file.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (user->outputBinary) {
    PetscPrintf(PETSC_COMM_WORLD, "Writing binary output...\n\n");
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: WriteIceFieldToFile
  Purpose: Writes the contents of user->ice to a human-readable .dat text file.
  Inputs:
    - filename: Path to output .dat file.
    - user: Application context containing the ice field.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user) {
  FILE *file;
  PetscFunctionBegin;
  PetscErrorCode ierr;
  IGAElement element;
  IGAPoint point;
  PetscInt indGP;
  PetscReal ice_val;

  file = fopen(filename, "w");
  if (!file) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening file for writing: %s", filename);
  }

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
      ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
      while (IGAElementNextPoint(element, point)) {
        indGP = point->index + point->count * point->parent->index;

        ice_val = user->ice[indGP];
        fprintf(file, "%g %g %g %g\n", point->mapX[0][0], point->mapX[0][1], point->mapX[0][2], ice_val);

        // user->ice[indGP] = 1.0;
      }
      ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  fclose(file);
  PetscPrintf(PETSC_COMM_WORLD, "Ice field successfully written to file: %s\n", filename);
  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: WriteBinaryFile
  Purpose: Writes a PETSc vector to a binary file using a PetscViewer.
  Inputs:
    - field: PETSc vector to write.
    - filename: Output file name.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteBinaryFile(Vec field, const char *filename) {
  PetscErrorCode ierr;
  PetscViewer viewer;
  PetscBool fileExists;

  PetscFunctionBegin;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

  ierr = VecView(field, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscTestFile(filename, 'r', &fileExists); CHKERRQ(ierr);
  if (!fileExists) {
    PetscPrintf(PETSC_COMM_WORLD, "Warning: Binary output file %s was not created!\n", filename);
  }
  
  PetscFunctionReturn(0);
}