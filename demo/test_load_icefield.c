#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"
#include <math.h>
#include <mpi.h>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

#define TEST_ICE_VALUE 0.8  // Example ice phase value for testing

// Structure to hold application-specific parameters and data
typedef struct {
  IGA iga;  // Isogeometric Analysis (IGA) object for managing the finite element discretization

  // **Material properties**
  PetscReal thcond_ice;   // Thermal conductivity of ice (W/mK)
  PetscReal thcond_air;   // Thermal conductivity of air (W/mK)
  PetscReal cp_ice;       // Specific heat capacity of ice (J/kgK)
  PetscReal cp_air;      // Specific heat capacity of air (J/kgK)
  PetscReal rho_ice;     // Density of ice (kg/m³)
  PetscReal rho_air;     // Density of air (kg/m³)

  // **Initial conditions**
  PetscReal temp0;            // Initial temperature (K)
  PetscReal grad_temp0[3];    // Initial temperature gradient (x, y, z components)
  PetscReal *ice;             // Ice phase variable
  Vec T_sol;                  // Solution vector for temperature


  // **Domain size and mesh resolution**
  PetscInt dim;         // 2D or 3D
  PetscReal Lx, Ly, Lz;  // Domain size
  PetscInt Nx, Ny, Nz;  // Mesh resolution
  PetscReal eps;        // Interface width parameter for phase field method
  PetscInt  p;         // Polynomial orders for basis functions
  PetscInt  C;         // Global continuity order

  // **Dirichlet BC (Fixed Temperature at Top & Bottom)**
  PetscReal T_top;    // Temperature at y = Ly

  // **Neumann BC (Constant flux at top or bottom)**
  PetscReal q_bottom; // Heat flux at y = 0
  PetscBool useFluxBottom; // If true, apply flux at bottom boundary

  // **Input options**
  char init_mode[256]; // Mode for initializing the ice field 

  // **Output options**
  PetscBool outputBinary; // Flag for binary output

} AppCtx;

PetscErrorCode LoadIceField(AppCtx *user, const char *iga_filename, const char *vec_filename) {
  PetscErrorCode ierr;
  IGA iga_input;  // IGA object from the saved file (NOT the simulation's IGA)
  Vec U;
  PetscScalar *array;

  PetscFunctionBegin;

  // Read the IGA structure from the saved file
  ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
  ierr = IGARead(iga_input, iga_filename); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Loaded IGA from %s\n", iga_filename);

  // Read the solution vector associated with this IGA
  ierr = IGACreateVec(iga_input, &U); CHKERRQ(ierr);
  ierr = IGAReadVec(iga_input, U, vec_filename); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Loaded vector from %s\n", vec_filename);

  // Access array from PETSc Vec
  ierr = VecGetArray(U, &array); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Accessed array from vector\n");

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

// Function to create a dummy IGA and solution vector
PetscErrorCode CreateTestIGA(const char *iga_filename, const char *vec_filename) {
  PetscErrorCode ierr;
  IGA iga_test;
  IGAAxis axisX, axisY;
  Vec U;
  PetscScalar *array;
  PetscInt num_points = 10;  // Small test grid
  PetscInt dof = 3;  // Assume same DOF count as real case
  PetscInt p = 2;    // Degree of basis functions

  PetscFunctionBegin;

  // Create and set up IGA
  ierr = IGACreate(PETSC_COMM_WORLD, &iga_test); CHKERRQ(ierr);
  ierr = IGASetDim(iga_test, 2); CHKERRQ(ierr);
  ierr = IGASetDof(iga_test, dof); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga_test, 0, "phaseice"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga_test, 1, "temperature"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga_test, 2, "vap_density"); CHKERRQ(ierr);

  // Get axes and set degrees
  ierr = IGAGetAxis(iga_test, 0, &axisX); CHKERRQ(ierr);
  ierr = IGAGetAxis(iga_test, 1, &axisY); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisX, p); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisY, p); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "Setting up IGA with degree %D\n", p);

  // Ensure the IGA is properly initialized
  ierr = IGASetUp(iga_test); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Test IGA setup complete\n");

  ierr = IGAWrite(iga_test, iga_filename); CHKERRQ(ierr);

  // Create a test vector with DOF
  ierr = IGACreateVec(iga_test, &U); CHKERRQ(ierr);
  ierr = VecGetArray(U, &array); CHKERRQ(ierr);

  // Assign test values (ice phase stored in DOF 0)
  for (PetscInt i = 0; i < num_points; i++) {
      array[i * dof] = TEST_ICE_VALUE; // Set ice phase
  }

  ierr = VecRestoreArray(U, &array); CHKERRQ(ierr);
  ierr = IGAWriteVec(iga_test, U, vec_filename); CHKERRQ(ierr);

  // Clean up
  ierr = VecDestroy(&U); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_test); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


// Main test function
int main(int argc, char **argv) {
    PetscErrorCode ierr;
    AppCtx user;
    PetscInt num_points = 10;

    const char *iga_pathname = "/Users/jacksonbaglino/SimulationResults/DrySed_Metamorphism/NASAv2/NASAv2_10G_2D_T-20.0_hum0.70_2025-03-13__14.20.59";
    char iga_filename[512];
    char vec_filename[512];
    snprintf(iga_filename, sizeof(iga_filename), "%s/igasol.dat", iga_pathname);
    snprintf(vec_filename, sizeof(vec_filename), "%s/sol_00440.dat", iga_pathname);

    // Initialize PETSc
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;

    PetscPrintf(PETSC_COMM_WORLD, "iga_filename: %s\n", iga_filename);
    PetscPrintf(PETSC_COMM_WORLD, "vec_filename: %s\n", vec_filename);

    // Allocate space for ice array in user struct
    user.ice = (PetscScalar *)malloc(num_points * sizeof(PetscScalar));
    if (!user.ice) {
        PetscPrintf(PETSC_COMM_WORLD, "Memory allocation failed for user.ice\n");
        return -1;
    }

    // Create test files
    ierr = CreateTestIGA(iga_filename, vec_filename); CHKERRQ(ierr);

    // Call function to test
    ierr = LoadIceField(&user, iga_filename, vec_filename); CHKERRQ(ierr);

    // Validate results
    PetscBool success = PETSC_TRUE;
    for (PetscInt i = 0; i < num_points; i++) {
        if (PetscAbsScalar(user.ice[i] - TEST_ICE_VALUE) > 1e-6) {
            success = PETSC_FALSE;
            PetscPrintf(PETSC_COMM_WORLD, "Test failed at index %D: Expected %g, got %g\n",
                        i, (double)TEST_ICE_VALUE, (double)user.ice[i]);
        }
    }

    if (success) {
        PetscPrintf(PETSC_COMM_WORLD, "✅ Test passed! Ice field values match expected results.\n");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "❌ Test failed. Check printed mismatches.\n");
    }

    // Clean up
    free(user.ice);
    ierr = PetscFinalize();
    return ierr;
}
