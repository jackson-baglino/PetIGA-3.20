#ifndef THERMAL_SOLVER_H
#define THERMAL_SOLVER_H

#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"
#include <math.h>
#include <mpi.h>

#include <petscsys.h>
#include <petscviewer.h>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

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
    Vec ice;             // Ice phase variable
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

// ========================= Function Declarations =========================
// ** Environment Setup and Initialization **
PetscErrorCode GetEnvironment(AppCtx *user);
PetscErrorCode FormInitialCondition(AppCtx *user);
PetscErrorCode InitializeFields(AppCtx *user, IGA iga);
void InitializeUserContext(AppCtx *user);

// ** Boundary Conditions and Solver Setup **
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user);
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user);
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx);

// ** Material Properties **
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, PetscScalar *dcond_ice);
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, PetscScalar *dcp_ice);
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, PetscScalar *drho_ice);

// ** IGA Setup and Solver Execution **
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga);
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga);

// ** File I/O and Data Storage **
PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K);
PetscErrorCode WriteBinaryFile(Vec field, const char *filename);
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename);
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user);

#endif /* THERMAL_SOLVER_H */