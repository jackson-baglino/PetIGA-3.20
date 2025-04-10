#ifndef USER_CONTEXT_H
#define USER_CONTEXT_H

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
    IGA iga_input;  // IGA object for reading input data

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

#endif // USER_CONTEXT_H