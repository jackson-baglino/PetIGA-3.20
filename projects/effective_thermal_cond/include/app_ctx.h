#ifndef APP_CTX_H
#define APP_CTX_H

#include "petiga.h"
#include <petsctime.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <regex.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

typedef struct {
  IGA iga;        /* IGA object for the homogenization problem          */
  IGA iga_input;  /* IGA object used to read the input phase-field data */

  /* Material properties */
  PetscReal thcond_ice;  /* Thermal conductivity of ice  (W/m·K) */
  PetscReal thcond_air;  /* Thermal conductivity of air  (W/m·K) */
  PetscReal cp_ice;      /* Specific heat capacity ice   (J/kg·K) */
  PetscReal cp_air;      /* Specific heat capacity air   (J/kg·K) */
  PetscReal rho_ice;     /* Density of ice               (kg/m³)  */
  PetscReal rho_air;     /* Density of air               (kg/m³)  */

  /* Domain and mesh */
  PetscInt  dim;           /* Problem dimension: 2 or 3              */
  PetscReal Lx, Ly, Lz;   /* Domain size in each direction (m)      */
  PetscInt  Nx, Ny, Nz;   /* Number of elements in each direction   */
  PetscReal eps;           /* Phase-field interface width (m)        */
  PetscInt  p;             /* Polynomial order for basis functions   */
  PetscInt  C;             /* Inter-element continuity order         */

  /* Phase field */
  PetscReal *ice;          /* Ice fraction at every Gauss point      */

  /* Temperature solution vector */
  Vec T_sol;               /* Solution of the homogenization cell problem */

  /* Boundary conditions */
  PetscReal T_top;         /* Top boundary temperature (K)            */
  PetscReal q_bottom;      /* Bottom boundary heat flux (W/m²)        */

  /* Initial condition / input options */
  char init_mode[256];     /* "circle" | "layered" | "file"           */
  char init_dir[256];      /* Directory containing sol_*.dat files    */

  /* Output options */
  PetscBool   outputBinary; /* Write binary t_vec.dat output           */
  PetscInt    sol_index;    /* -1 = loop all files; ≥1 = single index  */
  char        output_dir[PETSC_MAX_PATH_LEN]; /* Output directory path */

} AppCtx;

#endif /* APP_CTX_H */
