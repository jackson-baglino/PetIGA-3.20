#ifndef OPTIONS_HELPER_H
#define OPTIONS_HELPER_H

#include "NASA_types.h"

PetscErrorCode GetOptions(AppCtx *user,
                          PetscInt *Nx, PetscInt *Ny, PetscInt *Nz,
                          PetscReal *Lx, PetscReal *Ly, PetscReal *Lz,
                          PetscReal *delt_t, PetscReal *t_final, PetscInt *n_out,
                          PetscInt *dim, PetscInt *p, PetscInt *C,
                          PetscBool *output, PetscBool *monitor);

#endif // OPTIONS_HELPER_H
