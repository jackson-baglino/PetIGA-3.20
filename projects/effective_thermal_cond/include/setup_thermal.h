#ifndef SETUP_H
#define SETUP_H

#include "app_ctx.h"

/*
 * InitializeUserContext — fill AppCtx with default material properties
 * (thcond_ice, thcond_air, cp_*, rho_*) and discretization order (p, C).
 * Must be called before SetupIGA.
 */
void InitializeUserContext(AppCtx *user);

/*
 * SetupIGA — create and configure the IGA object for the homogenization
 * problem: dim DOFs (one per direction), fully periodic in all axes,
 * uniform mesh of size Nx × Ny [× Nz] over [0,Lx] × [0,Ly] [× [0,Lz]].
 * Writes igaice.dat to output_dir.
 */
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga);

#endif /* SETUP_H */
