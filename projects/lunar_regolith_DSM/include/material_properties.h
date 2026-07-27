#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H

#include "NASA_types.h"

/* Computes effective thermal conductivity and its derivative with respect to ice */
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, PetscScalar *dcond_ice);

/* Computes effective heat capacity and its derivative with respect to ice */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, PetscScalar *dcp_ice);

/* Computes effective density and its derivative with respect to ice */
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, PetscScalar *drho_ice);

/* Computes vapor diffusivity and its temperature derivative */
void VaporDiffus(AppCtx *user, PetscScalar tem, PetscScalar *difvap, PetscScalar *d_difvap);

/* Computes the saturation vapor density and its derivative */
void RhoVS_I(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs, PetscScalar *d_rhovs);

/* computes the phase dependent mobility*/
void Mobility(AppCtx *user, PetscScalar ice, PetscScalar *mob);

/* Computes sigma0 using logarithmic interpolation based on temperature */
void Sigma0(PetscScalar temp, PetscScalar *sigm0);

#endif // MATERIAL_PROPERTIES_H
