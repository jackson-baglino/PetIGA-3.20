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

/* ---- Three-phase (ice/sediment/air) variants — 3-phase model only (dof=4) ---- */
/* Effective properties over (ice, sed); air = 1 - ice - sed. Derivatives are
 * w.r.t. the two independent DOFs. */
void ThermalCond3(AppCtx *user, PetscScalar ice, PetscScalar sed,
                  PetscScalar *cond, PetscScalar *dcond_ice, PetscScalar *dcond_sed);
void HeatCap3(AppCtx *user, PetscScalar ice, PetscScalar sed,
              PetscScalar *cp, PetscScalar *dcp_ice, PetscScalar *dcp_sed);
void Density3(AppCtx *user, PetscScalar ice, PetscScalar sed,
              PetscScalar *rho, PetscScalar *drho_ice, PetscScalar *drho_sed);

/* Triple-well driving force g = dF^tri/dphi_i - dF^tri/dphi_a and its
 * derivatives w.r.t. (phi_i, phi_s). See material_properties.c for the form. */
void TripleWell(AppCtx *user, PetscScalar ice, PetscScalar sed,
                PetscScalar *g, PetscScalar *dg_ice, PetscScalar *dg_sed);

#endif // MATERIAL_PROPERTIES_H
