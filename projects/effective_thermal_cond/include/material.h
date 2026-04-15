#ifndef MATERIAL_H
#define MATERIAL_H

#include "app_ctx.h"

/*
 * ThermalCond — linear mixture rule for ice/air thermal conductivity.
 *
 * Computes:
 *   cond      = ice * thcond_ice + (1-ice) * thcond_air
 *   dcond_ice = thcond_ice - thcond_air   (derivative w.r.t. ice fraction)
 *
 * Pass NULL for any output you don't need.
 */
void ThermalCond(AppCtx *user, PetscScalar ice,
                 PetscScalar *cond, PetscScalar *dcond_ice);

#endif /* MATERIAL_H */
