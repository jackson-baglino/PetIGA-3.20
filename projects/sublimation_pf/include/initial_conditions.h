#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include "NASA_types.h"

/* 1D initial conditions */
PetscErrorCode FormInitialCondition1D(IGA iga, Vec U, AppCtx *user);

/* 2D Ice Slab — 1D-equivalent centered slab, uniform in y */
PetscErrorCode FormInitialIceSlab2D(IGA iga, Vec U, AppCtx *user);

/* Single ice grain centred in the domain (no sediment) */
PetscErrorCode FormInitialSingleIceGrain2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialSingleIceGrain1D(IGA iga, Vec U, AppCtx *user);

/* Two ice semicircles centered on the x=0 and x=Lx boundaries (Ostwald ripening test) */
PetscErrorCode FormInitialTwoIceGrainsBoundary2D(IGA iga, Vec U, AppCtx *user);

/* N ice grains (centers/radii from -ice_grain_cx/-ice_grain_cy/-ice_grain_R)
 * on an optionally multi-bump sediment geometry (-sed_grain_x/-sed_grain_R) */
PetscErrorCode FormInitialMultiGrains2D(IGA iga, Vec U, AppCtx *user);

/* 3-phase (dof=4): flat sediment slab (y < sed_slab_height) + one ice grain.
 * sed_slab_height <= 0 => phi_s=0 everywhere (2-phase validation mode). */
PetscErrorCode FormInitialSedSlabGrain2D(IGA iga, Vec U, AppCtx *user);

/* 3-phase (dof=4) 1D: sediment|ice|air stack; sed_slab_height<=0 => phi_s=0
 * centred grain (2-phase validation). */
PetscErrorCode FormInitialSedIce1D(IGA iga, Vec U, AppCtx *user);

#endif // INITIAL_CONDITIONS_H
