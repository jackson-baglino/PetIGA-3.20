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

#endif // INITIAL_CONDITIONS_H
