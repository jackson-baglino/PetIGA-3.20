#ifndef MONITORING_H
#define MONITORING_H

#include "enceladus_types.h"

/* Monitor function for TS solver */
PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Output monitor function for writing solution data */
PetscErrorCode OutputMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Interface-CFL timestep limiter: clamps the next dt from the last accepted
   step's max pointwise phase-change rate (see monitoring.c) */
PetscErrorCode InterfaceCFLMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* TS pre-step callback: consumes deferred bounds-rollback requests set by Monitor() */
PetscErrorCode BoundsRollbackPreStep(TS ts);

/* Post-setup memory report, and a fail-fast check against the allocation.
   Called once, after every large allocation exists and before the expensive
   loop starts, so an under-sized job dies in seconds instead of occupying
   cores for hours. No-op guard unless a batch scheduler told us the budget. */
PetscErrorCode MemoryBudgetCheck(AppCtx *user, const char *stage);

#endif // MONITORING_H