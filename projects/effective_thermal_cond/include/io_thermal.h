#ifndef IO_THERMAL_H
#define IO_THERMAL_H

#include "user_context.h"
#include "material_properties.h"

/* -------------------------- Function Declarations ------------------------- */
PetscErrorCode ReadSolutionVec(const char *iga_file, const char *vec_file, IGA *iga_out, 
                              AppCtx *user);
PetscErrorCode LoadIceField(AppCtx *user, const char *iga_filename, const char *vec_filename);
PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K);
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename);
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user);
PetscErrorCode WriteBinaryFile(Vec field, const char *filename);

#endif // ENV_CONFIG_H