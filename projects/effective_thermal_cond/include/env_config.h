#ifndef ENV_CONFIG_H
#define ENV_CONFIG_H

#include "app_ctx.h"

/*
 * ParseOptions — read all simulation parameters from the PETSc options
 * database (populated from -options_file or the command line) into user.
 *
 * Hardcoded defaults are applied first; every option can be overridden by
 * passing it explicitly on the command line or in an options file via
 *   -options_file path/to/run.opts
 *
 * Run with -help to see a full list of supported options.
 */
PetscErrorCode ParseOptions(AppCtx *user);

#endif /* ENV_CONFIG_H */
