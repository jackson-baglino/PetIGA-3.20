#ifndef IO_H
#define IO_H

#include "app_ctx.h"

/*
 * file_exists — return PETSC_TRUE if the named file is accessible.
 */
PetscBool file_exists(const char *filename);

/*
 * WriteBinaryFile — write a PETSc vector to a binary file (collective).
 * Must be called by all MPI ranks simultaneously.
 */
PetscErrorCode WriteBinaryFile(Vec field, const char *filename);

/*
 * WriteIceFieldToFile — write Gauss-point coordinates and ice values to a
 * human-readable text file.  Serial-only: guarded internally with a
 * MPI_Comm_size check; silently skipped in parallel runs.
 */
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user);

/*
 * WriteOutput — write user->T_sol to a binary file if outputBinary is set.
 * Must be called collectively by all ranks (PETSc binary I/O is collective).
 */
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename);

/*
 * WriteKeffToCSV — append the keff tensor (row-major, dim×dim) to a CSV file.
 * Writes the header on the first call if the file does not exist.
 * Call from rank 0 only.
 */
PetscErrorCode WriteKeffToCSV(AppCtx *user, const char *filename,
                               PetscInt dim, const PetscReal *keff);

#endif /* IO_H */
