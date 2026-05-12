#include "io_thermal.h"

/*-----------------------------------------------------------------------------
  file_exists — return PETSC_TRUE if the file is accessible for reading.
-----------------------------------------------------------------------------*/
PetscBool file_exists(const char *filename)
{
  struct stat buf;
  return (stat(filename, &buf) == 0) ? PETSC_TRUE : PETSC_FALSE;
}

/*-----------------------------------------------------------------------------
  WriteBinaryFile
  Write a PETSc vector to a binary file.

  PARALLEL NOTE: PetscViewerBinaryOpen + VecView are collective operations.
  This function must be called by ALL MPI ranks simultaneously.
-----------------------------------------------------------------------------*/
PetscErrorCode WriteBinaryFile(Vec field, const char *filename)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscBool      ok;

  PetscFunctionBegin;

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename,
                               FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
  ierr = VecView(field, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscTestFile(filename, 'r', &ok); CHKERRQ(ierr);
  if (!ok)
    PetscPrintf(PETSC_COMM_WORLD,
                "Warning: binary output file %s was not created.\n", filename);

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  WriteIceFieldToFile
  Write Gauss-point coordinates and ice-fraction values to a text file.

  SERIAL ONLY: each MPI rank owns a subset of local elements; if called in
  parallel all ranks would write to the same file, corrupting it.  This
  function silently does nothing when running with more than one rank.
  Use it as a diagnostic tool in single-process runs only.
-----------------------------------------------------------------------------*/
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscMPIInt    size;
  FILE          *file;
  IGAElement     element;
  IGAPoint       point;
  PetscInt       indGP;

  PetscFunctionBegin;

  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if (size > 1) {
    PetscPrintf(PETSC_COMM_WORLD,
                "Note: WriteIceFieldToFile skipped in parallel run "
                "(serial diagnostic only).\n");
    PetscFunctionReturn(0);
  }

  file = fopen(filename, "w");
  if (!file)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN,
            "Cannot open ice field output file: %s", filename);

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      indGP = point->index + point->count * point->parent->index;
      if (user->dim == 2) {
        fprintf(file, "%g %g %g\n",
                point->mapX[0][0], point->mapX[0][1], user->ice[indGP]);
      } else {
        fprintf(file, "%g %g %g %g\n",
                point->mapX[0][0], point->mapX[0][1], point->mapX[0][2],
                user->ice[indGP]);
      }
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  fclose(file);
  PetscPrintf(PETSC_COMM_WORLD, "Ice field written to: %s\n", filename);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  WriteOutput
  Write user->T_sol to a binary file if outputBinary is set.

  PARALLEL NOTE: calls WriteBinaryFile which is collective — call from all
  ranks, NOT inside an if (rank == 0) guard.
-----------------------------------------------------------------------------*/
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (user->outputBinary) {
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);   CHKERRQ(ierr);
    ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  WriteKeffToCSV
  Append the dim×dim keff tensor to a CSV file (rank-0 only).
  Writes the header row the first time the file is created.

  Column order: sol_index, k_00, k_01, ..., k_{d-1}{d-1}
-----------------------------------------------------------------------------*/
PetscErrorCode WriteKeffToCSV(AppCtx *user, const char *filename,
                               PetscInt dim, const PetscReal *keff)
{
  FILE     *fp;
  PetscBool exists = file_exists(filename);

  PetscFunctionBegin;

  fp = fopen(filename, "a");
  if (!fp)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
            "Cannot open CSV file: %s", filename);

  if (!exists) {
    fprintf(fp, "sol_index");
    for (PetscInt i = 0; i < dim; i++)
      for (PetscInt j = 0; j < dim; j++)
        fprintf(fp, ",k_%d%d", (int)i, (int)j);
    fprintf(fp, "\n");
  }

  fprintf(fp, "%d", (int)user->sol_index);
  for (PetscInt k = 0; k < dim * dim; k++)
    fprintf(fp, ",%.12e", keff[k]);
  fprintf(fp, "\n");

  fclose(fp);
  PetscFunctionReturn(0);
}
