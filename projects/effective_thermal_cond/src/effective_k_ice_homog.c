/*==============================================================================
  effective_k_ice_homog.c — main entry point

  Computes the effective thermal conductivity tensor of a heterogeneous ice/air
  medium via periodic homogenization.  For each input microstructure (phase
  field snapshot), the code:

    1. Sets up a periodic IGA discretisation.
    2. Reads or generates the ice phase-field on Gauss points.
    3. Solves the cell problem:  -div( k(x) grad t ) = div( k(x) I )
       with periodic BCs, yielding the correction temperature t.
    4. Integrates  k_eff[i,j] = < k (grad_t[i,j] + delta_ij) >  over the cell.
    5. Writes k_eff to k_eff.csv.

  Usage:
    ./effective_k_ice_homog -options_file inputs/default.opts [overrides]
    mpiexec -np N ./effective_k_ice_homog -options_file run.opts

  See inputs/default.opts for a full list of options.  Run with -help for
  on-line documentation.
==============================================================================*/

#include "app_ctx.h"
#include "env_config.h"
#include "setup_thermal.h"
#include "field_init.h"
#include "assembly.h"
#include "solver.h"
#include "io_thermal.h"

int main(int argc, char *argv[])
{
  AppCtx         user;
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  PetscLogDouble t_start, t_1, t_2, t_end;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  ierr = PetscTime(&t_start); CHKERRQ(ierr);

  /* ---- Parse all simulation parameters from options database ---- */
  ierr = ParseOptions(&user); CHKERRQ(ierr);
  InitializeUserContext(&user);

  PetscPrintf(PETSC_COMM_WORLD,
              "Domain : Lx=%.3g  Ly=%.3g  Lz=%.3g  (dim=%d)\n",
              user.Lx, user.Ly, user.Lz, user.dim);
  PetscPrintf(PETSC_COMM_WORLD,
              "Mesh   : Nx=%d  Ny=%d  Nz=%d  eps=%.3g\n",
              user.Nx, user.Ny, user.Nz, user.eps);
  PetscPrintf(PETSC_COMM_WORLD,
              "BCs    : T_top=%.2f K  flux_bottom=%.4g W/m^2\n",
              user.T_top, user.q_bottom);
  PetscPrintf(PETSC_COMM_WORLD,
              "Output : %s  (sol_index=%d, binary=%s)\n\n",
              user.output_dir, user.sol_index,
              user.outputBinary ? "yes" : "no");

  /* ---- Ensure output directory exists ---- */
  {
    struct stat st;
    if (rank == 0 && stat(user.output_dir, &st) != 0)
      mkdir(user.output_dir, 0755);
    MPI_Barrier(PETSC_COMM_WORLD);
  }

  /* ---- Set up IGA (shared across all file iterations) ---- */
  IGA iga;
  ierr = SetupIGA(&user, &iga); CHKERRQ(ierr);
  user.iga = iga;

  /* ---- Create solver objects once (reused across file iterations) ---- */
  Mat A;
  Vec b;
  KSP ksp;
  ierr = CreateSolverObjects(iga, &user, &A, &b, &ksp); CHKERRQ(ierr);

  /* ---- keff output buffer ---- */
  PetscReal keff[9]; /* dim×dim, max 3×3 */

  /* ================================================================
     Branch A: process a single solution file
     ================================================================ */
  if (user.sol_index >= 0) {
    ierr = PetscTime(&t_1); CHKERRQ(ierr);

    ierr = FormInitialCondition(&user); CHKERRQ(ierr);
    ierr = Solve(&user, iga, A, b, ksp); CHKERRQ(ierr);

    ierr = PetscTime(&t_2); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,
                "KSP solve completed in %.2f s.\n\n", t_2 - t_1);

    ierr = ComputeKeffective(iga, user.T_sol, keff, &user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Effective thermal conductivity:\n");
    for (PetscInt i = 0; i < user.dim; i++)
      for (PetscInt j = 0; j < user.dim; j++)
        PetscPrintf(PETSC_COMM_WORLD, "  k_eff[%d][%d] = %.12e\n",
                    i, j, keff[i * user.dim + j]);
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    /* Binary output is collective — call from all ranks */
    char t_vec_file[PETSC_MAX_PATH_LEN];
    snprintf(t_vec_file, sizeof(t_vec_file), "%s/t_vec.dat", user.output_dir);
    ierr = WriteOutput(&user, user.T_sol, t_vec_file); CHKERRQ(ierr);

    /* Text / CSV outputs — rank 0 only */
    if (rank == 0) {
      char ice_file[PETSC_MAX_PATH_LEN];
      snprintf(ice_file, sizeof(ice_file), "%s/ice_data.dat", user.output_dir);
      ierr = WriteIceFieldToFile(ice_file, &user); CHKERRQ(ierr);

      char csv_file[PETSC_MAX_PATH_LEN];
      snprintf(csv_file, sizeof(csv_file), "%s/k_eff.csv", user.output_dir);
      ierr = WriteKeffToCSV(&user, csv_file, user.dim, keff); CHKERRQ(ierr);
    }

    ierr = PetscFree(user.ice); CHKERRQ(ierr);

  /* ================================================================
     Branch B: loop over all sol_*.dat files in init_dir
     ================================================================ */
  } else {
    PetscPrintf(PETSC_COMM_WORLD,
                "sol_index < 0 — scanning all files in: %s\n\n",
                user.init_dir);

    /* Directory scan on rank 0; broadcast result to all ranks */
    PetscInt  file_indices[10000];
    PetscInt  num_files = 0;

    if (rank == 0) {
      DIR           *dir;
      struct dirent *entry;
      regex_t        regex;

      dir = opendir(user.init_dir);
      if (!dir) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "Error: cannot open init_dir: %s\n", user.init_dir);
        MPI_Abort(PETSC_COMM_WORLD, 1);
      }

      regcomp(&regex, "^sol_([0-9]+)\\.dat$", REG_EXTENDED);
      while ((entry = readdir(dir)) != NULL) {
        if (regexec(&regex, entry->d_name, 0, NULL, 0) == 0) {
          int index;
          sscanf(entry->d_name, "sol_%05d.dat", &index);
          if (num_files < 10000)
            file_indices[num_files++] = index;
        }
      }
      closedir(dir);
      regfree(&regex);
      PetscSortInt(num_files, file_indices);
    }

    /* Broadcast the file list to all ranks */
    MPI_Bcast(&num_files, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(file_indices, num_files, MPI_INT, 0, PETSC_COMM_WORLD);

    PetscPrintf(PETSC_COMM_WORLD,
                "Found %d solution files.\n\n", num_files);

    for (PetscInt fi = 0; fi < num_files; fi++) {
      ierr = PetscTime(&t_1); CHKERRQ(ierr);

      user.sol_index = file_indices[fi];
      PetscPrintf(PETSC_COMM_WORLD,
                  "--- Processing sol_%05d ---\n\n", user.sol_index);

      /* Free previous ice field before re-allocating */
      if (user.ice) { ierr = PetscFree(user.ice); CHKERRQ(ierr); }

      ierr = FormInitialCondition(&user); CHKERRQ(ierr);
      ierr = Solve(&user, iga, A, b, ksp); CHKERRQ(ierr);

      ierr = ComputeKeffective(iga, user.T_sol, keff, &user); CHKERRQ(ierr);
      for (PetscInt i = 0; i < user.dim; i++)
        for (PetscInt j = 0; j < user.dim; j++)
          PetscPrintf(PETSC_COMM_WORLD, "  k_eff[%d][%d] = %.12e\n",
                      i, j, keff[i * user.dim + j]);
      PetscPrintf(PETSC_COMM_WORLD, "\n");

      /* Binary output — collective */
      char t_vec_file[PETSC_MAX_PATH_LEN];
      snprintf(t_vec_file, sizeof(t_vec_file), "%s/t_vec_%05d.dat",
               user.output_dir, user.sol_index);
      ierr = WriteOutput(&user, user.T_sol, t_vec_file); CHKERRQ(ierr);

      /* CSV output — rank 0 */
      if (rank == 0) {
        char csv_file[PETSC_MAX_PATH_LEN];
        snprintf(csv_file, sizeof(csv_file), "%s/k_eff.csv", user.output_dir);
        ierr = WriteKeffToCSV(&user, csv_file, user.dim, keff); CHKERRQ(ierr);
      }

      ierr = PetscTime(&t_2); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD,
                  "sol_%05d done in %.2f s\n"
                  "-------------------------------------------------\n\n",
                  user.sol_index, t_2 - t_1);
    }

    if (user.ice) { ierr = PetscFree(user.ice); CHKERRQ(ierr); }
  }

  /* ---- Cleanup ---- */
  ierr = DestroySolverObjects(&A, &b, &ksp); CHKERRQ(ierr);
  ierr = VecDestroy(&user.T_sol); CHKERRQ(ierr);
  ierr = IGADestroy(&iga); CHKERRQ(ierr);

  ierr = PetscTime(&t_end); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,
              "Total elapsed time: %.2f s\nProgram completed successfully.\n",
              t_end - t_start);

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
