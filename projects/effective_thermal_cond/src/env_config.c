#include "env_config.h"

/*-----------------------------------------------------------------------------
  ParseOptions
  Read all simulation parameters from the PETSc options database into user.
  Hardcoded defaults are applied first; any option can be overridden by a
  -options_file or by explicit command-line flags.

  Typical usage:
    ./effective_k_ice_homog -options_file inputs/default.opts \
        -Nx 128 -Ny 128 -init_mode layered
-----------------------------------------------------------------------------*/
PetscErrorCode ParseOptions(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  /* ---- Hardcoded defaults ---- */
  user->dim          = 2;
  user->Nx           = 64;
  user->Ny           = 64;
  user->Nz           = 1;
  user->Lx           = 1.0e-3;
  user->Ly           = 1.0e-3;
  user->Lz           = 1.0e-3;
  user->eps          = 1.0e-4;
  user->T_top        = 243.15;   /* 273.15 - 30 K */
  user->q_bottom     = -1.0;     /* W/m² */
  user->outputBinary = PETSC_TRUE;
  user->sol_index    = -1;       /* -1 = loop over all files */
  user->output_dir[0] = '\0';
  strcpy(user->init_mode, "layered");
  user->init_dir[0]  = '\0';

  /* ---- Read from options database ---- */
  PetscOptionsBegin(PETSC_COMM_WORLD, "",
                    "Effective Thermal Conductivity (Homogenization)", "");

  ierr = PetscOptionsInt("-dim",
      "Problem dimension (2 or 3)",
      __FILE__, user->dim, &user->dim, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-Nx",
      "Number of elements in x",
      __FILE__, user->Nx, &user->Nx, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-Ny",
      "Number of elements in y",
      __FILE__, user->Ny, &user->Ny, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-Nz",
      "Number of elements in z (3-D only)",
      __FILE__, user->Nz, &user->Nz, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-Lx",
      "Domain length in x (m)",
      __FILE__, user->Lx, &user->Lx, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-Ly",
      "Domain length in y (m)",
      __FILE__, user->Ly, &user->Ly, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-Lz",
      "Domain length in z (m, 3-D only)",
      __FILE__, user->Lz, &user->Lz, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-eps",
      "Phase-field interface width (m)",
      __FILE__, user->eps, &user->eps, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-temp_top",
      "Temperature at the top boundary (K)",
      __FILE__, user->T_top, &user->T_top, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-flux_bottom",
      "Heat flux at the bottom boundary (W/m²)",
      __FILE__, user->q_bottom, &user->q_bottom, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-output_binary",
      "Write binary t_vec.dat solution output",
      __FILE__, user->outputBinary, &user->outputBinary, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-sol_index",
      "Solution file index to process (-1 = all files)",
      __FILE__, user->sol_index, &user->sol_index, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-output_dir",
      "Directory for output files",
      __FILE__, user->output_dir, user->output_dir,
      sizeof(user->output_dir), NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-init_mode",
      "Ice field initialisation mode: circle | layered | file",
      __FILE__, user->init_mode, user->init_mode,
      sizeof(user->init_mode), NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-init_dir",
      "Directory containing sol_*.dat input files (file mode)",
      __FILE__, user->init_dir, user->init_dir,
      sizeof(user->init_dir), NULL); CHKERRQ(ierr);

  PetscOptionsEnd();

  /* ---- Basic validation ---- */
  if (user->dim != 2 && user->dim != 3)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "-dim must be 2 or 3 (got %d)", (int)user->dim);

  if (user->Nx < 1 || user->Ny < 1 || (user->dim == 3 && user->Nz < 1))
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "Mesh resolution must be >= 1 in every direction");

  if (user->eps <= 0.0)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
            "-eps must be positive (got %g)", (double)user->eps);

  if (user->output_dir[0] == '\0')
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
            "-output_dir is required (no output directory specified)");

  PetscFunctionReturn(0);
}
