#include "NASA_main.h" // Need to change name later

int main(int argc, char *argv[]) {
    /* Petsc Initialization */
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

    /* Start timer */
    PetscLogDouble itim;
    ierr = PetscTime(&itim); CHKERRQ(ierr);

    /* Get number of processes (number of cores used) */
    PetscInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    PetscPrintf(PETSC_COMM_WORLD, "Running on %d processes.\n\n\n", size);

    /* Define simulation specific parameters */
    AppCtx user;                         /* User-defined application context */
    PetscMemzero(&user, sizeof(AppCtx)); /* Initialize user context to zero */
    PetscBool flag_BC_Tfix;              /* fix temperature at boundaries */
    PetscBool flag_BC_rhovfix;           /* fix vapor density at boundaries */

    user.Lambd      = 1.0;      /* Model parameter Lambda (triple-junction penalty) */
    user.air_lim    = 1.0e-6;   /* Air phase fraction floor */

    user.lat_sub    = 2.83e6;   /* Latent heat of sublimation */

    user.thcond_ice = 2.29;     /* Thermal conductivity of ice */
    user.thcond_air = 0.02;     /* Thermal conductivity of air */

    user.cp_ice     = 1.96e3;   /* Specific heat capacity of ice */
    user.cp_air     = 1.044e3;  /* Specific heat capacity of air */

    user.rho_ice    = 919.0;    /* Density of ice */
    user.rho_air    = 1.341;    /* Density of air */

    user.dif_vap    = 2.178e-5; /* Vapor diffusivity in air */

    user.T_melt     = 0.0;      /* Melting temperature of ice */

    user.flag_tIC   = 0;              /* IC variant: 0=centered slab, 2=flat interface */
    user.readFlag   = PETSC_FALSE;    /* read initial field from file */
    user.flag_Tdep  = PETSC_FALSE;    /* temperature-dependent material properties */

    user.phase_lo   = -0.05;   /* lower bound: phi below this → abort */
    user.phase_hi   =  1.05;   /* upper bound: phi above this → abort */

    user.d0_sub0    = 1.0e-9;    /* capillary length scale (physical) */
    user.beta_sub0  = 1.4e5;     /* kinetic coefficient (physical) */

    /* Surface energy parameters of the double-well free energy [J/m²]:
     *   F_dub(phi_i) = C*phi_i^2(1-phi_i)^2,  C = (Sigma_i+Sigma_a)/2 + Lambda
     * Sigma_i = ice-side surface energy, Sigma_a = air-side surface energy. */
    PetscReal Sigma_i = 0.109; /* ice surface energy [J/m²] */
    PetscReal Sigma_a = 0.132; /* air surface energy [J/m²] */

    /* Define common variables (can be overridden by PETSc options) */
    PetscInt  p   = 1;          /* Polynomial order */
    PetscInt  C   = 0;          /* Global continuity order */

    PetscInt  dof = 3;          /* Degrees of freedom per node (ice, temperature, vapor) */
    PetscInt  dim = 2;          /* Problem dimension (2D or 3D) */

    PetscInt  Nx  = 64;         /* Number of elements in x direction */
    PetscInt  Ny  = 64;         /* Number of elements in y direction */
    PetscInt  Nz  = 64;         /* Number of elements in z direction */

    PetscReal Lx  = 1.0e-3;     /* Domain length in x direction */
    PetscReal Ly  = 1.0e-3;     /* Domain length in y direction */
    PetscReal Lz  = 1.0e-3;     /* Domain length in z direction */

    PetscReal delt_t = 1.0e-4;  /* Time step size */
    PetscReal t_final = 0.0;    /* Final simulation time (does not advance if 0) */

    PetscInt  n_out   = 10;     /* Number of outputs */

    PetscReal humidity = 0.95;  /* Initial humidity */
    PetscReal temp     = -20.0; /* Initial temperature */

    PetscReal grad_temp0[3] = {0.0, 0.0, 0.0}; /* Initial temperature gradient */

    PetscReal eps = 9.0e-7;     /* Interface width parameter */

    /* Define grain parameters (can be overridden by PETSc options) */
    user.NCice       = 50;      /* Number of ice grains */
    user.RCice       = 0.3e-4;  /* Mean radius */
    user.RCice_dev   = 0.55;    /* Std dev of radius */

    /* Define boundary condition flags (can be overridden by PETSc options).
     * Default is INSULATING for both: no Dirichlet condition is registered,
     * no surface form is assembled in Residual (see `pnt->atboundary` early
     * return in assembly.c), so the natural BC ∂u/∂n = 0 is enforced — i.e.
     * zero flux through the boundary for both T and rho_v. Override only if
     * you actually want fixed-value Dirichlet conditions. */
    user.periodic    = 0;       /* Periodic boundary condition flag */
    flag_BC_Tfix     = PETSC_FALSE; /* insulating T (zero heat flux) — natural Neumann */
    flag_BC_rhovfix  = PETSC_FALSE; /* insulating rho_v (zero vapor flux) — natural Neumann */

    /* Define output parameters (can be overridden by PETSc options) */
    user.outp        = 0;       /* Output control flag (0: output according to t_interv) */
    user.t_out       = 0.0;     /* Next output time */
    if (n_out > 1) {
        user.t_interv = t_final / (n_out - 1); /* Output interval */
    } else {
        user.t_interv = t_final;
    }

    /* Adaptive time stepping parameters (can be overridden by PETSc options) */
    PetscInt  adap    = 1;                             /* Adaptive time stepping flag */
    PetscInt  NRmin   = 3;                             /* Minimum Newton-Raphson iterations */
    PetscInt  NRmax   = 5;                             /* Maximum Newton-Raphson iterations */
    PetscReal factor  = pow(10.0, 1.0 / 8.0);          /* Time step adjustment factor */
    PetscReal dtmin   = 0.0;                           /* Minimum time step size */
    PetscReal dtmax   = 0.0;                           /* Maximum time step size */
    PetscInt  max_rej = 10;                            /* Maximum number of rejected steps */

    /* Get simulation parameters from CLI .txt file (PETSc options) */
    PetscBool output  = PETSC_TRUE;                    /* Output flag */
    PetscBool monitor = PETSC_TRUE;                    /* Monitor flag */
    char      initial[PETSC_MAX_PATH_LEN] = {0};       /* Initial condition file */
    char      PFgeom[PETSC_MAX_PATH_LEN]  = {0};       /* Initial ice geometry file */
    char      geom_file[PETSC_MAX_PATH_LEN] = {0};     /* igakit-generated IGA geometry (.dat), overrides axis setup */
    char      ic_type[64]                 = "two_ice_grains_boundary"; /* IC geometry selector */
    PetscBool flg = PETSC_FALSE;

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Permafrost options", "IGA");
    /* --- Geometry & discretization --------------------------------------- */

    ierr = PetscOptionsInt("-dof", "Degrees of freedom per node", "", dof, &dof, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-dim", "Problem dimension (2 or 3)", "", dim, &dim, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-Nx", "Number of elements in x direction", "", Nx, &Nx, NULL); CHKERRQ(ierr);
    if (dim >= 2) {
        ierr = PetscOptionsInt("-Ny", "Number of elements in y direction", "", Ny, &Ny, NULL); CHKERRQ(ierr);
    }
    if (dim == 3) {
        ierr = PetscOptionsInt("-Nz", "Number of elements in z direction", "", Nz, &Nz, NULL); CHKERRQ(ierr);
    }
    PetscInt ngrad = dim; /* Number of grad_temp0 components to read */
    ierr = PetscOptionsReal("-Lx", "Domain length in x direction", "", Lx, &Lx, NULL); CHKERRQ(ierr);
    if (dim >= 2) {
        ierr = PetscOptionsReal("-Ly", "Domain length in y direction", "", Ly, &Ly, NULL); CHKERRQ(ierr);
    }
    if (dim == 3) {
        ierr = PetscOptionsReal("-Lz", "Domain length in z direction", "", Lz, &Lz, NULL); CHKERRQ(ierr);
    }
    ierr = PetscOptionsInt("-p", "Polynomial order", "", p, &p, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-C", "Global continuity order", "", C, &C, NULL); CHKERRQ(ierr);

    /* --- Time stepping & output cadence ---------------------------------- */
    ierr = PetscOptionsReal("-delt_t", "Time step size", "", delt_t, &delt_t, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-t_final", "Final simulation time", "", t_final, &t_final, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-n_out", "Number of outputs", "", n_out, &n_out, NULL); CHKERRQ(ierr);

    /* --- Initial conditions: environment --------------------------------- */
    ierr = PetscOptionsReal("-humidity", "Initial humidity", "", humidity, &humidity, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-temp", "Initial temperature", "", temp, &temp, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsRealArray("-grad_temp0",
                                 "Initial temperature gradient [dT/dx dT/dy dT/dz]",
                                 "",
                                 grad_temp0, &ngrad, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-eps", "Interface width parameter", "", eps, &eps, NULL); CHKERRQ(ierr);

    /* --- Grain geometry: ice --------------------------------------------- */
    ierr = PetscOptionsInt("-NCice", "Number of ice grains", "", user.NCice, &user.NCice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice", "Mean radius of ice grains", "", user.RCice, &user.RCice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice_dev", "Std dev of radius of ice grains", "", user.RCice_dev, &user.RCice_dev, NULL); CHKERRQ(ierr);

    /* --- Per-grain radii (two-grain boundary IC) -------------------------- */
    user.RCice0 = user.RCice;
    user.RCice1 = user.RCice;
    ierr = PetscOptionsReal("-RCice0", "Radius of grain 0 (x=0 boundary, two_ice_grains_boundary IC)", "", user.RCice0, &user.RCice0, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice1", "Radius of grain 1 (x=Lx boundary, two_ice_grains_boundary IC)", "", user.RCice1, &user.RCice1, NULL); CHKERRQ(ierr);

    /* --- Sediment-grain bump geometry (must match -geom_file) ------------- */
    user.geom_bump_R = 0.0;
    ierr = PetscOptionsReal("-geom_bump_R", "Sediment-grain bump half-width/height (must match build_geometry_sediment_grain.py's R_sed; 0 = flat domain)", "", user.geom_bump_R, &user.geom_bump_R, NULL); CHKERRQ(ierr);

    /* --- Multi-grain geometry: sediment bumps + ice grains ---------------- */
    user.n_sed_grains = 0;
    {
        PetscInt n = MAX_SED_GRAINS;
        ierr = PetscOptionsRealArray("-sed_grain_x",
                 "Sediment bump center x-positions [m]; the bottom edge of "
                 "-geom_file is the sum of SedimentBump() humps at these "
                 "centers (must match build_geometry_multi_grain.py's "
                 "SEDIMENT_GRAINS). Overrides -geom_bump_R single-bump IC.",
                 "", user.sed_grain_x, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            user.n_sed_grains = n;
            PetscInt nr = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-sed_grain_R",
                     "Sediment bump half-widths [m], one per -sed_grain_x entry",
                     "", user.sed_grain_R, &nr, NULL); CHKERRQ(ierr);
            if (nr != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-sed_grain_x and -sed_grain_R must have the same length (%d vs %d)",
                        (int)n, (int)nr);
            /* -sed_grain_h: peak heights [m]. Defaults to R if not provided. */
            PetscInt nh = MAX_SED_GRAINS;
            PetscBool hflg;
            ierr = PetscOptionsRealArray("-sed_grain_h",
                     "Sediment bump peak heights [m], one per -sed_grain_x entry "
                     "(defaults to sed_grain_R if omitted)",
                     "", user.sed_grain_h, &nh, &hflg); CHKERRQ(ierr);
            if (!hflg) {
                for (PetscInt k = 0; k < n; k++) user.sed_grain_h[k] = user.sed_grain_R[k];
            } else if (nh != n) {
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-sed_grain_x and -sed_grain_h must have the same length (%d vs %d)",
                        (int)n, (int)nh);
            }
        }
    }

    /* --- Top-wall (ceiling) bump geometry ---------------------------------- */
    user.n_top_grains = 0;
    {
        PetscInt n = MAX_SED_GRAINS;
        ierr = PetscOptionsRealArray("-top_grain_x",
                 "Ceiling bump center x-positions [m]; bumps push DOWN from Ly "
                 "(must match TOP_GRAINS in build_geometry_multi_grain.py)",
                 "", user.top_grain_x, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            user.n_top_grains = n;
            PetscInt nr = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-top_grain_R",
                     "Ceiling bump half-widths [m], one per -top_grain_x entry",
                     "", user.top_grain_R, &nr, NULL); CHKERRQ(ierr);
            if (nr != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-top_grain_x and -top_grain_R must have the same length (%d vs %d)",
                        (int)n, (int)nr);
            PetscInt nh = MAX_SED_GRAINS;
            PetscBool hflg;
            ierr = PetscOptionsRealArray("-top_grain_h",
                     "Ceiling bump peak heights [m] (defaults to top_grain_R if omitted)",
                     "", user.top_grain_h, &nh, &hflg); CHKERRQ(ierr);
            if (!hflg) {
                for (PetscInt k = 0; k < n; k++) user.top_grain_h[k] = user.top_grain_R[k];
            } else if (nh != n) {
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-top_grain_x and -top_grain_h must have the same length (%d vs %d)",
                        (int)n, (int)nh);
            }
        }
    }

    /* --- Ice "shell" capping a floor bump at constant thickness ------------ */
    user.n_ice_shells = 0;
    {
        PetscInt n = MAX_SED_GRAINS;
        ierr = PetscOptionsRealArray("-ice_shell_x",
                 "Ice-shell lateral center x-positions [m]; a smooth band of "
                 "constant vertical thickness sitting on SedimentBumpField(x), "
                 "windowed to [x-R, x+R] so it only covers the bump under it "
                 "(added on top of -ice_grain_* ellipses, not a replacement)",
                 "", user.ice_shell_x, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            user.n_ice_shells = n;
            PetscInt nr = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-ice_shell_R",
                     "Ice-shell lateral half-widths [m], one per -ice_shell_x entry",
                     "", user.ice_shell_R, &nr, NULL); CHKERRQ(ierr);
            if (nr != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-ice_shell_x and -ice_shell_R must have the same length (%d vs %d)",
                        (int)n, (int)nr);
            PetscInt nt = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-ice_shell_thickness",
                     "Ice-shell constant vertical thickness [m], one per -ice_shell_x entry",
                     "", user.ice_shell_thickness, &nt, NULL); CHKERRQ(ierr);
            if (nt != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-ice_shell_x and -ice_shell_thickness must have the same length (%d vs %d)",
                        (int)n, (int)nt);
        }
    }

    /* --- Flat ice layer encapsulating a floor bump -------------------------- */
    user.n_ice_flats = 0;
    {
        PetscInt n = MAX_SED_GRAINS;
        ierr = PetscOptionsRealArray("-ice_flat_x",
                 "Flat-ice-layer lateral center x-positions [m]; ice fills "
                 "everything below the ABSOLUTE height -ice_flat_height (not "
                 "relative to the bump's own surface like -ice_shell_thickness), "
                 "windowed to [x-R, x+R] -- gives a flat, non-rounded ice-air "
                 "interface burying the bump instead of a domed/conformal cap",
                 "", user.ice_flat_x, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            user.n_ice_flats = n;
            PetscInt nr = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-ice_flat_R",
                     "Flat-ice-layer lateral half-widths [m], one per -ice_flat_x entry",
                     "", user.ice_flat_R, &nr, NULL); CHKERRQ(ierr);
            if (nr != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-ice_flat_x and -ice_flat_R must have the same length (%d vs %d)",
                        (int)n, (int)nr);
            PetscInt nh = MAX_SED_GRAINS;
            ierr = PetscOptionsRealArray("-ice_flat_height",
                     "Flat-ice-layer absolute top height [m], one per -ice_flat_x entry",
                     "", user.ice_flat_height, &nh, NULL); CHKERRQ(ierr);
            if (nh != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-ice_flat_x and -ice_flat_height must have the same length (%d vs %d)",
                        (int)n, (int)nh);
        }
    }

    /* --- Multi-grain ice IC (-ic_type multi_grains) ------------------------ */
    user.n_act = 0;
    {
        PetscInt  n = 200;
        PetscReal cx[200];
        ierr = PetscOptionsRealArray("-ice_grain_cx",
                 "Ice grain center x-positions [m] (-ic_type multi_grains)",
                 "", cx, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            PetscInt  ncy = 200, nr = 200;
            PetscReal cy[200], rr[200];
            ierr = PetscOptionsRealArray("-ice_grain_cy", "Ice grain center y-positions [m]", "", cy, &ncy, NULL); CHKERRQ(ierr);
            ierr = PetscOptionsRealArray("-ice_grain_R",  "Ice grain radii [m]",               "", rr, &nr,  NULL); CHKERRQ(ierr);
            if (ncy != n || nr != n)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-ice_grain_cx, -ice_grain_cy, -ice_grain_R must have the same length (%d, %d, %d)",
                        (int)n, (int)ncy, (int)nr);
            user.n_act = n;
            for (PetscInt k = 0; k < n; k++) {
                user.cent[0][k] = cx[k];
                user.cent[1][k] = cy[k];
                user.radius[k]  = rr[k];
            }
        }
    }

    /* --- Elliptical ice grain semi-axes (-ice_grain_ax / -ice_grain_ay) -- */
    {
        PetscInt  nax = 200, nay = 200;
        PetscBool axflg, ayflg;
        PetscReal tmp[200];
        ierr = PetscOptionsRealArray("-ice_grain_ax",
                 "Ice grain semi-axis in x [m] (elliptical grains; defaults to -ice_grain_R)",
                 "", tmp, &nax, &axflg); CHKERRQ(ierr);
        if (axflg) {
            for (PetscInt k = 0; k < user.n_act; k++)
                user.ice_grain_ax[k] = (k < nax) ? tmp[k] : user.radius[k];
        } else {
            for (PetscInt k = 0; k < user.n_act; k++)
                user.ice_grain_ax[k] = user.radius[k];
        }
        ierr = PetscOptionsRealArray("-ice_grain_ay",
                 "Ice grain semi-axis in y [m] (elliptical grains; defaults to -ice_grain_R)",
                 "", tmp, &nay, &ayflg); CHKERRQ(ierr);
        if (ayflg) {
            for (PetscInt k = 0; k < user.n_act; k++)
                user.ice_grain_ay[k] = (k < nay) ? tmp[k] : user.radius[k];
        } else {
            for (PetscInt k = 0; k < user.n_act; k++)
                user.ice_grain_ay[k] = user.radius[k];
        }
    }

    /* --- Boundary conditions & physics flags ----------------------------- */
    ierr = PetscOptionsInt("-periodic", "Periodic boundary condition flag", "", user.periodic, &user.periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_Tfix",    "Fix temperature at boundaries",                    "", flag_BC_Tfix,    &flag_BC_Tfix,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_rhovfix", "Fix vapor density at boundaries",                  "", flag_BC_rhovfix, &flag_BC_rhovfix, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_Tdep",       "Temperature-dependent Gibbs-Thomson parameters",   "", user.flag_Tdep,  &user.flag_Tdep,  NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_tIC", "1D IC variant (0=centered slab, 2=flat interface)", "", user.flag_tIC, &user.flag_tIC, NULL); CHKERRQ(ierr);
    /* --- Thermophysical properties --------------------------------------- */
    ierr = PetscOptionsReal("-thcond_ice", "Thermal conductivity of ice", "", user.thcond_ice, &user.thcond_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thcond_air", "Thermal conductivity of air", "", user.thcond_air, &user.thcond_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-cp_ice", "Specific heat capacity of ice", "", user.cp_ice, &user.cp_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cp_air", "Specific heat capacity of air", "", user.cp_air, &user.cp_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-rho_ice", "Density of ice", "", user.rho_ice, &user.rho_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rho_air", "Density of air", "", user.rho_air, &user.rho_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-Sigma_i", "Ice-side surface energy in the double-well free energy [J/m^2]", "", Sigma_i, &Sigma_i, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Sigma_a", "Air-side surface energy in the double-well free energy [J/m^2]", "", Sigma_a, &Sigma_a, NULL); CHKERRQ(ierr);

    /* --- Output control -------------------------------------------------- */
    ierr = PetscOptionsInt("-outp", "Output control flag", "", user.outp, &user.outp, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-t_interv", "Output interval", "", user.t_interv, &user.t_interv, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-Permafrost_output", "Enable output files", "", output, &output, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-Permafrost_monitor", "Monitor the solution", "", monitor, &monitor, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-output_path", "Output file path", "", user.output_path, user.output_path, sizeof(user.output_path), NULL); CHKERRQ(ierr);

    /* --- Adaptive time stepping ----------------------------------------- */
    ierr = PetscOptionsInt("-adap", "Adaptive time stepping flag", "", adap, &adap, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-NRmin", "Minimum Newton-Raphson iterations", "", NRmin, &NRmin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-NRmax", "Maximum Newton-Raphson iterations", "", NRmax, &NRmax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-factor", "Time step adjustment factor", "", factor, &factor, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-dtmin", "Minimum time step size", "", dtmin, &dtmin, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-dtmax", "Maximum time step size", "", dtmax, &dtmax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-max_rej", "Maximum number of rejected steps", "", max_rej, &max_rej, NULL); CHKERRQ(ierr);

    /* --- Restart / initialization files --------------------------------- */
    ierr = PetscOptionsString("-initial_cond", "Load initial solution from file", "", initial, initial, sizeof(initial), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-initial_PFgeom", "Load initial ice geometry from file", "", PFgeom, PFgeom, sizeof(PFgeom), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-geom_file",
             "Load an igakit-generated IGA geometry (.dat) via IGARead, "
             "overriding -p/-C/-Nx/-Ny/-Nz axis setup with the geometry's own",
             "", geom_file, geom_file, sizeof(geom_file), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-ic_type",
             "Initial condition geometry (two_ice_grains_boundary|ice_slab|single_ice|multi_grains)",
             "permafrost2.c", ic_type, ic_type, sizeof(ic_type),
             NULL); CHKERRQ(ierr);

    /* --- Capillarly neck parameters ------------------------------------- */
    ierr = PetscOptionsReal("-R1", "Radius of capillary neck", "", user.R1, &user.R1, NULL); CHKERRQ(ierr);

    /* --- Penalty parameters --------------------------------------------- */
    ierr = PetscOptionsReal("-Lambda",
             "Triple-junction penalty strength in the free energy "
             "(larger values suppress spurious phases at binary interfaces; default 1.0)",
             "", user.Lambd, &user.Lambd, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-phase_lo",
             "Lower bound for phase fields phi_ice, phi_air "
             "(simulation aborts if any phi falls below this; default -0.05)",
             "", user.phase_lo, &user.phase_lo, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-phase_hi",
             "Upper bound for phase fields phi_ice, phi_air "
             "(simulation aborts if any phi exceeds this; default 1.05)",
             "", user.phase_hi, &user.phase_hi, NULL); CHKERRQ(ierr);

    PetscOptionsEnd();

    /* Assign parameters to user context */
    user.p = p;
    user.C = C;
    user.dim = dim;
    user.dof = dof;
    user.Nx = Nx;
    user.Ny = Ny;
    user.Nz = Nz;
    user.Lx = Lx;
    user.Ly = Ly;
    user.Lz = Lz;
    user.eps = eps;
    user.temp0 = temp;
    user.grad_temp0[0] = grad_temp0[0];
    user.grad_temp0[1] = grad_temp0[1];
    user.grad_temp0[2] = grad_temp0[2];
    user.hum0 = humidity;
    user.npoints = Nx * Ny * Nz; /* Total number of grid points (for allocating arrays in user context) */
    PetscStrncpy(user.initial_cond, initial, PETSC_MAX_PATH_LEN);
    PetscStrncpy(user.initial_PFgeom, PFgeom, PETSC_MAX_PATH_LEN);

    // Print simulation parameters
    PetscPrintf(PETSC_COMM_WORLD, "----- Simulation Parameters ----- \n");
    PetscPrintf(PETSC_COMM_WORLD, "Dimension:                    %dD \n", dim);
    PetscPrintf(PETSC_COMM_WORLD, "Domain size:                  Lx = %.2e m", Lx);
    if (dim >= 2) PetscPrintf(PETSC_COMM_WORLD, ", Ly = %.2e m", Ly);
    if (dim == 3) PetscPrintf(PETSC_COMM_WORLD, ", Lz = %.2e m", Lz);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    /* Elements/degree/continuity are printed after IGASetUp() below, where
     * the true values are known (geom_file overrides CLI -Nx/-Ny/-p/-C). */
    PetscPrintf(PETSC_COMM_WORLD, "Time stepping:                delt_t = %.2e s, t_final = %.2e s, n_out = %d \n", delt_t, t_final, n_out);
    PetscPrintf(PETSC_COMM_WORLD, "Initial conditions:           T = %.2f °C, humidity = %.2f \n", temp, humidity);
    PetscPrintf(PETSC_COMM_WORLD, "Interface width:              eps = %.2e m \n", eps);
    PetscPrintf(PETSC_COMM_WORLD, "Temperature gradient:         dT/dx = %.2e °C/m", grad_temp0[0]);
    if (dim >= 2) PetscPrintf(PETSC_COMM_WORLD, ", dT/dy = %.2e °C/m", grad_temp0[1]);
    if (dim == 3) PetscPrintf(PETSC_COMM_WORLD, ", dT/dz = %.2e °C/m", grad_temp0[2]);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Vapor diffusivity:            D_v = %.5e m²/s \n", user.dif_vap);
    PetscPrintf(PETSC_COMM_WORLD, "Lambda:                       %.2e \n", user.Lambd);

    /* Compute saturation vapor density and its derivative based on initial temperature */
    PetscReal rho_rhovs;   /* Ratio of ice density to saturation vapor density */
    PetscReal rhoI_vs;     /* Saturation vapor density over ice (at given temperature) */
    PetscReal d_rhovs;     /* Derivative of saturation vapor density with respect to temperature */
    RhoVS_I(&user, user.temp0, &rhoI_vs, &d_rhovs);
    rho_rhovs = user.rho_ice / rhoI_vs; /* Compute ratio */

    /* Adjust boundary condition flags for periodic case */
    if (user.periodic == 1 && flag_BC_Tfix)    flag_BC_Tfix    = PETSC_FALSE;
    if (user.periodic == 1 && flag_BC_rhovfix) flag_BC_rhovfix = PETSC_FALSE;

    /* Resolved boundary-condition summary. Without IGASetBoundaryValue() and
     * with the `if (pnt->atboundary) return 0;` guards in Residual_A1, the
     * natural Neumann BC ∂u/∂n = 0 is what's enforced for any DOF that's not
     * explicitly Dirichlet-fixed. The ice phase field is always left at
     * natural Neumann (no Dirichlet for it anywhere in this code). */
    PetscPrintf(PETSC_COMM_WORLD, "----- Boundary conditions ----- \n");
    PetscPrintf(PETSC_COMM_WORLD, "  ice (phi_i):   natural Neumann (zero gradient)\n");
    PetscPrintf(PETSC_COMM_WORLD, "  temperature:   %s\n",
                flag_BC_Tfix    ? "Dirichlet (fixed value at boundary)"
                                : "natural Neumann (insulating, zero heat flux)");
    PetscPrintf(PETSC_COMM_WORLD, "  vapor (rho_v): %s\n",
                flag_BC_rhovfix ? "Dirichlet (fixed value at boundary)"
                                : "natural Neumann (insulating, zero vapor flux)");
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    /* Time stepping parameters */
    if (n_out > 1) {
        user.t_interv = t_final / (n_out - 1); /* Output interval */
    } else {
        user.t_interv = t_final;
    }

    if (dtmin <= 0.0) dtmin = 0.01 * delt_t;
    if (dtmax <= 0.0) dtmax = 0.5 * user.t_interv;

    /* If dtmax > 0.5*t_interv, print error message */
    if (dtmax > 0.5 * user.t_interv) {
        PetscPrintf(PETSC_COMM_WORLD, "OUTPUT DATA ERROR: Reduce maximum time step, or increase t_interval \n\n");
    }

    /* If adaptive time stepping is enabled, print parameters */
    if (adap == 1) {
        PetscPrintf(PETSC_COMM_WORLD, "Adaptive time stepping is ON \n");
        PetscPrintf(PETSC_COMM_WORLD, "    NRmin = %d, NRmax = %d, factor = %.4f \n", NRmin, NRmax, factor);
        PetscPrintf(PETSC_COMM_WORLD, "    dtmin = %.2e sec, dtmax = %.2e sec \n\n", dtmin, dtmax);
    }

    /* Gibbs-Thomson kinetic parameters */
    user.diff_sub = 0.5 * (user.thcond_air / user.rho_air / user.cp_air + user.thcond_ice / user.rho_ice / user.cp_ice);

    user.Etai = Sigma_i; /* Ice surface energy in the double-well free energy */
    user.Etaa = Sigma_a; /* Air surface energy in the double-well free energy */

    /* Allow CLI override of the physical attachment-kinetics coefficient
     * beta_sub0 via -beta_sub0 <value> (default 1.4e5, set above). Unlike
     * overriding -alph_sub directly (which only rescales the phase-change
     * source term and leaves mob_sub, the AC interfacial-relaxation rate,
     * untouched), beta_sub0 feeds into tau_sub BEFORE both alph_sub and
     * mob_sub are derived below, so a smaller beta_sub0 (faster physical
     * deposition/sublimation kinetics) raises alph_sub and mob_sub together
     * and preserves alph_sub/mob_sub = 3*lambda_sub/eps exactly (that ratio
     * doesn't depend on beta_sub0 at all -- only the overall rate does).
     * That keeps the phase-change source and the interface-profile
     * relaxation in the same proportion that the Karma-Plapp matched
     * asymptotics calibrated, so speeding up sublimation this way should
     * not distort the equilibrium diffuse profile the way overshooting
     * d0_GT did (job64420343/64420350: 1e-7 fattened the profile because
     * that term only enters the curvature-dependent forcing, decoupled from
     * mob_sub's restoring strength). */
    {
        PetscReal  beta_sub0_cli = -1.0;
        PetscBool  flg           = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-beta_sub0",
                                   &beta_sub0_cli, &flg); CHKERRQ(ierr);
        if (flg && beta_sub0_cli > 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -beta_sub0 override: %.4e -> %.4e (factor %.2f)\n",
                        user.beta_sub0, beta_sub0_cli, user.beta_sub0 / beta_sub0_cli);
            user.beta_sub0 = beta_sub0_cli;
        }
    }

    PetscReal a1 = 5.0, a2 = 0.1581; /* Constants for GT relation */
    PetscReal d0_sub;                /* Capillary length parameter */
    PetscReal beta_sub;              /* Kinetic coefficient */
    PetscReal lambda_sub;            /* Gibbs-Thomson parameter */
    PetscReal tau_sub;               /* Gibbs-Thomson parameter */

    d0_sub = user.d0_sub0 / rho_rhovs;
    beta_sub = user.beta_sub0 / rho_rhovs;
    lambda_sub = a1 * user.eps / d0_sub;
    tau_sub = user.eps * lambda_sub * (beta_sub / a1 + a2 * user.eps / user.diff_sub + a2 * user.eps / user.dif_vap);
    user.mob_sub = 1 * user.eps / 3.0 / tau_sub; /* Mobility parameter for sublimation */
    user.alph_sub = lambda_sub / tau_sub;  /* Phase change rate parameter, eq.(9) Moure & Fu (2024) SI */

    /* Gibbs-Thomson capillary length for the curvature correction to the
     * saturation vapor density (rhoI_vs_eff = rhoI_vs*(1 + d0_GT*kappa) in
     * assembly.c). Defaults to d0_sub0, the bare microscopic capillary
     * length constant (NOT d0_sub = d0_sub0/rho_rhovs used just above for
     * the kinetic-coefficient derivation -- that rescaling by the huge
     * rho_ice/rho_vs(T) density ratio belongs to the Karma-Plapp asymptotic
     * matching for alph_sub/mob_sub, not to the standalone Kelvin equation,
     * and would make the curvature correction ~1e-6x too small to have any
     * visible effect; d0_sub0 ~ 1e-9 m is close to the literature ice
     * capillary length, ~9.6e-10 m). Without this term every interface
     * shares one flat rhoI_vs(T) and there is no driving force for Ostwald
     * ripening (confirmed directly: job64415277's VI-solver fix removed the
     * bound-violation artifact that was masquerading as curvature-driven
     * sublimation, and ripening stopped). */
    user.d0_GT = user.d0_sub0;

    /* Allow CLI override, e.g. -d0_GT 0 to disable GT for diagnostics while
     * leaving the kinetic-coefficient derivation above untouched. */
    {
        PetscReal  d0_GT_cli = -1.0;
        PetscBool  flg       = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-d0_GT",
                                   &d0_GT_cli, &flg); CHKERRQ(ierr);
        if (flg && d0_GT_cli >= 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -d0_GT override: %.4e -> %.4e%s\n",
                        user.d0_GT, d0_GT_cli,
                        (d0_GT_cli == 0.0) ? " (Gibbs-Thomson DISABLED)" : "");
            user.d0_GT = d0_GT_cli;
        }
    }
    PetscPrintf(PETSC_COMM_WORLD, "d0_GT (capillary length): %.4e m  %s\n",
                user.d0_GT, (user.d0_GT == 0.0) ? "[GT DISABLED]" : "[GT active]");

    /* Allow per-test override of mob_sub via -mob_sub <value>. Tests with
     * very stiff geometries (touching/merging grains in 2D) can reduce
     * mob_sub by ~10x to trade kinetics speed for AC stability. The
     * physical value computed above is the default. */
    {
        PetscReal  mob_sub_cli = -1.0;
        PetscBool  flg         = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-mob_sub",
                                   &mob_sub_cli, &flg); CHKERRQ(ierr);
        if (flg && mob_sub_cli > 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -mob_sub override: %.2e -> %.2e (factor %.1f)\n",
                        user.mob_sub, mob_sub_cli, mob_sub_cli / user.mob_sub);
            user.mob_sub = mob_sub_cli;
        }
    }

    /* Allow CLI override of the phase-change rate alph_sub via -alph_sub <value>.
     * Setting -alph_sub 0 fully decouples the phase-field equations from the
     * vapor and temperature equations (S_sub = 0 everywhere), which is a useful
     * diagnostic for isolating AC dynamics from coupling effects. Negative
     * values are ignored. */
    {
        PetscReal  alph_sub_cli = -1.0;
        PetscBool  flg          = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-alph_sub",
                                   &alph_sub_cli, &flg); CHKERRQ(ierr);
        if (flg && alph_sub_cli >= 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -alph_sub override: %.2e -> %.2e%s\n",
                        user.alph_sub, alph_sub_cli,
                        (alph_sub_cli == 0.0) ? " (phase-change coupling DISABLED)" : "");
            user.alph_sub = alph_sub_cli;
        }
    }

    PetscPrintf(PETSC_COMM_WORLD, "Mobility terms: mob_sub = %.2e m^3/s \n", user.mob_sub);

    if (!user.flag_Tdep) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "FIXED PARAMETERS: tau %.4e  lambda %.4e  M0 %.4e  alpha %.4e \n\n",
                    tau_sub, lambda_sub, user.mob_sub, user.alph_sub);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "TEMPERATURE DEPENDENT G-T PARAMETERS \n\n");
    }

    /* Create IGA and set up problem */
    IGA iga;
    ierr = IGACreate(PETSC_COMM_WORLD, &iga); CHKERRQ(ierr);
    ierr = IGASetDim(iga, dim); CHKERRQ(ierr);
    ierr = IGASetDof(iga, dof); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 0, "phaseice"); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 1, "temperature"); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 2, "vap_density"); CHKERRQ(ierr);

    /* Set up axes: either from a custom igakit geometry (-geom_file), which
     * defines its own dim/degree/knots/control-net and overrides -p/-C/-N*,
     * or from a uniform tensor-product Cartesian patch (default). */
    if (geom_file[0] != '\0') {
        PetscPrintf(PETSC_COMM_WORLD, "Reading IGA geometry from %s\n", geom_file);
        ierr = IGARead(iga, geom_file); CHKERRQ(ierr);
    } else {
        IGAAxis axis0, axis1, axis2;
        ierr = IGAGetAxis(iga, 0, &axis0); CHKERRQ(ierr);
        if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis0, PETSC_TRUE); CHKERRQ(ierr); }
        ierr = IGAAxisSetDegree(axis0, p); CHKERRQ(ierr);
        ierr = IGAAxisInitUniform(axis0, Nx, 0.0, Lx, C); CHKERRQ(ierr);

        if (dim >= 2) {
            ierr = IGAGetAxis(iga, 1, &axis1); CHKERRQ(ierr);
            if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis1, PETSC_TRUE); CHKERRQ(ierr); }
            ierr = IGAAxisSetDegree(axis1, p); CHKERRQ(ierr);
            ierr = IGAAxisInitUniform(axis1, Ny, 0.0, Ly, C); CHKERRQ(ierr);
        }

        if (dim == 3) {
            ierr = IGAGetAxis(iga, 2, &axis2); CHKERRQ(ierr);
            if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis2, PETSC_TRUE); CHKERRQ(ierr); }
            ierr = IGAAxisSetDegree(axis2, p); CHKERRQ(ierr);
            ierr = IGAAxisInitUniform(axis2, Nz, 0.0, Lz, C); CHKERRQ(ierr);
        }
    }

    ierr = IGASetFromOptions(iga); CHKERRQ(ierr);
    ierr = IGASetUp(iga); CHKERRQ(ierr);
    user.iga = iga;

    /* Number of quadrature points on this rank (used for alph and mob arrays).
     * Read the per-axis degree from the IGA itself (rather than trusting the
     * CLI -p) since -geom_file can set a different degree than -p. */
    PetscInt p_axis[3] = {p, p, p};
    for (PetscInt d = 0; d < dim; d++) {
        IGAAxis ax;
        ierr = IGAGetAxis(iga, d, &ax); CHKERRQ(ierr);
        ierr = IGAAxisGetDegree(ax, &p_axis[d]); CHKERRQ(ierr);
    }

    /* When -geom_file is used, IGARead() sets the actual mesh size from the
     * .dat file and overrides the CLI -Nx/-Ny/-Nz values.  Read the true
     * element counts back from the IGA (node_sizes[d] = N_elements + p) so
     * that user.Nx/Ny/Nz, user.npoints, and the printed header are correct. */
    if (geom_file[0] != '\0') {
        Nx = iga->node_sizes[0] - p_axis[0];
        Ny = (dim >= 2) ? iga->node_sizes[1] - p_axis[1] : 1;
        Nz = (dim == 3) ? iga->node_sizes[2] - p_axis[2] : 1;
        user.Nx = Nx;  user.Ny = Ny;  user.Nz = Nz;
        user.npoints = Nx * Ny * Nz;
        p = p_axis[0];
    }

    PetscPrintf(PETSC_COMM_WORLD, "Elements:                     Nx = %d", Nx);
    if (dim >= 2) PetscPrintf(PETSC_COMM_WORLD, ", Ny = %d", Ny);
    if (dim == 3) PetscPrintf(PETSC_COMM_WORLD, ", Nz = %d", Nz);
    if (geom_file[0] != '\0')
        PetscPrintf(PETSC_COMM_WORLD, "  (from -geom_file; CLI -Nx/-Ny ignored)");
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Polynomial order:             p = %d \n", p);
    PetscPrintf(PETSC_COMM_WORLD, "Continuity order:             C = %d \n", C);

    PetscInt nmb;
    if (dim == 1) {
        nmb = iga->elem_width[0] * (p_axis[0] + 1);
    } else if (dim == 2) {
        nmb = iga->elem_width[0] * iga->elem_width[1] * (p_axis[0] + 1) * (p_axis[1] + 1);
    } else {
        nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2]
              * (p_axis[0] + 1) * (p_axis[1] + 1) * (p_axis[2] + 1);
    }
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.alph);    CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.mob);     CHKERRQ(ierr);
    ierr = PetscMemzero(user.alph,    sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    ierr = PetscMemzero(user.mob,     sizeof(PetscReal) * nmb); CHKERRQ(ierr);

    /* Residual and Jacobian setup */
    ierr = IGASetFormIFunction(iga, Residual, &user); CHKERRQ(ierr);
    ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
    // ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);

    /* Boundary conditions (could 'functionalize' this at some point) */
    // Set vapor density BCs
    if (flag_BC_rhovfix) {
        PetscReal rho0_vs;
        RhoVS_I(&user, user.temp0, &rho0_vs, NULL);
        for (PetscInt l = 0; l < dim; l++) {
            for (PetscInt m = 0; m < 2; m++) {
                ierr = IGASetBoundaryValue(iga, l, m, 2, user.hum0 * rho0_vs); CHKERRQ(ierr);
            }
        }
    }

    // Set temperature BCs
    if (flag_BC_Tfix) {
        PetscReal T_BC[3][2] = {{0}};
        PetscReal LL[3]      = {user.Lx, user.Ly, user.Lz};
        for (PetscInt l = 0; l < dim; l++) {
            for (PetscInt m = 0; m < 2; m++) {
                T_BC[l][m] = user.temp0 + (2.0 * m - 1.0) * user.grad_temp0[l] * LL[l] / 2.0;
                ierr = IGASetBoundaryValue(iga, l, m, 1, T_BC[l][m]); CHKERRQ(ierr);
            }
        }
    }

    /* Set up TS */
    TS ts;
    ierr = IGACreateTS(iga, &ts); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, t_final); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, delt_t); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSALPHA); CHKERRQ(ierr);
    ierr = TSAlphaSetRadius(ts, 0.5); CHKERRQ(ierr);
    if (monitor) { ierr = TSMonitorSet(ts, Monitor, &user, NULL); CHKERRQ(ierr); }
    if (output) { ierr = TSMonitorSet(ts, OutputMonitor, &user, NULL); CHKERRQ(ierr); }

    /* Application context — so BoundsRollbackPreStep can fetch user via
     * TSGetApplicationContext to consume deferred bounds-rollback requests. */
    ierr = TSSetApplicationContext(ts, &user); CHKERRQ(ierr);
    ierr = TSSetPreStep(ts, BoundsRollbackPreStep); CHKERRQ(ierr);

    ierr = TSSetFromOptions(ts); CHKERRQ(ierr);

    ts->adap               = adap;
    ts->NRmin              = NRmin;
    ts->NRmax              = NRmax;
    ts->factor             = factor;
    ts->dtmax              = dtmax;
    ts->dtmin              = dtmin;
    ts->max_reject         = max_rej;
    ts->max_snes_failures  = -1;

    /* Set up SNES non-linear convergence test */
    SNES nonlin;
    ierr = TSGetSNES(ts, &nonlin); CHKERRQ(ierr);
    ierr = SNESSetConvergenceTest(nonlin, SNESDOFConvergence, &user, NULL); CHKERRQ(ierr);

    /* Cache the SNES handle on user so Residual() can call
     * SNESSetFunctionDomainError() when a trial iterate has phi out of bounds.
     * That tells SNES the current line-search trial is invalid; line search
     * backtracks and Newton tries a smaller step. Catches dt-induced AC
     * instabilities while they are still recoverable, before the resulting
     * bad state is committed to ts->vec_sol. */
    user.snes = nonlin;

    /* Bound-constrained Newton solve: enforce -0.1 <= ice <= 1.1 directly
     * on the DOF vector via a variational-inequality SNES (-snes_type
     * vinewtonssls in solver.opts). Field values at quadrature points are a
     * convex combination of nearby DOFs, so bounding the DOFs themselves
     * also bounds the field everywhere.
     *
     * TEMPORARY, pre-conference workaround (2026-06-21): strict [0,1] bounds
     * combined with the per-DOF-block SNES convergence test (see
     * snes_convergence.c) trivially satisfy ABS(atol) almost every step --
     * confirmed by A/B testing atol=1e-6 vs 1e-8 (bit-for-bit identical
     * results) and atol=1e-20 (Newton stagnates, never converges; residual
     * floors are below float64-meaningful precision for some DOF blocks).
     * That's a real per-DOF-tolerance-design issue, not a one-line fix.
     * Loosening the hard VI bounds back toward the old soft-bound/rollback
     * regime's tolerance gives the Newton step slack to move at all instead
     * of getting trivially clamped/declared-converged at the strict bound --
     * the same slack that let Ostwald ripening show up before the VI
     * switch.
     *
     * Widened from -0.05/1.05 to -0.1/1.1 (2026-06-21, still same day):
     * even with -dtmax lowered to 1.0e5, every step still converges via
     * trivial ABS(atol) in 1 iteration (confirmed: 1052/1052 steps on
     * job64440694), and small ice-cap features still pulse (shrink several
     * steps, jump back up one step, repeat) and plateau at a small nonzero
     * residual instead of fully sublimating -- a DOF sitting exactly at the
     * -0.05/1.05 bound becomes VI-"active" and its complementarity
     * treatment can pin it there rather than letting it keep evolving
     * toward wherever the (numerically noisy) single-Newton-step dynamics
     * actually want to take it. Widening the bound just gives more room
     * before that pinning kicks in.
     *
     * Tightened back to -0.05/1.05 (2026-06-22): with the smaller -dtmax
     * (1.0e4) and surgical/no-left-grains geometry changes, runs stay
     * well-behaved without needing the extra -0.1/1.1 slack.
     *
     * Technically incorrect (allows larger unphysical excursions outside
     * [0,1]); intentional tradeoff to get a visibly-evolving result for the
     * 2026-06-23 conference. Revisit the real fix (per-DOF-block atol
     * matched to each field's natural scale) after the conference. */
    Vec Xl, Xu;
    ierr = IGACreateVec(iga, &Xl); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &Xu); CHKERRQ(ierr);
    ierr = VecStrideSet(Xl, 0, -0.05);           CHKERRQ(ierr);
    ierr = VecStrideSet(Xu, 0, 1.05);            CHKERRQ(ierr);
    ierr = VecStrideSet(Xl, 1, PETSC_NINFINITY); CHKERRQ(ierr);
    ierr = VecStrideSet(Xu, 1, PETSC_INFINITY);  CHKERRQ(ierr);
    ierr = VecStrideSet(Xl, 2, PETSC_NINFINITY); CHKERRQ(ierr);
    ierr = VecStrideSet(Xu, 2, PETSC_INFINITY);  CHKERRQ(ierr);
    ierr = SNESVISetVariableBounds(nonlin, Xl, Xu); CHKERRQ(ierr);

    /* Create solution vector (ice, temperature, vapor) */
    Vec U;
    ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);
    ierr = VecZeroEntries(U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Setting up initial conditions... \n");

    if (dim == 1) {
        /* --- 1D Initial Conditions — selected by -ic_type ----------------- */
        PetscPrintf(PETSC_COMM_WORLD, "IC type: %s (1D)\n", ic_type);
        if (strcmp(ic_type, "single_ice") == 0) {
            ierr = FormInitialSingleIceGrain1D(iga, U, &user); CHKERRQ(ierr);
        } else {
            /* Default: centered slab or flat interface, variant via -flag_tIC */
            ierr = FormInitialCondition1D(iga, U, &user); CHKERRQ(ierr);
        }

    } else {
        /* --- 2D / 3D Initial Conditions — selected by -ic_type ------------ */
        PetscPrintf(PETSC_COMM_WORLD, "IC type: %s\n", ic_type);

        if (strcmp(ic_type, "two_ice_grains_boundary") == 0) {
            ierr = FormInitialTwoIceGrainsBoundary2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "ice_slab") == 0) {
            ierr = FormInitialIceSlab2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "single_ice") == 0) {
            ierr = FormInitialSingleIceGrain2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "multi_grains") == 0) {
            ierr = FormInitialMultiGrains2D(iga, U, &user); CHKERRQ(ierr);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Unknown -ic_type. Valid: two_ice_grains_boundary ice_slab single_ice multi_grains");
        }
    }

    /* Solve the system */
    ierr = TSSolve(ts, U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Solution completed. \n");

    /* Cleanup Resources */
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = VecDestroy(&Xl); CHKERRQ(ierr);
    ierr = VecDestroy(&Xu); CHKERRQ(ierr);
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.alph); CHKERRQ(ierr);
    ierr = PetscFree(user.mob); CHKERRQ(ierr);
    /* End Timer */
    PetscLogDouble ltim, tim;
    ierr = PetscTime(&ltim); CHKERRQ(ierr);
    tim = ltim - itim;
    PetscPrintf(PETSC_COMM_WORLD, "Setup time %e sec  =  %.2f min \n\n", tim, tim / 60.0);

    /* Finalize PETSc */
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
