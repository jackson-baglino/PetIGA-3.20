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

    user.xi_v       = 1.0e-4; // 1.0e-5;   /* Time scaling parameter for vapor */
    user.xi_T       = 1.0e-2; // 1.0e-4;   /* Time scaling parameter for temperature */

    user.Lambd      = 1.0;      /* Model parameter Lambda */
    user.air_lim    = 1.0e-6;   /* Air phase fraction */
    user.n_relax    = 0;            /* AC relaxation steps before physics (0 = none) */
    user.flag_relax = PETSC_FALSE;  /* set after option parsing */
    user.t_sed_freeze    = 1.0;        /* 3-phase duration (s); 0 = start immediately in 2-phase */
    user.flag_sed_frozen = PETSC_FALSE; /* set below from t_sed_freeze */
    user.flag_avenue     = 2;          /* default: Avenue 2 */

    user.lat_sub    = 2.83e6;   /* Latent heat of sublimation */

    user.thcond_ice = 2.29;     /* Thermal conductivity of ice */
    user.thcond_sed = 36.0;     /* Thermal conductivity of metal (UPDATE!) */
    user.thcond_air = 0.02;     /* Thermal conductivity of air */

    user.cp_ice     = 1.96e3;   /* Specific heat capacity of ice */
    user.cp_sed     = 4.86e2;   /* Specific heat capacity of metal (UPDATE!) */
    user.cp_air     = 1.044e3;  /* Specific heat capacity of air */

    user.rho_ice    = 919.0;    /* Density of ice */
    user.rho_sed    = 7753.0;   /* Density of metal (UPDATE!) */
    user.rho_air    = 1.341;    /* Density of air */

    user.dif_vap    = 2.178e-5; /* Vapor diffusivity in air */

    user.T_melt     = 0.0;      /* Melting temperature of ice */

    user.flag_tIC   = 0;              /* IC variant: 0=centered slab, 2=flat interface */
    user.readFlag   = PETSC_FALSE;    /* read initial field from file */
    user.flag_Tdep  = PETSC_FALSE;    /* temperature-dependent material properties */

    /* Penalty parameter defaults — match assembly.c hardcoded values.
     * eps is not yet final here; will be overridden after PetscOptionsEnd. */
    user.difvap_pen = 1.0e-5;  /* factor (dimensionless): D_pen = difvap_pen * difvap */
    user.k_pen      = 1.0e7;   /* vapour interface equilibrium stiffness */
    user.k_sed_pen  = -1.0;    /* sentinel: computed from 1e-4/eps² after options */
    user.phase_lo   = -0.25;   /* lower bound: phi below this → abort */
    user.phase_hi   =  1.25;   /* upper bound: phi above this → abort */
    user.d0_sub0    = 1.0e-9;   /* Parameter d0 for substrate */
    user.beta_sub0  = 1.4e5;    /* Parameter beta for substrate */

    PetscReal gamma_im = 0.033; /* Surface energies for ice-metal interface */
    PetscReal gamma_iv = 0.109; /* Surface energies for ice-vapor interface */
    PetscReal gamma_mv = 0.056; /* Surface energies for metal-vapor interface */

    // PetscReal gamma_im = 0.030; /* Surface energies for ice-metal interface */
    // PetscReal gamma_iv = 0.109; /* Surface energies for ice-vapor interface */
    // PetscReal gamma_mv = 0.300; /* Surface energies for metal-vapor interface */

    /* Define common variables (can be overridden by PETSc options) */
    PetscInt  p   = 1;          /* Polynomial order */
    PetscInt  C   = 0;          /* Global continuity order */

    PetscInt  dof = 4;          /* Degrees of freedom per node (ice, temperature, vapor, sediment) */
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
    user.NCsed       = 30;      /* Number of sediment grains */
    user.RCsed       = 0.2e-4;  /* Mean radius */
    user.RCsed_dev   = 0.55;    /* Std dev of radius */

    user.NCice       = 50;      /* Number of ice grains */
    user.RCice       = 0.3e-4;  /* Mean radius */
    user.RCice_dev   = 0.55;    /* Std dev of radius */

    /* Define boundary condition flags (can be overridden by PETSc options) */
    user.periodic    = 0;       /* Periodic boundary condition flag */
    flag_BC_Tfix     = PETSC_TRUE;  /* fix temperature at boundaries by default */
    flag_BC_rhovfix  = PETSC_FALSE;

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
    char      ic_type[64]                 = "enclosed"; /* IC geometry selector */

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

    /* --- Grain geometry: sediment ---------------------------------------- */
    ierr = PetscOptionsInt("-NCsed", "Number of sediment grains", "", user.NCsed, &user.NCsed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCsed", "Mean radius of sediment grains", "", user.RCsed, &user.RCsed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCsed_dev", "Std dev of radius of sediment grains", "", user.RCsed_dev, &user.RCsed_dev, NULL); CHKERRQ(ierr);

    /* --- Grain geometry: ice --------------------------------------------- */
    ierr = PetscOptionsInt("-NCice", "Number of ice grains", "", user.NCice, &user.NCice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice", "Mean radius of ice grains", "", user.RCice, &user.RCice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice_dev", "Std dev of radius of ice grains", "", user.RCice_dev, &user.RCice_dev, NULL); CHKERRQ(ierr);

    /* --- Per-grain radii and separation (enclosed grain pair IC) ----------- */
    user.RCice0    = user.RCice;
    user.RCice1    = user.RCice;
    user.RCsed0    = user.RCsed;
    user.RCsed1    = user.RCsed;
    user.grain_sep = 0.0;
    ierr = PetscOptionsReal("-RCice0",    "Outer ice radius of grain 0 (enclosed IC)",       "", user.RCice0,    &user.RCice0,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCice1",    "Outer ice radius of grain 1 (enclosed IC)",       "", user.RCice1,    &user.RCice1,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCsed0",    "Sediment core radius of grain 0 (enclosed IC)",   "", user.RCsed0,    &user.RCsed0,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-RCsed1",    "Sediment core radius of grain 1 (enclosed IC)",   "", user.RCsed1,    &user.RCsed1,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-grain_sep", "Air gap between outer ice surfaces (enclosed IC)", "", user.grain_sep, &user.grain_sep, NULL); CHKERRQ(ierr);

    /* --- Boundary conditions & physics flags ----------------------------- */
    ierr = PetscOptionsInt("-periodic", "Periodic boundary condition flag", "", user.periodic, &user.periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_Tfix",    "Fix temperature at boundaries",                    "", flag_BC_Tfix,    &flag_BC_Tfix,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_rhovfix", "Fix vapor density at boundaries",                  "", flag_BC_rhovfix, &flag_BC_rhovfix, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_Tdep",       "Temperature-dependent Gibbs-Thomson parameters",   "", user.flag_Tdep,  &user.flag_Tdep,  NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_tIC", "1D IC variant (0=centered slab, 2=flat interface)", "", user.flag_tIC, &user.flag_tIC, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-n_relax",
             "Allen-Cahn relaxation steps before full physics; phase fields settle with no "
             "temperature/vapor coupling (0 = disabled, default)",
             "", user.n_relax, &user.n_relax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-t_sed_freeze",
             "Duration of 3-phase period (s). Set to 0 to start immediately in 2-phase "
             "(frozen sediment). Set > 0 to run 3-phase until this time, then switch.",
             "", user.t_sed_freeze, &user.t_sed_freeze, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_avenue",
             "Residual formulation: 1=AC+vapor penalty+freeze→zero, 2=AC+vapor penalty+freeze→penalty (default), 3=Cahn-Hilliard (no penalties)",
             "", user.flag_avenue, &user.flag_avenue, NULL); CHKERRQ(ierr);

    /* --- Thermophysical properties --------------------------------------- */
    ierr = PetscOptionsReal("-thcond_ice", "Thermal conductivity of ice", "", user.thcond_ice, &user.thcond_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thcond_sed", "Thermal conductivity of inert phase", "", user.thcond_sed, &user.thcond_sed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-thcond_air", "Thermal conductivity of air", "", user.thcond_air, &user.thcond_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-cp_ice", "Specific heat capacity of ice", "", user.cp_ice, &user.cp_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cp_sed", "Specific heat capacity of inert phase", "", user.cp_sed, &user.cp_sed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-cp_air", "Specific heat capacity of air", "", user.cp_air, &user.cp_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-rho_ice", "Density of ice", "", user.rho_ice, &user.rho_ice, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rho_sed", "Density of inert phase", "", user.rho_sed, &user.rho_sed, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rho_air", "Density of air", "", user.rho_air, &user.rho_air, NULL); CHKERRQ(ierr);

    ierr = PetscOptionsReal("-gamma_im", "Surface energy ice-inert phase", "", gamma_im, &gamma_im, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-gamma_iv", "Surface energy ice-vapor", "", gamma_iv, &gamma_iv, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-gamma_mv", "Surface energy inert phase-vapor", "", gamma_mv, &gamma_mv, NULL); CHKERRQ(ierr);

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
    ierr = PetscOptionsString("-ic_type",
             "Initial condition geometry (enclosed|capillary|layered|"
             "random_enclosed|random_packed|ice_cap)",
             "permafrost2.c", ic_type, ic_type, sizeof(ic_type),
             NULL); CHKERRQ(ierr);

    /* --- Capillarly neck parameters ------------------------------------- */
    ierr = PetscOptionsReal("-R1", "Radius of capillary neck", "", user.R1, &user.R1, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-mob_sed", "Sediment phase mobility (0 = inert)", "", user.mob_sed, &user.mob_sed, NULL); CHKERRQ(ierr);

    /* --- Penalty parameters --------------------------------------------- */
    ierr = PetscOptionsReal("-difvap_pen",
             "Multiplicative factor for vapour diffusivity in air: D_pen = difvap_pen * difvap "
             "(dimensionless; default 1e-5)",
             "", user.difvap_pen, &user.difvap_pen, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-k_pen",
             "Vapour interface equilibrium stiffness "
             "(default: difvap_pen/eps²; -1 = use default)",
             "", user.k_pen, &user.k_pen, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-k_sed_pen",
             "Sediment shape-restoring stiffness "
             "(default: 1e-3/eps²; -1 = use default)",
             "", user.k_sed_pen, &user.k_sed_pen, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Lambda",
             "Triple-junction penalty strength in the free energy "
             "(larger values suppress spurious phases at binary interfaces; default 1.0)",
             "", user.Lambd, &user.Lambd, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-phase_lo",
             "Lower bound for phase fields phi_ice, phi_sed, phi_air "
             "(simulation aborts if any phi falls below this; default -0.25)",
             "", user.phase_lo, &user.phase_lo, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-phase_hi",
             "Upper bound for phase fields phi_ice, phi_sed, phi_air "
             "(simulation aborts if any phi exceeds this; default 1.25)",
             "", user.phase_hi, &user.phase_hi, NULL); CHKERRQ(ierr);

    /* --- Flags ---------------------------------------------------------- */

    PetscOptionsEnd();

    /* Resolve sentinel defaults using the final eps value */
    if (user.k_sed_pen < 0.0) user.k_sed_pen = 1.0e-4 / (eps * eps);

    /* Relaxation: if n_relax > 0, start in AC-only mode; monitoring.c flips flag_relax off at step n_relax */
    user.flag_relax = (user.n_relax > 0) ? PETSC_TRUE : PETSC_FALSE;

    /* Determine initial freeze state from t_sed_freeze:
     *   t_sed_freeze <= 0 → start immediately in 2-phase (sediment frozen from t = 0)
     *   t_sed_freeze >  0 → start in 3-phase; monitoring.c flips flag_sed_frozen at t = t_sed_freeze
     * Note: if relaxation is active, the freeze transition is deferred until after relaxation. */
    user.flag_sed_frozen = (user.t_sed_freeze <= 0.0) ? PETSC_TRUE : PETSC_FALSE;

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
    PetscPrintf(PETSC_COMM_WORLD, "Elements:                     Nx = %d", Nx);
    if (dim >= 2) PetscPrintf(PETSC_COMM_WORLD, ", Ny = %d", Ny);
    if (dim == 3) PetscPrintf(PETSC_COMM_WORLD, ", Nz = %d", Nz);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Time stepping:                delt_t = %.2e s, t_final = %.2e s, n_out = %d \n", delt_t, t_final, n_out);
    PetscPrintf(PETSC_COMM_WORLD, "Initial conditions:           T = %.2f °C, humidity = %.2f \n", temp, humidity);
    PetscPrintf(PETSC_COMM_WORLD, "Interface width:              eps = %.2e m \n", eps);
    PetscPrintf(PETSC_COMM_WORLD, "Temperature gradient:         dT/dx = %.2e °C/m", grad_temp0[0]);
    if (dim >= 2) PetscPrintf(PETSC_COMM_WORLD, ", dT/dy = %.2e °C/m", grad_temp0[1]);
    if (dim == 3) PetscPrintf(PETSC_COMM_WORLD, ", dT/dz = %.2e °C/m", grad_temp0[2]);
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    PetscPrintf(PETSC_COMM_WORLD, "Vapor diffusivity:            D_v = %.5e m²/s \n", user.dif_vap);
    PetscPrintf(PETSC_COMM_WORLD, "k_sed_pen:                    %.2e \n", user.k_sed_pen);

    /* Compute saturation vapor density and its derivative based on initial temperature */
    PetscReal rho_rhovs;   /* Ratio of ice density to saturation vapor density */
    PetscReal rhoI_vs;     /* Saturation vapor density over ice (at given temperature) */
    PetscReal d_rhovs;     /* Derivative of saturation vapor density with respect to temperature */
    RhoVS_I(&user, user.temp0, &rhoI_vs, &d_rhovs);
    rho_rhovs = user.rho_ice / rhoI_vs; /* Compute ratio */

    /* Adjust boundary condition flags for periodic case */
    if (user.periodic == 1 && flag_BC_Tfix)    flag_BC_Tfix    = PETSC_FALSE;
    if (user.periodic == 1 && flag_BC_rhovfix) flag_BC_rhovfix = PETSC_FALSE;

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

    user.Etai = gamma_iv + gamma_im - gamma_mv; /* Surface energy differences for ice */
    user.Etam = gamma_mv + gamma_im - gamma_iv; /* Surface energy differences for metal */
    user.Etaa = gamma_iv + gamma_mv - gamma_im; /* Surface energy differences for air */

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
    user.alph_sub = 10 * lambda_sub / tau_sub;  /* Phase change rate parameter */
    user.mob_sed = user.mob_sub;  /* Sediment is inert by default; override with -mob_sed */
    user.mob_air = user.mob_sub;

    PetscPrintf(PETSC_COMM_WORLD, "Mobility terms: mob_sub = %.2e m^3/s, mob_sed = %.2e m^3/s \n", user.mob_sub, user.mob_sed);

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
    ierr = IGASetFieldName(iga, 3, "sediment"); CHKERRQ(ierr);

    /* Set up axes */
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

    ierr = IGASetFromOptions(iga); CHKERRQ(ierr);
    ierr = IGASetUp(iga); CHKERRQ(ierr);
    user.iga = iga;

    /* Number of quadrature points on this rank (used for Phi_sed and mob arrays) */
    PetscInt nmb;
    if (dim == 1) {
        nmb = iga->elem_width[0] * (p + 1);
    } else if (dim == 2) {
        nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(p + 1);
    } else {
        nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2] * CU(p + 1);
    }
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.alph);    CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.mob);     CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.Phi_sed0); CHKERRQ(ierr);
    ierr = PetscMemzero(user.alph,    sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    ierr = PetscMemzero(user.mob,     sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    ierr = PetscMemzero(user.Phi_sed0, sizeof(PetscReal) * nmb); CHKERRQ(ierr);

    /* Residual and Jacobian setup */
    ierr = IGASetFormIFunction(iga, Residual, &user); CHKERRQ(ierr);
    // ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);
    ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);

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

    /* Create solution vector (ice, temperature, vapor, sediment) */
    Vec U;
    ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);
    ierr = VecZeroEntries(U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Setting up initial conditions... \n");

    if (dim == 1) {
        /* --- 1D Initial Conditions — selected by -ic_type ----------------- */
        PetscPrintf(PETSC_COMM_WORLD, "IC type: %s (1D)\n", ic_type);
        if (strcmp(ic_type, "enclosed") == 0) {
            ierr = FormInitialEnclosed1D(iga, U, &user); CHKERRQ(ierr);
        } else {
            /* Default: centered slab or flat interface, variant via -flag_tIC */
            ierr = FormInitialCondition1D(iga, U, &user); CHKERRQ(ierr);
        }

    } else {
        /* --- 2D / 3D Initial Conditions — selected by -ic_type ------------ */
        PetscPrintf(PETSC_COMM_WORLD, "IC type: %s\n", ic_type);

        if (strcmp(ic_type, "capillary") == 0) {
            ierr = FormIC_grain_ana(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "layered") == 0) {
            ierr = FormInitialLayeredPermafrost2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "enclosed") == 0) {
            ierr = FormInitialEnclosedPermafrost2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "contact_sed") == 0) {
            ierr = FormInitialContactSedPermafrost2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "random_enclosed") == 0) {
            ierr = FormInitialRandomEnclosedPermafrost2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "random_packed") == 0) {
            ierr = FormInitialRandomPackedPermafrost2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "ice_cap") == 0) {
            ierr = FormInitialFlatSedIceCap2D(iga, U, &user); CHKERRQ(ierr);
        } else if (strcmp(ic_type, "ice_slab") == 0) {
            ierr = FormInitialIceSlab2D(iga, U, &user); CHKERRQ(ierr);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Unknown -ic_type. Valid: enclosed contact_sed capillary layered "
                    "random_enclosed random_packed ice_cap ice_slab");
        }
    }

    /* Snapshot the initial sediment field at every quadrature point.
     * Phi_sed0[indGP] is indexed identically to alph[] and mob[] so that
     * the Residual function can look it up with the same indGP expression. */
    {
        Vec            localU;
        const PetscScalar *arrayU;
        IGAElement     element;
        IGAPoint       point;
        PetscScalar   *UU;
        PetscInt       indd = 0;

        ierr = IGAGetLocalVecArray(iga, U, &localU, &arrayU); CHKERRQ(ierr);
        ierr = IGABeginElement(iga, &element);                CHKERRQ(ierr);
        while (IGANextElement(iga, element)) {
            ierr = IGAElementGetValues(element, arrayU, &UU); CHKERRQ(ierr);
            ierr = IGAElementBeginPoint(element, &point);     CHKERRQ(ierr);
            while (IGAElementNextPoint(element, point)) {
                PetscScalar sol[4];
                ierr = IGAPointFormValue(point, UU, &sol[0]); CHKERRQ(ierr);
                user.Phi_sed0[indd] = PetscRealPart(sol[3]);  /* DOF 3 = sediment */
                indd++;
            }
            ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        }
        ierr = IGAEndElement(iga, &element);                          CHKERRQ(ierr);
        ierr = IGARestoreLocalVecArray(iga, U, &localU, &arrayU);    CHKERRQ(ierr);
    }

    /* Solve the system */
    ierr = TSSolve(ts, U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Solution completed. \n");

    /* Cleanup Resources */
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.alph); CHKERRQ(ierr);
    ierr = PetscFree(user.mob); CHKERRQ(ierr);
    ierr = PetscFree(user.Phi_sed0); CHKERRQ(ierr);
    /* End Timer */
    PetscLogDouble ltim, tim;
    ierr = PetscTime(&ltim); CHKERRQ(ierr);
    tim = ltim - itim;
    PetscPrintf(PETSC_COMM_WORLD, "Setup time %e sec  =  %.2f min \n\n", tim, tim / 60.0);

    /* Finalize PETSc */
    ierr = PetscFinalize(); CHKERRQ(ierr);
    return 0;
}
