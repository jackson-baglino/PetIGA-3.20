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
    PetscInt flag_sedgrav;               /* Flag for sediment gravity */
    PetscInt flag_BC_Tfix;               /* Flag for temperature boundary condition */
    PetscInt flag_BC_rhovfix;            /* Flag for fixed rho_v boundary condition */

    user.xi_v       = 1.0e-5;   /* Time scaling parameter for vapor */
    user.xi_T       = 1.0e-4;   /* Time scaling parameter for temperature */
    user.flag_xiT   = 1;        /* Flag for temperature */

    user.Lambd      = 1.0;      /* Model parameter Lambda */
    user.air_lim    = 1.0e-6;   /* Air phase fraction */
    user.nsteps_IC  = 10;       /* Number of initial condition steps (???) */

    user.lat_sub    = 2.83e6;   /* Latent heat of sublimation */

    user.thcond_ice = 2.29;     /* Thermal conductivity of ice */
    user.thcond_met = 36.0;     /* Thermal conductivity of metal (UPDATE!) */
    user.thcond_air = 0.02;     /* Thermal conductivity of air */

    user.cp_ice     = 1.96e3;   /* Specific heat capacity of ice */
    user.cp_met     = 4.86e2;   /* Specific heat capacity of metal (UPDATE!) */
    user.cp_air     = 1.044e3;  /* Specific heat capacity of air */

    user.rho_ice    = 919.0;    /* Density of ice */
    user.rho_met    = 7753.0;   /* Density of metal (UPDATE!) */
    user.rho_air    = 1.341;    /* Density of air */

    user.dif_vap    = 2.178e-5; /* Vapor diffusivity in air */

    user.T_melt     = 0.0;      /* Melting temperature of ice */

    user.flag_it0   = 1;        /* Initialization flag */
    user.flag_tIC   = 0;        /* Initial condition flag */
    user.readFlag   = 0;        /* Flag to read ice grains from file (UPDATE IMPLEMENTATION!) */

    user.flag_Tdep  = 1;        /* Temperature-dependent Gibbs-Thomson parameters */
    user.d0_sub0    = 1.0e-9;   /* Parameter d0 for substrate */
    user.beta_sub0  = 1.4e5;    /* Parameter beta for substrate */

    PetscReal gamma_im = 0.033; /* Surface energies for ice-metal interface */
    PetscReal gamma_iv = 0.109; /* Surface energies for ice-vapor interface */
    PetscReal gamma_mv = 0.056; /* Surface energies for metal-vapor interface */

    /* Define common variables (can be overridden by PETSc options) */
    PetscInt  p   = 1;          /* Polynomial order */
    PetscInt  C   = 0;          /* Global continuity order */

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
    flag_sedgrav     = 0;       /* Sediment gravity flag */
    user.NCsed       = 30;      /* Number of sediment grains */
    user.RCsed       = 0.2e-4;  /* Mean radius */
    user.RCsed_dev   = 0.55;    /* Std dev of radius */

    user.NCice       = 50;      /* Number of ice grains */
    user.RCice       = 0.5e-4;  /* Mean radius */
    user.RCice_dev   = 0.55;    /* Std dev of radius */

    /* Define boundary condition flags (can be overridden by PETSc options) */
    user.periodic    = 0;       /* Periodic boundary condition flag */
    flag_BC_Tfix     = 1;       /* Temperature BC flag */
    flag_BC_rhovfix  = 0;       /* Vapor density BC flag */

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

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Permafrost options", "IGA");
    /* --- Geometry & discretization --------------------------------------- */
    ierr = PetscOptionsInt("-dim", "Problem dimension (2 or 3)", "", dim, &dim, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-Nx", "Number of elements in x direction", "", Nx, &Nx, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-Ny", "Number of elements in y direction", "", Ny, &Ny, NULL); CHKERRQ(ierr);
    if (dim == 3) {
        ierr = PetscOptionsInt("-Nz", "Number of elements in z direction", "", Nz, &Nz, NULL); CHKERRQ(ierr);
    }
    PetscInt ngrad = dim; /* Number of grad_temp0 components to read */
    ierr = PetscOptionsReal("-Lx", "Domain length in x direction", "", Lx, &Lx, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-Ly", "Domain length in y direction", "", Ly, &Ly, NULL); CHKERRQ(ierr);
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

    /* --- Boundary conditions & physics flags ----------------------------- */
    ierr = PetscOptionsInt("-periodic", "Periodic boundary condition flag", "", user.periodic, &user.periodic, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_sedgrav", "Sediment gravity flag", "", flag_sedgrav, &flag_sedgrav, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_BC_Tfix", "Temperature BC flag", "", flag_BC_Tfix, &flag_BC_Tfix, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-flag_BC_rhovfix", "Vapor density BC flag", "", flag_BC_rhovfix, &flag_BC_rhovfix, NULL); CHKERRQ(ierr);

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

    /* --- Capillarly neck parameters ------------------------------------- */
    ierr = PetscOptionsReal("-R1", "Radius of capillary neck", "", user.R1, &user.R1, NULL); CHKERRQ(ierr);

    PetscOptionsEnd();

    /* Assign parameters to user context */
    user.p = p;
    user.C = C;
    user.dim = dim;
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
    PetscStrncpy(user.initial_cond, initial, PETSC_MAX_PATH_LEN);
    PetscStrncpy(user.initial_PFgeom, PFgeom, PETSC_MAX_PATH_LEN);

    // Print simulation parameters
    PetscPrintf(PETSC_COMM_WORLD, "----- Simulation Parameters ----- \n");
    PetscPrintf(PETSC_COMM_WORLD, "Dimension: %dD \n", dim);
    PetscPrintf(PETSC_COMM_WORLD, "Domain size: Lx = %.2e m, Ly = %.2e m", Lx, Ly);
    if (dim == 3) {
        PetscPrintf(PETSC_COMM_WORLD, ", Lz = %.2e m", Lz);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\nElements: Nx = %d, Ny = %d", Nx, Ny);
    if (dim == 3) {
        PetscPrintf(PETSC_COMM_WORLD, ", Nz = %d", Nz);
    }
    PetscPrintf(PETSC_COMM_WORLD, "\nTime step: delt_t = %.2e sec, t_final = %.2e sec, n_out = %d \n", delt_t, t_final, n_out);
    PetscPrintf(PETSC_COMM_WORLD, "Initial conditions: temp0 = %.2f C, humidity = %.2f \n", temp, humidity);
    PetscPrintf(PETSC_COMM_WORLD, "Interface width: eps = %.2e m \n", eps);
    PetscPrintf(PETSC_COMM_WORLD, "Temperature gradient: dT/dx = %.2e C/m, dT/dy = %.2e C/m", grad_temp0[0], grad_temp0[1]);
    if (dim == 3) {
        PetscPrintf(PETSC_COMM_WORLD, ", dT/dz = %.2e C/m", grad_temp0[2]);
    }

    /* Compute saturation vapor density and its derivative based on initial temperature */
    PetscReal rho_rhovs;   /* Ratio of ice density to saturation vapor density */
    PetscReal rhoI_vs;     /* Saturation vapor density over ice (at given temperature) */
    PetscReal d_rhovs;     /* Derivative of saturation vapor density with respect to temperature */
    RhoVS_I(&user, user.temp0, &rhoI_vs, &d_rhovs);
    rho_rhovs = user.rho_ice / rhoI_vs; /* Compute ratio */

    /* Adjust boundary condition flags for periodic case */
    if (user.periodic == 1 && flag_BC_Tfix == 1) flag_BC_Tfix = 0;
    if (user.periodic == 1 && flag_BC_rhovfix == 1) flag_BC_rhovfix = 0;

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

    if (user.flag_Tdep == 0) {
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
    ierr = IGASetDof(iga, 3); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 0, "phaseice"); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 1, "temperature"); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 2, "vap_density"); CHKERRQ(ierr);

    /* Set up axes */
    IGAAxis axis0, axis1, axis2;
    ierr = IGAGetAxis(iga, 0, &axis0); CHKERRQ(ierr);
    if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis0, PETSC_TRUE); CHKERRQ(ierr); }
    ierr = IGAAxisSetDegree(axis0, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis0, Nx, 0.0, Lx, C); CHKERRQ(ierr);

    ierr = IGAGetAxis(iga, 1, &axis1); CHKERRQ(ierr);
    if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis1, PETSC_TRUE); CHKERRQ(ierr); }
    ierr = IGAAxisSetDegree(axis1, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis1, Ny, 0.0, Ly, C); CHKERRQ(ierr);

    if (dim == 3) {
        ierr = IGAGetAxis(iga, 2, &axis2); CHKERRQ(ierr);
        if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis2, PETSC_TRUE); CHKERRQ(ierr); }
        ierr = IGAAxisSetDegree(axis2, p); CHKERRQ(ierr);
        ierr = IGAAxisInitUniform(axis2, Nz, 0.0, Lz, C); CHKERRQ(ierr);
    }

    ierr = IGASetFromOptions(iga); CHKERRQ(ierr);
    ierr = IGASetUp(iga); CHKERRQ(ierr);
    user.iga = iga;

    /* Setup Initial Conditions */
    PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(p + 1);
    if (dim == 3) {
        nmb = nmb * iga->elem_width[2] * (p + 1);
    }
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.Phi_sed); CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.alph);    CHKERRQ(ierr);
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.mob);     CHKERRQ(ierr);
    ierr = PetscMemzero(user.Phi_sed, sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    ierr = PetscMemzero(user.alph,    sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    ierr = PetscMemzero(user.mob,     sizeof(PetscReal) * nmb); CHKERRQ(ierr);

    /* Residual and Jacobian setup */
    ierr = IGASetFormIFunction(iga, Residual, &user); CHKERRQ(ierr);
    ierr = IGASetFormIJacobian(iga, Jacobian, &user); CHKERRQ(ierr);

    /* Boundary conditions (could 'functionalize' this at some point) */
    // Set vapor density BCs
    if (flag_BC_rhovfix == 1) {
        PetscReal rho0_vs;
        RhoVS_I(&user, user.temp0, &rho0_vs, NULL);
        for (PetscInt l = 0; l < dim; l++) {
            for (PetscInt m = 0; m < 2; m++) {
                ierr = IGASetBoundaryValue(iga, l, m, 2, user.hum0 * rho0_vs); CHKERRQ(ierr);
            }
        }
    }

    // Set temperature BCs
    if (flag_BC_Tfix == 1) {
        PetscReal T_BC[dim][2], LL[dim];
        LL[0] = user.Lx;
        LL[1] = user.Ly;
        if (dim == 3) {LL[2] = user.Lz; }
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

    /* Setup sediment phase */
    IGA igaS;
    ierr = IGACreate(PETSC_COMM_WORLD, &igaS); CHKERRQ(ierr);
    ierr = IGASetDim(igaS, dim); CHKERRQ(ierr);
    ierr = IGASetDof(igaS, 1); CHKERRQ(ierr);
    // ierr = IGASetFieldName(igaS, 0, "Sediment Phase"); CHKERRQ(ierr);

    // Set up axes for sediment phase
    IGAAxis axis0S, axis1S, axis2S;
    ierr = IGAGetAxis(igaS, 0, &axis0S);
    if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis0S, PETSC_TRUE); CHKERRQ(ierr); }
    ierr = IGAAxisSetDegree(axis0S, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis0S, Nx, 0.0, Lx, C); CHKERRQ(ierr);
    ierr = IGAGetAxis(igaS, 1, &axis1S);
    if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis1S, PETSC_TRUE); CHKERRQ(ierr); }
    ierr = IGAAxisSetDegree(axis1S, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis1S, Ny, 0.0, Ly, C); CHKERRQ(ierr);
    if (dim == 3) {
        ierr = IGAGetAxis(igaS, 2, &axis2S);
        if (user.periodic == 1) { ierr = IGAAxisSetPeriodic(axis2S, PETSC_TRUE); CHKERRQ(ierr); }
        ierr = IGAAxisSetDegree(axis2S, p); CHKERRQ(ierr);
        ierr = IGAAxisInitUniform(axis2S, Nz, 0.0, Lz, C); CHKERRQ(ierr);
    }
    ierr = IGASetFromOptions(igaS); CHKERRQ(ierr);
    ierr = IGASetUp(igaS); CHKERRQ(ierr);

    /* Create solution vectors for ice, temp, and air (U) and sediment phase (S) */
    Vec U, S;
    ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);
    ierr = IGACreateVec(igaS, &S); CHKERRQ(ierr);
    ierr = VecZeroEntries(U); CHKERRQ(ierr);
    ierr = VecZeroEntries(S); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Setting up initial conditions... \n");
    // PetscPrintf(PETSC_COMM_WORLD, "Initial condition type: %s \n", ICGeomTypeNames[ic_type_opt]);

    // PetscPrintf(PETSC_COMM_WORLD,
    //         "IC type: capillary  (using analytic capillary neck geometry)\n");
    // ierr = FormIC_grain_ana(iga, U, igaS, S, &user); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "IC type: layered  (using layered geometry)\n");
    ierr = FormInitialLayeredPermafrost2D(iga, igaS, U, S, &user); CHKERRQ(ierr);

    /* Optional statements for varying input types (e.g., capillary, layered, etc.) */

    /* Write initial sediment phase to file */
    char filename_sed[PETSC_MAX_PATH_LEN], filevect_sed[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(filename_sed, sizeof(filename_sed), "%s/igasoil.dat", user.output_path); CHKERRQ(ierr);
    ierr = IGAWrite(igaS, filename_sed); CHKERRQ(ierr);

    ierr = PetscSNPrintf(filevect_sed, sizeof(filevect_sed), "%s/soil.dat", user.output_path); CHKERRQ(ierr);
    ierr = IGAWriteVec(igaS, S, filevect_sed); CHKERRQ(ierr);

    /* Solve the system */
    ierr = TSSolve(ts, U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Solution completed. \n");

    /* Cleanup Resources */
    // Destroy vectors and TS
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = VecDestroy(&S); CHKERRQ(ierr);
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = IGADestroy(&igaS); CHKERRQ(ierr);

    // Free allocated memory for sediment phase and parameters
    ierr = PetscFree(user.Phi_sed); CHKERRQ(ierr);
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
