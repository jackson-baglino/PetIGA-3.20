#include "NASA_main.h" // Need to change name later

/* SNESSetFunctionDomainError() is logically collective, but Residual()
 * (assembly.c) flags it per quadrature point — i.e. only on the rank that
 * owns the offending point. PETSc's line searches branch on the raw flag
 * BEFORE any norm collective (SNESLineSearchApply_Basic returns early with
 * SNES_LINESEARCH_FAILED_DOMAIN), so a rank-local flag splits the ranks'
 * control flow and deadlocks the run: observed 2026-07-09 on the 6-rank
 * two-grain case at step 52, one rank stuck in MPI_Barrier inside
 * SNESNEWTONLSCheckLocalMin_Private, the rest in VecNormEnd/MPI_Allreduce
 * inside the line search. (PETSc's own VecSetInf mitigation in
 * SNESComputeFunction only helps paths that compute a norm before
 * branching, which the basic line search does not.)
 *
 * Wrap the TS-provided SNES residual so every evaluation ends with an
 * allreduce that makes the flag uniform across ranks; the wrapper runs
 * inside SNESComputeFunction, so its VecSetInf tagging then also fires
 * consistently on all ranks. */
static PetscErrorCode SNESTSFormFunction_DomainErrSync(SNES snes, Vec X, Vec F, void *ctx)
{
    PetscErrorCode ierr;
    PetscBool derr_loc = PETSC_FALSE, derr_glob = PETSC_FALSE;
    PetscFunctionBegin;
    ierr = SNESTSFormFunction(snes, X, F, ctx); CHKERRQ(ierr);
    ierr = SNESGetFunctionDomainError(snes, &derr_loc); CHKERRQ(ierr);
    ierr = MPIU_Allreduce(&derr_loc, &derr_glob, 1, MPIU_BOOL, MPI_LOR,
                          PetscObjectComm((PetscObject)snes)); CHKERRQ(ierr);
    if (derr_glob && !derr_loc) { ierr = SNESSetFunctionDomainError(snes); CHKERRQ(ierr); }
    PetscFunctionReturn(0);
}

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
    user.rhovfix_axis = -1;     /* -1 = pin vapor on every face (legacy) */

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

    /* Interface-CFL timestep limiter (InterfaceCFLMonitor) */
    user.flag_dtCFL   = PETSC_TRUE;   /* on by default */
    user.cfl_dphimax  = 0.2;          /* max pointwise |dphi| per step */
    user.cfl_U_prev   = NULL;
    user.cfl_t_prev   = 0.0;

    user.axisym = PETSC_FALSE;        /* axisymmetric r-z mode (see NASA_types.h) */
    user.ic_grain_union = PETSC_FALSE; /* multi_grains IC: additive (see NASA_types.h) */
    user.ssa_view = NULL;              /* SSA_evo.dat viewer, opened lazily in Monitor() */
    user.decouple_phase_change = PETSC_FALSE;  /* see NASA_types.h / assembly.c */

    user.phase_lo   = -0.05;   /* lower bound: phi below this → abort */
    user.phase_hi   =  1.05;   /* upper bound: phi above this → abort */

    /* Temporal-scaling factors (M&F 2024 §3.1, eqs. 25-26): slow the fast
     * T / vapor diffusion timescales by 1/xi while keeping the quasi-steady
     * fields (and thus interface velocity) physical, permitting large dt.
     * Each xi must scale the diffusion term AND its phase-change source
     * together so xi cancels in the quasi-steady balance. M&F values:
     * xi_T = 1 (no solidification), xi_v = 1e-3. */
    user.xi_T = 1.0;    /* thermal: scales conduction + latent heat in R_tem  */
    user.xi_v = 1e-3;   /* vapor:   scales diffusion + rho_ice source in R_vap */

    user.d0_sub0    = 1e-7; // 9.6e-10;   /* capillary length d0 = gamma*Vm/(R*T) at -5°C [m] */
    user.beta_sub0  = 9.9e5;     /* beta0 = (1/alpha_c)*sqrt(2pi*m/kT)/(rho_vs/rho_i)
                                  * at alpha_c=2e-3 (Libbrecht 2017), T=-5°C [s/m] */

    /* Surface energy parameters of the double-well free energy [J/m²]:
     *   F_dub(phi_i) = C*phi_i^2(1-phi_i)^2,  C = (Sigma_i+Sigma_a)/2 + Lambda
     * Sigma_i = ice-side surface energy, Sigma_a = air-side surface energy. */
    PetscReal Sigma_i = 0.109; /* ice surface energy [J/m²] */
    PetscReal Sigma_a = 0.132; /* air surface energy [J/m²] */

    /* Prescribed contact angle at the regolith wall. The substrate is the
     * domain boundary in this two-phase model, so its wetting behaviour is set
     * by the three surface energies through Young's equation,
     *     cos(theta) = (gamma_as - gamma_is)/gamma_ia.
     * gamma_is = gamma_as (the default) gives cos(theta) = 0, i.e. theta = 90
     * deg, which is exactly the natural Neumann wall the solver has always had
     * -- so leaving these alone changes nothing. gamma_ia defaults to Sigma_i
     * below, once -Sigma_i has been read. */
    PetscReal gamma_ia = -1.0;   /* < 0 => "unset", take Sigma_i */
    PetscReal gamma_is = 0.0;
    PetscReal gamma_as = 0.0;
    PetscReal contact_angle_deg = 0.0;
    PetscBool contact_angle_set = PETSC_FALSE;
    char      wall_faces[64] = "";
    PetscBool test_wall_measure = PETSC_FALSE;
    PetscBool test_wall_jacobian = PETSC_FALSE;

    /* Define common variables (can be overridden by PETSc options) */
    PetscInt  p   = 2;          /* Polynomial order */
    PetscInt  C   = 1;          /* Global continuity order */

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
    PetscReal t_start = 0.0;    /* Clock value to resume from (-initial_cond restarts) */

    PetscInt  n_out   = 10;     /* Number of outputs */

    PetscReal humidity = 0.95;  /* Initial humidity */
    PetscReal temp     = -20.0; /* Initial temperature */

    /* Temperature the mesh/eps were sized for by comp_eps.py (set by generated
     * geometry .opts). Guards against running at a -temp inconsistent with the
     * eps/mesh resolution (see the check after PetscOptionsEnd). */
    PetscReal eps_valid_temp      = 0.0;
    PetscBool eps_valid_temp_set  = PETSC_FALSE;
    PetscBool eps_temp_override   = PETSC_FALSE;

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

    /* Stall detector: 500 consecutive bit-identical ||U|| values. Generous
     * enough that no genuine transient trips it, small enough that a frozen run
     * dies in minutes instead of burning its whole wall-clock allocation. */
    user.stall_norm  = -1.0;
    user.stall_count = 0;
    user.stall_limit = 500;
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

    PetscOptionsBegin(PETSC_COMM_WORLD, "", "Phase-field options", "IGA");
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
    ierr = PetscOptionsReal("-eps_valid_temp",
                            "Temperature [C] the mesh/eps were sized for (comp_eps.py)",
                            "", eps_valid_temp, &eps_valid_temp, &eps_valid_temp_set); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-eps_temp_override",
                            "Proceed even if -temp disagrees with -eps_valid_temp",
                            "", eps_temp_override, &eps_temp_override, NULL); CHKERRQ(ierr);

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

    /* --- Affine wall baselines (tapered / wedge domains) ------------------- *
     * The bumps above have compact support and cannot express a linear ramp,
     * so a wall that rises or falls across the whole domain needs these. They
     * MUST match build_geometry_multi_grain.py's --bot-y0/--bot-slope/
     * --top-y0/--top-slope: the mesh is cut from the same two curves, and the
     * IC is seeded by mapping (u,v) through them, so a mismatch puts the ice
     * where the mesh is not. Defaults reproduce a flat [0,Ly] channel. */
    user.wall_bot_y0    = 0.0;
    user.wall_bot_slope = 0.0;
    user.wall_top_y0    = Ly;
    user.wall_top_slope = 0.0;
    ierr = PetscOptionsReal("-wall_bot_y0",
             "Bottom wall height at x=0 [m] (affine baseline under the bumps)",
             "", user.wall_bot_y0, &user.wall_bot_y0, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wall_bot_slope",
             "Bottom wall dy/dx [-]; negative = floor drops to the right",
             "", user.wall_bot_slope, &user.wall_bot_slope, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wall_top_y0",
             "Top wall height at x=0 [m] (defaults to Ly)",
             "", user.wall_top_y0, &user.wall_top_y0, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wall_top_slope",
             "Top wall dy/dx [-]; positive = ceiling rises to the right",
             "", user.wall_top_slope, &user.wall_top_slope, NULL); CHKERRQ(ierr);

    /* --- Wedge-bridging ice band ------------------------------------------ *
     * Annulus about the wedge apex. An apex-centred arc is perpendicular to
     * every ray from the apex, hence to both wedge walls, so this meets both
     * at the natural 90-degree contact angle -- which a circle cannot do. */
    user.wedge_apex_x   = 0.0;
    user.wedge_apex_y   = 0.0;
    user.n_wedge_bands  = 0;
    ierr = PetscOptionsReal("-wedge_apex_x", "Wedge apex x [m] (virtual, outside the domain)",
             "", user.wedge_apex_x, &user.wedge_apex_x, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-wedge_apex_y", "Wedge apex y [m]",
             "", user.wedge_apex_y, &user.wedge_apex_y, NULL); CHKERRQ(ierr);
    {
        PetscInt n1 = MAX_SED_GRAINS, n2 = MAX_SED_GRAINS;
        PetscBool f1;
        ierr = PetscOptionsRealArray("-wedge_band_r1",
                 "Inner radii of the ice bands from the apex [m], one per band",
                 "", user.wedge_band_r1, &n1, &f1); CHKERRQ(ierr);
        if (f1) {
            user.n_wedge_bands = n1;
            ierr = PetscOptionsRealArray("-wedge_band_r2",
                     "Outer radii of the ice bands from the apex [m], one per band",
                     "", user.wedge_band_r2, &n2, NULL); CHKERRQ(ierr);
            if (n2 != n1)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                        "-wedge_band_r1 and -wedge_band_r2 must have the same "
                        "length (%d vs %d)", (int)n1, (int)n2);
            for (PetscInt k = 0; k < n1; k++) {
                if (user.wedge_band_r1[k] < 0.0)
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                            "-wedge_band_r1[%d] must be >= 0 (got %g)",
                            (int)k, (double)user.wedge_band_r1[k]);
                if (user.wedge_band_r2[k] <= user.wedge_band_r1[k])
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                            "band %d: r2 (%g) must exceed r1 (%g)", (int)k,
                            (double)user.wedge_band_r2[k], (double)user.wedge_band_r1[k]);
                if (k > 0 && user.wedge_band_r1[k] <= user.wedge_band_r2[k-1])
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                            "bands %d and %d overlap (r1[%d]=%g <= r2[%d]=%g); "
                            "they must be listed in increasing radius and be "
                            "disjoint", (int)k-1, (int)k, (int)k,
                            (double)user.wedge_band_r1[k], (int)k-1,
                            (double)user.wedge_band_r2[k-1]);
            }
        }
    }

    /* --- Frozen capillary bridges spanning a pore throat ------------------- *
     * Region outside two axis-centred meniscus circles and between their
     * centres, giving CONCAVE (waisted) menisci at 90 degrees to both walls.
     * The radii must be solved against the actual wall shape -- generate these
     * with preprocess/build_geometry_two_throat.py, do not hand-write them. */
    user.n_bridges = 0;
    {
        PetscInt n = MAX_SED_GRAINS;
        ierr = PetscOptionsRealArray("-bridge_cxL",
                 "Capillary-bridge LEFT meniscus circle centres, x [m]; one per bridge",
                 "", user.bridge_cxL, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            user.n_bridges = n;
            struct { const char *name; PetscReal *arr; const char *desc; } req[] = {
                {"-bridge_rL",  user.bridge_rL,  "LEFT meniscus radii [m]"},
                {"-bridge_cxR", user.bridge_cxR, "RIGHT meniscus circle centres, x [m]"},
                {"-bridge_rR",  user.bridge_rR,  "RIGHT meniscus radii [m]"},
                {"-bridge_cy",  user.bridge_cy,  "channel-axis y of each bridge [m]"},
            };
            for (PetscInt q = 0; q < 4; q++) {
                PetscInt nq = MAX_SED_GRAINS;
                ierr = PetscOptionsRealArray(req[q].name, req[q].desc, "",
                         req[q].arr, &nq, NULL); CHKERRQ(ierr);
                if (nq != n)
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                            "-bridge_cxL and %s must have the same length (%d vs %d)",
                            req[q].name, (int)n, (int)nq);
            }
            for (PetscInt k = 0; k < n; k++) {
                if (user.bridge_rL[k] <= 0.0 || user.bridge_rR[k] <= 0.0)
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                            "bridge %d: meniscus radii must be > 0", (int)k);
                /* Intersection of the two discs must be non-empty: their
                 * centres must be closer together than the sum of the radii. */
                if (PetscAbsReal(user.bridge_cxR[k] - user.bridge_cxL[k]) >=
                    user.bridge_rL[k] + user.bridge_rR[k])
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                            "bridge %d: meniscus discs do not overlap (|dcx|=%g >= "
                            "rL+rR=%g) -- the plug would be empty", (int)k,
                            (double)PetscAbsReal(user.bridge_cxR[k] - user.bridge_cxL[k]),
                            (double)(user.bridge_rL[k] + user.bridge_rR[k]));
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

    /* --- Ice-grain array capacity ------------------------------------------ */
    user.n_grain_max = 2000;
    ierr = PetscOptionsInt("-n_grain_max",
             "Capacity of the ice-grain centre/radius arrays", "",
             user.n_grain_max, &user.n_grain_max, NULL); CHKERRQ(ierr);
    if (user.n_grain_max < 1)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "-n_grain_max must be >= 1 (got %d)", (int)user.n_grain_max);
    for (PetscInt d = 0; d < 3; d++) {
        ierr = PetscCalloc1(user.n_grain_max, &user.cent[d]); CHKERRQ(ierr);
    }
    ierr = PetscCalloc1(user.n_grain_max, &user.radius);       CHKERRQ(ierr);
    ierr = PetscCalloc1(user.n_grain_max, &user.ice_grain_ax); CHKERRQ(ierr);
    ierr = PetscCalloc1(user.n_grain_max, &user.ice_grain_ay); CHKERRQ(ierr);

    /* --- Multi-grain ice IC (-ic_type multi_grains) ------------------------ */
    user.n_act = 0;
    {
        PetscInt   n = user.n_grain_max;
        PetscReal *cx;
        ierr = PetscMalloc1(user.n_grain_max, &cx); CHKERRQ(ierr);
        ierr = PetscOptionsRealArray("-ice_grain_cx",
                 "Ice grain center x-positions [m] (-ic_type multi_grains)",
                 "", cx, &n, &flg); CHKERRQ(ierr);
        if (flg) {
            PetscInt   ncy = user.n_grain_max, nr = user.n_grain_max;
            PetscReal *cy, *rr;
            ierr = PetscMalloc1(user.n_grain_max, &cy); CHKERRQ(ierr);
            ierr = PetscMalloc1(user.n_grain_max, &rr); CHKERRQ(ierr);
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
            ierr = PetscFree(cy); CHKERRQ(ierr);
            ierr = PetscFree(rr); CHKERRQ(ierr);
        }
        ierr = PetscFree(cx); CHKERRQ(ierr);
    }

    /* --- Multi-grain ice IC from file (-ic_type multi_grains_file) ---------
     * Whitespace-delimited grain list, one grain per line, produced by
     * preprocess/generate_packing.py:
     *
     *     Lx Ly            <- optional 2-token header (domain size, metres)
     *     x  y  r          <- one row per grain, metres
     *
     * A 4-token row (x y z r) from the legacy 3D-flavoured generators is
     * also accepted: z is read and discarded, since the model is 2D. Mixing
     * row widths inside one file is rejected rather than guessed at -- the
     * legacy reader inferred the width from the first line only and left the
     * radius uninitialised on short rows, which silently produced garbage
     * grains. Parsed on rank 0 and broadcast; every rank must agree. */
    ierr = PetscOptionsString("-grains_file",
             "Path to a grain list for -ic_type multi_grains_file", "",
             user.grains_file, user.grains_file,
             sizeof(user.grains_file), &flg); CHKERRQ(ierr);

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
    user.thin_iface_corr = PETSC_TRUE;   /* ON by default since 2026-09-13 */
    ierr = PetscOptionsBool("-thin_iface_corr",
             "Include the Karma thin-interface counter-terms in tau_sub. Default 0: "
             "for a one-sided vapour diffusivity the O(eps) kinetic contribution they "
             "compensate is identically zero, so including them inflates the realised "
             "beta above -beta_sub0 (1.22x at the -20 C wedge parameters). See "
             "docs/gt_deficit/",
             "", user.thin_iface_corr, &user.thin_iface_corr, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_Tfix",    "Fix temperature at boundaries",                    "", flag_BC_Tfix,    &flag_BC_Tfix,    NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-flag_BC_rhovfix", "Fix vapor density at boundaries",                  "", flag_BC_rhovfix, &flag_BC_rhovfix, NULL); CHKERRQ(ierr);
    /* Default both faces to -humidity, so omitting these reproduces the old
     * single-value behaviour exactly. -humidity is parsed well above this. */
    user.rhovfix_lo = humidity;
    user.rhovfix_hi = humidity;
    ierr = PetscOptionsReal("-rhovfix_lo",
             "Vapor reservoir on the m=0 face (x=0), as a MULTIPLE of rho_vs(temp0). "
             "Defaults to -humidity. Specify at the 1e-6 level: a grain's own "
             "Gibbs-Thomson equilibrium is within ~5e-6 of rho_vs, so a "
             "humidity-style 0.99 swamps the curvature physics",
             "", user.rhovfix_lo, &user.rhovfix_lo, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-rhovfix_hi",
             "Vapor reservoir on the m=1 face (x=Lx), as a multiple of rho_vs(temp0)",
             "", user.rhovfix_hi, &user.rhovfix_hi, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-rhovfix_axis",
             "With -flag_BC_rhovfix: pin vapor only on this axis's two faces "
             "(0=x, 1=y, 2=z; -1 = every face, the legacy default). A pore "
             "channel wants 0 -- its top/bottom are solid wall, not reservoirs",
             "", user.rhovfix_axis, &user.rhovfix_axis, NULL); CHKERRQ(ierr);
    if (user.rhovfix_axis >= dim)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "-rhovfix_axis %d is out of range for dim=%d", (int)user.rhovfix_axis, (int)dim);
    ierr = PetscOptionsBool("-flag_Tdep",       "Temperature-dependent Gibbs-Thomson parameters",   "", user.flag_Tdep,  &user.flag_Tdep,  NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-dtCFL",           "Interface-CFL timestep limiter",                   "", user.flag_dtCFL, &user.flag_dtCFL, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-axisym",          "Axisymmetric r-z mode (x=axis, y=radius; grains on y=0)", "", user.axisym, &user.axisym, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-ic_grain_union",  "multi_grains IC from the union signed distance (eps-independent phi=0.5 contour) instead of summed tanh profiles", "", user.ic_grain_union, &user.ic_grain_union, NULL); CHKERRQ(ierr);
    /* Vapor diffusivity override: molecular D_v is the default; larger
     * values model convectively enhanced chamber transport (an effective
     * Sherwood-number correction) — see the 2026-07-12 Molaro validation
     * campaign, where the vapor-diffusion-limited neck rate fell ~3x below
     * experiment with every model-side mechanism eliminated. NOTE: the
     * kinetic derivation (tau_sub via M&F SI Eq. 9) uses this value too,
     * consistently. */
    ierr = PetscOptionsReal("-dif_vap",         "Vapor diffusivity in air [m^2/s]",                 "", user.dif_vap, &user.dif_vap, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-dtCFL_dphimax",   "Max pointwise |dphi| per step for the CFL limiter","", user.cfl_dphimax, &user.cfl_dphimax, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-decouple_phase_change", "Zero ice_t-driven source terms in R_tem/R_vap too (not just S_sub in R_ice)", "", user.decouple_phase_change, &user.decouple_phase_change, NULL); CHKERRQ(ierr);
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

    /* --- Prescribed contact angle at the regolith wall ------------------- */
    ierr = PetscOptionsReal("-gamma_ia", "Ice-air surface energy for Young's equation [J/m^2] (default: -Sigma_i)", "", gamma_ia, &gamma_ia, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-gamma_is", "Ice-regolith surface energy [J/m^2]", "", gamma_is, &gamma_is, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-gamma_as", "Air-regolith surface energy [J/m^2]", "", gamma_as, &gamma_as, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-contact_angle_deg", "DEBUG: set theta directly, bypassing Young's equation", "", contact_angle_deg, &contact_angle_deg, &contact_angle_set); CHKERRQ(ierr);
    ierr = PetscOptionsString("-wall_faces", "Domain faces that are regolith, e.g. \"y0,y1\" (default: none)", "", wall_faces, wall_faces, sizeof(wall_faces), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsInt("-stall_limit", "Abort if ||U|| is bit-identical for this many consecutive steps (0 disables)", "", user.stall_limit, &user.stall_limit, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-test_wall_measure", "DEBUG: assemble the wall term on a uniform phi=1/2 field, check its surface measure, and exit", "", test_wall_measure, &test_wall_measure, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-test_wall_jacobian", "DEBUG: check the analytic Jacobian against finite differences on the initial condition, and exit", "", test_wall_jacobian, &test_wall_jacobian, NULL); CHKERRQ(ierr);

    /* --- Output control -------------------------------------------------- */
    ierr = PetscOptionsInt("-outp", "Output control flag", "", user.outp, &user.outp, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsReal("-t_interv", "Output interval", "", user.t_interv, &user.t_interv, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-pf_output", "Enable output files", "", output, &output, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsBool("-pf_monitor", "Monitor the solution", "", monitor, &monitor, NULL); CHKERRQ(ierr);
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
    ierr = PetscOptionsReal("-t_start", "Clock value to resume from when restarting via -initial_cond [s]", "", t_start, &t_start, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-geom_file",
             "Load an igakit-generated IGA geometry (.dat) via IGARead, "
             "overriding -p/-C/-Nx/-Ny/-Nz axis setup with the geometry's own",
             "", geom_file, geom_file, sizeof(geom_file), NULL); CHKERRQ(ierr);
    ierr = PetscOptionsString("-ic_type",
             "Initial condition geometry (two_ice_grains_boundary|ice_slab|single_ice|multi_grains|multi_grains_file)",
             "src/<project>_main.c", ic_type, ic_type, sizeof(ic_type),
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

    /* --- Mesh/run temperature-consistency guard --------------------------
     * eps, and hence the mesh element count Nx = ceil(Lx*sqrt(2)/eps), is sized
     * by comp_eps.py through the TEMPERATURE-dependent kinetic bound. A geometry
     * .opts generated for temperature T sets -eps_valid_temp T. Running that
     * mesh at a different -temp means eps/mesh are inconsistent with the actual
     * kinetics — the diffuse interface is under- or over-resolved and the physics
     * is silently wrong. Fail loudly unless -eps_temp_override is set on purpose. */
    if (eps_valid_temp_set) {
        PetscReal dT = PetscAbsReal(temp - eps_valid_temp);
        if (dT > 1.0 && !eps_temp_override) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                    "\n*** TEMPERATURE MISMATCH ***\n"
                    "  run -temp          = %g C\n"
                    "  mesh -eps_valid_temp = %g C  (eps=%.4e sized here by comp_eps.py)\n"
                    "eps and the mesh element count are inconsistent with the kinetics at "
                    "%g C; the diffuse interface will be under/over-resolved. Regenerate the "
                    "geometry with comp_eps.py for %g C (see the builder's printed command), "
                    "or pass -eps_temp_override 1 to proceed intentionally.",
                    (double)temp, (double)eps_valid_temp, (double)eps,
                    (double)temp, (double)temp);
        } else if (dT > 1.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                    "  *** WARNING: -temp %g C != mesh -eps_valid_temp %g C; "
                    "-eps_temp_override active. eps/mesh may under/over-resolve the "
                    "interface. ***\n", (double)temp, (double)eps_valid_temp);
        }
    }

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

    /* Cap output volume. With -outp N, every N-th accepted step writes a
     * snapshot, and t_final/dtmax is a LOWER bound on the total step count
     * (the adaptive ramp only adds steps). If even that bound exceeds 1000
     * files, per-step output is a disk hazard (each sol_*.dat is ~0.5 MB and
     * the vtk conversion doubles it): switch to time-uniform output with
     * exactly 1000 snapshots. Runs short enough to stay under 1000 files
     * keep the -outp per-step behavior (naturally log-spaced under the
     * adaptive dt). */
    if (user.outp >= 1 && dtmax > 0.0 &&
        t_final / (dtmax * (PetscReal)user.outp) > 1000.0) {
        user.outp     = 0;
        n_out         = 1000;
        user.t_interv = t_final / (PetscReal)(n_out - 1);
        PetscPrintf(PETSC_COMM_WORLD,
                    "Output cap: t_final/dtmax = %.0f exceeds 1000 snapshots; "
                    "switching to %d time-uniform outputs (t_interv = %.3e s)\n",
                    (double)(t_final / dtmax), (int)n_out,
                    (double)user.t_interv);
    }

    /* Gibbs-Thomson kinetic parameters */
    user.diff_sub = 0.5 * (user.thcond_air / user.rho_air / user.cp_air + user.thcond_ice / user.rho_ice / user.cp_ice);

    user.Etai = Sigma_i; /* Ice surface energy in the double-well free energy */
    user.Etaa = Sigma_a; /* Air surface energy in the double-well free energy */

    /* ---- Prescribed contact angle: resolve gamma's -> cos(theta) ---------
     * gamma_ia defaults to (Sigma_i + Sigma_a)/2, NOT to Sigma_i.
     *
     * Sigma_i is not the ice-air surface energy except in the degenerate case
     * Sigma_a = Sigma_i. In the phase-field convention these projects inherit,
     *     Sigma_i = gamma_ia + gamma_is - gamma_as
     *     Sigma_a = gamma_ia + gamma_as - gamma_is
     * so Sigma_i + Sigma_a = 2*gamma_ia and the ice-air energy is their MEAN.
     * This solver's own double-well coefficient says the same thing a few lines
     * up: C = (Sigma_i + Sigma_a)/2 + Lambda. With the lunar defaults
     * (0.109, 0.132) that is 0.1205 J/m^2; taking Sigma_i alone understates it
     * by 10.6 %, which propagates straight into theta through
     * cos(theta) = (gamma_as - gamma_is)/gamma_ia.
     *
     * It happened to be harmless in enceladus_DSM, which sets Sigma_a = Sigma_i
     * deliberately, and in every run so far, which passed -gamma_ia explicitly.
     * Pass -gamma_ia when the substrate matters; the default is only a
     * fallback. */
    if (gamma_ia <= 0.0) gamma_ia = 0.5 * (Sigma_i + Sigma_a);
    user.gamma_ia = gamma_ia;
    user.gamma_is = gamma_is;
    user.gamma_as = gamma_as;
    user.costhet_direct = contact_angle_set;

    if (contact_angle_set) {
        /* Debug override: theta straight from the CLI, Young bypassed. */
        if (contact_angle_deg < 0.0 || contact_angle_deg > 180.0)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                    "-contact_angle_deg %g is outside [0, 180]",
                    (double)contact_angle_deg);
        user.costhet = PetscCosReal(contact_angle_deg * PETSC_PI / 180.0);
    } else {
        const PetscReal dgam = gamma_as - gamma_is;
        if (PetscAbsReal(dgam) > gamma_ia)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                    "Young's equation has no solution: |gamma_as - gamma_is| = "
                    "|%g - %g| = %g exceeds gamma_ia = %g. A contact angle "
                    "exists only for |gamma_as - gamma_is| <= gamma_ia; outside "
                    "that range one phase wets the regolith completely.",
                    (double)gamma_as, (double)gamma_is,
                    (double)PetscAbsReal(dgam), (double)gamma_ia);
        user.costhet = dgam / gamma_ia;
    }

    /* ---- Which faces are regolith ---------------------------------------
     * -wall_faces "y0,y1": axis letter x/y/z followed by side 0 (low
     * coordinate) or 1 (high). Faces left out keep the natural Neumann
     * dphi/dn = 0 they have always had, so omitting the flag is a no-op. */
    ierr = PetscMemzero(user.wall_face, sizeof(user.wall_face)); CHKERRQ(ierr);
    user.wall_any = PETSC_FALSE;
    if (wall_faces[0] != '\0') {
        for (const char *c = wall_faces; *c != '\0'; c++) {
            PetscInt axis, side;
            if (*c == ',' || *c == ' ') continue;
            if      (*c == 'x' || *c == 'X') axis = 0;
            else if (*c == 'y' || *c == 'Y') axis = 1;
            else if (*c == 'z' || *c == 'Z') axis = 2;
            else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                         "-wall_faces \"%s\": expected an axis letter x/y/z, got '%c'. "
                         "Use a comma-separated list like \"y0,y1\".", wall_faces, *c);
            c++;
            if      (*c == '0') side = 0;
            else if (*c == '1') side = 1;
            else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                         "-wall_faces \"%s\": axis letter must be followed by side "
                         "0 (low) or 1 (high).", wall_faces);
            if (axis >= dim)
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                        "-wall_faces \"%s\" names axis %d but -dim is %d.",
                        wall_faces, (int)axis, (int)dim);
            user.wall_face[axis][side] = PETSC_TRUE;
            user.wall_any = PETSC_TRUE;
        }
    }
    if (user.wall_any && user.periodic == 1)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                "-wall_faces is set but -periodic 1: a periodic domain has no "
                "walls to apply a contact angle to. Use -periodic 0.");

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
     * not distort the equilibrium diffuse profile the way a forcing term
     * decoupled from mob_sub's restoring strength would. */
    {
        PetscReal  beta_sub0_cli = -1.0;
        PetscBool  set_beta      = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-beta_sub0",
                                   &beta_sub0_cli, &set_beta); CHKERRQ(ierr);
        if (set_beta && beta_sub0_cli > 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -beta_sub0 override: %.4e -> %.4e (factor %.2f)\n",
                        user.beta_sub0, beta_sub0_cli, user.beta_sub0 / beta_sub0_cli);
            user.beta_sub0 = beta_sub0_cli;
        }
    }
    {
        PetscReal  d0_sub0_cli = -1.0;
        PetscBool  set_d0      = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-d0_sub0",
                                   &d0_sub0_cli, &set_d0); CHKERRQ(ierr);
        if (set_d0 && d0_sub0_cli > 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -d0_sub0 override: %.4e -> %.4e (factor %.2f)\n",
                        user.d0_sub0, d0_sub0_cli, d0_sub0_cli / user.d0_sub0);
            user.d0_sub0 = d0_sub0_cli;
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
    /* Thin-interface counter-terms, ON by default (-thin_iface_corr).
     *
     *   tau_sub = eps^2*beta/d0  +  a1*a2*(eps^3/d0)*(1/D_therm + 1/D_v)
     *
     * DEFAULT CHANGED 2026-09-13, OFF -> ON. The argument for OFF (below, and
     * docs/gt_deficit/) bounds only the counter-terms PROPORTIONAL TO v_n, and
     * shows those vanish here because the one-sided D_v*phi_a makes the inner
     * deviation of sigma identically null. It does not bound the spurious
     * SURFACE DIFFUSION that a one-sided model produces without an
     * anti-trapping current -- gt_deficit.tex says so explicitly, and notes
     * that operator "is not observable in this measurement". This model has no
     * anti-trapping current, so that term is live and uncontrolled, and the
     * beta-agreement evidence that motivated OFF simply cannot see it.
     *
     * Cost of being wrong in this direction is bounded and known: beta is
     * realised 1.03-1.12x larger than requested over the eps range in use
     * (correction/tau = 2.9% at eps = 0.75 um, 5.9% at 1.50, 11.7% at 3.00 --
     * it scales as eps). Cost of being wrong the other way is an uncontrolled
     * spurious operator on every interface.
     *
     * CAVEAT, unresolved: the thermal counter-term is the HISTORICAL form.
     * a2*eps/diff_sub treats latent heat as driving sigma at unit strength,
     * whereas temperature reaches sigma only through rho_vs(T) and so carries a
     * Clausius-Clapeyron factor; gt_deficit.tex 198-200 gives the correct
     * coefficient as (d rho_vs/dT)*L_sub/k, smaller by 11x to 1300x. So
     * tau_therm here is over-weighted. At eps = 0.75 um that is 45.7 s of a
     * 2261 s tau_sub, and correcting it would move tau_sub by ~1.8%; the vapor
     * term (18.7 s) is unaffected. Worth fixing on its own terms.
     *
     * Historical note on the OFF rationale:
     *
     * These inflate tau_sub by the spurious O(eps) kinetic contribution that the
     * sharp-interface asymptotics are expected to subtract back off, so that the
     * realised kinetic coefficient equals the requested -beta_sub0. For this
     * model that contribution is zero: the vapour diffusivity is one-sided
     * (D_v*phi_a), which makes the inner deviation of sigma identically null, so
     * nothing is subtracted and the inflation survives as an error in beta --
     * 1.22x at the -20 C wedge parameters, 1.49x at the Molaro ones. Measured
     * over a 40x sweep in -beta_sub0, beta_fit/beta_bare = 1.0009 +/- 0.0007.
     *
     * -thin_iface_corr 1 restores the terms in their historical form, for
     * comparison against earlier runs. (In lunar this is not bit-for-bit: D_v is
     * now taken at temp0 rather than at 0 C, which moves tau_sub by 0.7%.) It is
     * not a corrected form: the
     * thermal term additionally omits the Clausius-Clapeyron factor by which
     * temperature acts on sigma, and is over-weighted by 11-1300x as a result.
     * Since the physically correct correction for this model is ~0, a "fixed"
     * ON branch would be indistinguishable from OFF and is not provided.
     *
     * See docs/gt_deficit/. */
    /* D_v is taken at temp0, matching enceladus and the pointwise path; the base
     * 0 C value used here previously differed by 13% in that term at -20 C. */
    {
        PetscScalar dv_T0;
        VaporDiffus(&user, (PetscScalar)user.temp0, &dv_T0, NULL);
        PetscReal c_vap = 0.0, c_therm = 0.0;
        if (user.thin_iface_corr) {
            c_therm = a2 * user.eps / user.diff_sub;
            c_vap   = a2 * user.eps / PetscRealPart(dv_T0);
        }
        user.tau_kin   = user.eps * lambda_sub * (beta_sub / a1);
        user.tau_therm = user.eps * lambda_sub * c_therm;
        user.tau_vap   = user.eps * lambda_sub * c_vap;
        tau_sub = user.tau_kin + user.tau_therm + user.tau_vap;
    }
    user.mob_sub = 1 * user.eps / 3.0 / tau_sub; /* Mobility parameter for sublimation */
    user.alph_sub = lambda_sub / tau_sub;  /* Phase change rate parameter, eq.(9) Moure & Fu (2024) SI */

    /* Allow per-test override of mob_sub via -mob_sub <value>. Tests with
     * very stiff geometries (touching/merging grains in 2D) can reduce
     * mob_sub by ~10x to trade kinetics speed for AC stability. The
     * physical value computed above is the default. */
    {
        PetscReal  mob_sub_cli = -1.0;
        PetscBool  set_mob     = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-mob_sub",
                                   &mob_sub_cli, &set_mob); CHKERRQ(ierr);
        if (set_mob && mob_sub_cli > 0.0) {
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
        PetscBool  set_alph     = PETSC_FALSE;
        ierr = PetscOptionsGetReal(NULL, NULL, "-alph_sub",
                                   &alph_sub_cli, &set_alph); CHKERRQ(ierr);
        if (set_alph && alph_sub_cli >= 0.0) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  -alph_sub override: %.2e -> %.2e%s\n",
                        user.alph_sub, alph_sub_cli,
                        (alph_sub_cli == 0.0) ? " (phase-change coupling DISABLED)" : "");
            user.alph_sub = alph_sub_cli;
        }
    }

    /* CLI overrides for xi_T and xi_v. */
    ierr = PetscOptionsGetReal(NULL, NULL, "-xi_T", &user.xi_T, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-xi_v", &user.xi_v, NULL); CHKERRQ(ierr);

    /* Create IGA and set up problem */
    IGA iga;
    ierr = IGACreate(PETSC_COMM_WORLD, &iga); CHKERRQ(ierr);
    ierr = IGASetDim(iga, dim); CHKERRQ(ierr);
    /* Guard: the residual/Jacobian in assembly.c and the Field struct in
     * NASA_types.h assume exactly this many DOF per node (currently 3: ice,
     * temperature, vapor). -dof is exposed for future multi-field work (e.g.
     * an explicit sediment phase), but changing it without also updating Field
     * and assembly.c silently misinterprets the solution vector. Fail loudly. */
    {
        const PetscInt dof_fields = (PetscInt)(sizeof(Field) / sizeof(PetscScalar));
        if (dof != dof_fields)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                    "-dof %d does not match the %d-field solver (Field struct in "
                    "NASA_types.h and the [%d] arrays in assembly.c must change "
                    "together).", (int)dof, (int)dof_fields, (int)dof_fields);
    }
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

    /* Regolith faces get a boundary form so PetIGA visits their quadrature
     * points: the prescribed-contact-angle wall term is assembled there (see
     * the pnt->atboundary branches in assembly.c). Faces not listed are never
     * visited and keep the natural Neumann dphi/dn = 0. Enabling the form with
     * cos(theta) = 0 is still an exact no-op -- the branch returns immediately
     * -- so a theta = 90 deg run reproduces the old behaviour bit for bit. */
    for (PetscInt l = 0; l < dim; l++) {
        for (PetscInt m = 0; m < 2; m++) {
            if (!user.wall_face[l][m]) continue;
            ierr = IGASetBoundaryForm(iga, l, m, PETSC_TRUE); CHKERRQ(ierr);
        }
    }
    // ierr = IGASetFormIJacobian(iga, IGAFormIJacobianFD, &user); CHKERRQ(ierr);

    /* Boundary conditions (could 'functionalize' this at some point) */

    /* Record of what is actually pinned on each face, so the per-wall report
     * further down prints the real configuration rather than re-deriving it
     * from the flags. The carve-outs below (-rhovfix_axis, the gradient-
     * parallel skip, the axisymmetric y=0 exemption) make "which faces ended
     * up Dirichlet" genuinely awkward to restate, and a report that restates
     * it is a second implementation waiting to disagree with this one.
     * Indexed [axis][side][dof]; side 0 is the low-coordinate face. */
    PetscBool bc_dirichlet[3][2][3];
    PetscReal bc_value[3][2][3];
    ierr = PetscMemzero(bc_dirichlet, sizeof(bc_dirichlet)); CHKERRQ(ierr);
    ierr = PetscMemzero(bc_value,     sizeof(bc_value));     CHKERRQ(ierr);

    // Set vapor density BCs
    if (flag_BC_rhovfix) {
        PetscReal rho0_vs;
        RhoVS_I(&user, user.temp0, &rho0_vs, NULL);
        for (PetscInt l = 0; l < dim; l++) {
            /* -rhovfix_axis restricts the reservoir to ONE axis's two faces.
             * In a pore-channel domain only the open ENDS are reservoirs; the
             * top/bottom are solid regolith wall. Pinning vapor on the wall
             * faces would feed the ice directly at its contact line, driving
             * growth along the whole wall and corrupting the contact angle —
             * so a wall-bounded channel wants -rhovfix_axis 0 (x-faces only),
             * not the legacy all-faces behaviour (-1). */
            if (user.rhovfix_axis >= 0 && l != user.rhovfix_axis) continue;
            for (PetscInt m = 0; m < 2; m++) {
                /* Axisymmetric mode: the y = 0 face is the symmetry AXIS —
                 * interior space of the 3D problem, usually with ice sitting
                 * on it — never a reservoir. Keep it natural Neumann (the
                 * exact axis condition) and pin vapor only on the true
                 * outer boundaries. */
                if (user.axisym && l == 1 && m == 0) continue;
                /* m=0 is the low-coordinate face, m=1 the high one. */
                PetscReal frac = (m == 0) ? user.rhovfix_lo : user.rhovfix_hi;
                ierr = IGASetBoundaryValue(iga, l, m, 2, frac * rho0_vs); CHKERRQ(ierr);
                bc_dirichlet[l][m][2] = PETSC_TRUE;
                bc_value[l][m][2]     = frac * rho0_vs;
            }
        }
    }

    // Set temperature BCs
    if (flag_BC_Tfix) {
        PetscReal T_BC[3][2] = {{0}};
        PetscReal LL[3]      = {user.Lx, user.Ly, user.Lz};
        /* A temperature gradient is sustained by pinning T (Dirichlet) on the
         * two faces PERPENDICULAR to it and leaving the faces PARALLEL to it
         * insulating (natural Neumann). IGASetBoundaryValue writes a single
         * constant per face, which is exact for the perpendicular faces (T is
         * constant along them) but WRONG for a parallel face: a purely
         * transverse gradient (e.g. dT/dx with the grains on the y=0 edge)
         * has T varying along the top/bottom edges, so pinning them to a
         * uniform temp0 would flatten the gradient right where the grains sit
         * and erase the driving ΔT. So skip any face whose axis carries no
         * gradient component -- UNLESS there is no gradient at all, in which
         * case every face is pinned to temp0 (a uniform isothermal reservoir,
         * the original behaviour). */
        PetscBool has_grad = (PetscBool)(user.grad_temp0[0] != 0.0 ||
                                         user.grad_temp0[1] != 0.0 ||
                                         user.grad_temp0[2] != 0.0);
        for (PetscInt l = 0; l < dim; l++) {
            /* Face parallel to the gradient (no component along axis l): leave
             * it insulating so the transverse gradient survives on that edge. */
            if (has_grad && user.grad_temp0[l] == 0.0) continue;
            for (PetscInt m = 0; m < 2; m++) {
                /* Axisym: skip the y=0 axis face (interior space, not a
                 * thermal reservoir) — same guard as the vapor BC above. */
                if (user.axisym && l == 1 && m == 0) continue;
                T_BC[l][m] = user.temp0 + (2.0 * m - 1.0) * user.grad_temp0[l] * LL[l] / 2.0;
                ierr = IGASetBoundaryValue(iga, l, m, 1, T_BC[l][m]); CHKERRQ(ierr);
                bc_dirichlet[l][m][1] = PETSC_TRUE;
                bc_value[l][m][1]     = T_BC[l][m];
            }
        }
    }

    /* Set up TS */
    TS ts;
    ierr = IGACreateTS(iga, &ts); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, t_final); CHKERRQ(ierr);
    if (t_start > 0.0) { ierr = TSSetTime(ts, t_start); CHKERRQ(ierr); }
    ierr = TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts, delt_t); CHKERRQ(ierr);
    ierr = TSSetType(ts, TSALPHA); CHKERRQ(ierr);
    ierr = TSAlphaSetRadius(ts, 0.5); CHKERRQ(ierr);
    if (monitor) { ierr = TSMonitorSet(ts, Monitor, &user, NULL); CHKERRQ(ierr); }
    if (output) { ierr = TSMonitorSet(ts, OutputMonitor, &user, NULL); CHKERRQ(ierr); }
    /* Interface-CFL limiter: registered unconditionally (cheap — one vector
     * diff + stride norm per accepted step); disable with -dtCFL 0. */
    ierr = TSMonitorSet(ts, InterfaceCFLMonitor, &user, NULL); CHKERRQ(ierr);

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

    /* Re-route the SNES residual through the domain-error-syncing wrapper
     * (see SNESTSFormFunction_DomainErrSync above). TSGetSNES() has already
     * installed SNESTSFormFunction with ctx=ts, and TSSetUp() only installs
     * it when unset, so overriding here sticks. */
    ierr = SNESSetFunction(nonlin, NULL, SNESTSFormFunction_DomainErrSync, ts); CHKERRQ(ierr);

    /* Bound-constrained Newton solve: enforce 0 <= ice <= 1 directly
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
     * Tightened to strict [0,1] (2026-06-22, same day, third trial): testing
     * whether the smaller -dtmax (1.0e4) alone is now enough to keep this
     * physically-correct bound well-behaved, without any VI slack at all.
     *
     * Loosened slightly to -0.01/1.01 (2026-06-22, same day, fourth trial):
     * a small amount of slack between strict [0,1] and the earlier
     * -0.05/1.05, to see where the smaller-dtmax regime actually needs the
     * line drawn.
     *
     * Technically incorrect (allows larger unphysical excursions outside
     * [0,1]); intentional tradeoff to get a visibly-evolving result for the
     * 2026-06-23 conference. Revisit the real fix (per-DOF-block atol
     * matched to each field's natural scale) after the conference.
     *
     * CLI-togglable for direct A/B comparison against the pre-VI solver
     * (2026-06-23): -vi_bounds 0 skips this whole block, so -snes_type can
     * be set back to plain newtonls on the command line to reproduce the
     * original unbounded Newton solve exactly (no VI machinery at all, not
     * just unenforced bounds). -vi_lo/-vi_hi override the bound values
     * themselves (default -0.01/1.01, the current production setting) so
     * the same binary can also run with strict [0,1] bounds for comparison
     * without a rebuild. */
    PetscBool vi_bounds = PETSC_TRUE;
    PetscReal vi_lo = -0.01, vi_hi = 1.01;
    ierr = PetscOptionsGetBool(NULL, NULL, "-vi_bounds", &vi_bounds, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-vi_lo", &vi_lo, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(NULL, NULL, "-vi_hi", &vi_hi, NULL); CHKERRQ(ierr);
    if (vi_bounds) {
        Vec Xl, Xu;
        ierr = IGACreateVec(iga, &Xl); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Xu); CHKERRQ(ierr);
        ierr = VecStrideSet(Xl, 0, vi_lo);           CHKERRQ(ierr);
        ierr = VecStrideSet(Xu, 0, vi_hi);           CHKERRQ(ierr);
        ierr = VecStrideSet(Xl, 1, PETSC_NINFINITY); CHKERRQ(ierr);
        ierr = VecStrideSet(Xu, 1, PETSC_INFINITY);  CHKERRQ(ierr);
        ierr = VecStrideSet(Xl, 2, PETSC_NINFINITY); CHKERRQ(ierr);
        ierr = VecStrideSet(Xu, 2, PETSC_INFINITY);  CHKERRQ(ierr);
        ierr = SNESVISetVariableBounds(nonlin, Xl, Xu); CHKERRQ(ierr);
        ierr = VecDestroy(&Xl); CHKERRQ(ierr);
        ierr = VecDestroy(&Xu); CHKERRQ(ierr);
    }

    /* ========================================================================
     * Comprehensive parameter summary — printed here so all true values are
     * available: Nx/Ny/p/C (after IGASetUp + geom_file override), derived
     * kinetic params (tau_sub/lambda_sub/mob_sub/alph_sub), and VI bounds
     * (vi_lo/vi_hi). Override notifications above document any CLI changes.
     * ======================================================================== */
    PetscPrintf(PETSC_COMM_WORLD,
        "\n================================================================================\n"
        " PHASE-FIELD SIMULATION PARAMETERS\n"
        "================================================================================\n");

    /* --- Mesh & discretization -------------------------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n MESH & DISCRETIZATION\n");
    if (user.axisym)
        PetscPrintf(PETSC_COMM_WORLD,
                    "   AXISYMMETRIC r-z mode: x = axis (z), y = radius (r); "
                    "integrands weighted by 2*pi*r\n");
    if (dim == 1) {
        PetscPrintf(PETSC_COMM_WORLD, "   Nx = %d%s\n",
                    Nx, geom_file[0] ? "  [from -geom_file]" : "");
        PetscPrintf(PETSC_COMM_WORLD, "   Lx = %.4e m\n", Lx);
        PetscPrintf(PETSC_COMM_WORLD, "   dx = %.4e m\n", Lx / Nx);
    } else if (dim == 2) {
        PetscPrintf(PETSC_COMM_WORLD, "   Nx = %d,  Ny = %d%s\n",
                    Nx, Ny, geom_file[0] ? "  [from -geom_file]" : "");
        PetscPrintf(PETSC_COMM_WORLD, "   Lx = %.4e m,  Ly = %.4e m\n", Lx, Ly);
        PetscPrintf(PETSC_COMM_WORLD, "   dx = %.4e m,  dy = %.4e m\n", Lx / Nx, Ly / Ny);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "   Nx = %d,  Ny = %d,  Nz = %d%s\n",
                    Nx, Ny, Nz, geom_file[0] ? "  [from -geom_file]" : "");
        PetscPrintf(PETSC_COMM_WORLD, "   Lx = %.4e m,  Ly = %.4e m,  Lz = %.4e m\n", Lx, Ly, Lz);
        PetscPrintf(PETSC_COMM_WORLD, "   dx = %.4e m,  dy = %.4e m,  dz = %.4e m\n",
                    Lx / Nx, Ly / Ny, Lz / Nz);
    }
    {
        const char *pname = (p == 1) ? "linear" : (p == 2) ? "quadratic" : (p == 3) ? "cubic" : "order";
        PetscPrintf(PETSC_COMM_WORLD, "   p = %d (%s),  C = %d\n", p, pname, C);
    }

    /* --- Phase-field interface -------------------------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n PHASE-FIELD INTERFACE\n");
    PetscPrintf(PETSC_COMM_WORLD, "   eps      =  %.4e m\n", user.eps);
    PetscPrintf(PETSC_COMM_WORLD, "   Sigma_i  =  %.4e J/m²   (ice surface energy)\n", user.Etai);
    PetscPrintf(PETSC_COMM_WORLD, "   Sigma_a  =  %.4e J/m²   (air surface energy)\n", user.Etaa);
    PetscPrintf(PETSC_COMM_WORLD, "   Lambda   =  %.4e\n", user.Lambd);

    /* --- Environment & initial conditions --------------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n ENVIRONMENT & INITIAL CONDITIONS\n");
    PetscPrintf(PETSC_COMM_WORLD, "   T0       = %7.2f °C,  humidity = %.4f\n", temp, humidity);
    if (dim == 1)
        PetscPrintf(PETSC_COMM_WORLD, "   grad_T   =  (%.4e) °C/m\n", grad_temp0[0]);
    else if (dim == 2)
        PetscPrintf(PETSC_COMM_WORLD, "   grad_T   =  (%.4e, %.4e) °C/m\n",
                    grad_temp0[0], grad_temp0[1]);
    else
        PetscPrintf(PETSC_COMM_WORLD, "   grad_T   =  (%.4e, %.4e, %.4e) °C/m\n",
                    grad_temp0[0], grad_temp0[1], grad_temp0[2]);

    /* --- Time stepping ----------------------------------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n TIME STEPPING\n");
    PetscPrintf(PETSC_COMM_WORLD, "   dt0      =  %.4e s\n", delt_t);
    PetscPrintf(PETSC_COMM_WORLD, "   t_final  =  %.4e s,  n_out = %d\n", t_final, n_out);
    if (adap == 1) {
        PetscPrintf(PETSC_COMM_WORLD,
                    "   adaptive    ON  (NRmin = %d,  NRmax = %d,  factor = %.4f)\n",
                    NRmin, NRmax, factor);
        PetscPrintf(PETSC_COMM_WORLD, "   dtmin    =  %.4e s,  dtmax = %.4e s\n", dtmin, dtmax);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "   adaptive    OFF (fixed dt)\n");
    }

    /* --- Transport & thermophysical properties ----------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n TRANSPORT & THERMOPHYSICAL PROPERTIES\n");
    PetscPrintf(PETSC_COMM_WORLD, "   D_v      =  %.4e m²/s    (vapor diffusivity)\n", user.dif_vap);
    PetscPrintf(PETSC_COMM_WORLD, "   k_ice    =  %.4e W/m/K,  k_air  = %.4e W/m/K\n",
                user.thcond_ice, user.thcond_air);
    PetscPrintf(PETSC_COMM_WORLD, "   rho_ice  =  %.4e kg/m³,  rho_air = %.4e kg/m³\n",
                user.rho_ice, user.rho_air);
    PetscPrintf(PETSC_COMM_WORLD, "   cp_ice   =  %.4e J/kg/K, cp_air = %.4e J/kg/K\n",
                user.cp_ice, user.cp_air);
    PetscPrintf(PETSC_COMM_WORLD, "   lat_sub  =  %.4e J/kg    (latent heat of sublimation)\n",
                user.lat_sub);

    /* --- Phase-change kinetics -------------------------------------------- */
    PetscPrintf(PETSC_COMM_WORLD, "\n PHASE-CHANGE KINETICS\n");
    PetscPrintf(PETSC_COMM_WORLD, "   rho_ice/rho_vs   = %.4e   (density ratio at T0)\n", rho_rhovs);
    PetscPrintf(PETSC_COMM_WORLD, "   d0_sub0          = %.4e m   (capillary length; M&F use 1e-7)\n",
                user.d0_sub0);
    PetscPrintf(PETSC_COMM_WORLD, "   beta_sub (K&P β₀, M&F β_sub, UNSCALED) = %.4e s/m"
                "   [M&F range: 2e4–2e6]\n", user.beta_sub0);
    PetscPrintf(PETSC_COMM_WORLD, "   beta_sub (SCALED = β₀·ρ_vs/ρ_ice = β_HK) = %.4e s/m\n",
                beta_sub);
    if (!user.flag_Tdep) {
        PetscPrintf(PETSC_COMM_WORLD, "   lambda   =  %.4e\n", lambda_sub);
        PetscPrintf(PETSC_COMM_WORLD, "   tau_sub  =  %.4e s\n", tau_sub);
        /* dtmax must scale with tau_sub, which goes as eps^2. A dtmax that is
         * safe at one resolution is not safe at a finer one, and the failure is
         * silent: batch 2026-09-12's eps = 0.75 um run inherited a dtmax sized
         * for eps = 1.50 um from a SHARED experiment file, reached
         * dtmax/tau_sub = 0.91 -- one step spanning the whole interface
         * relaxation time -- drove the phase field out of [phase_lo, phase_hi]
         * at step 383, and then froze. The clock ran on to t_final at dtmax and
         * the run reported success; only counting distinct sol_*.dat
         * fingerprints revealed it. Warn loudly instead.
         *
         * Threshold from measurement, not from that batch: dtmax_study.sh swept
         * dtmax/tau_sub from 0.051 to 0.815 and theta_inf moved 0.0028 deg,
         * with phi never leaving [0,1] and the CFL limiter never firing. The
         * 0.91 stall that originally motivated a tight threshold was the
         * wall-term sign flip, not dt; it is fixed by the clamp. So warn only
         * above 1.0, where a single step would exceed the entire interface
         * relaxation time. */
        if (dtmax > 0.0 && tau_sub > 0.0) {
            const PetscReal ratio = dtmax / tau_sub;
            PetscPrintf(PETSC_COMM_WORLD,
                "   dtmax/tau_sub = %.4f%s\n", (double)ratio,
                (ratio > 1.0) ? "   <-- SEE WARNING BELOW" : "   (ok, <= 1.0)");
            if (ratio > 1.0)
                PetscPrintf(PETSC_COMM_WORLD,
                    "\n   *** WARNING: dtmax = %.3e s is %.2f x tau_sub = %.3e s.\n"
                    "       A step longer than the whole interface relaxation time cannot\n"
                    "       and the solve can stall while the clock keeps advancing --\n"
                    "       a run that LOOKS successful but whose solution stopped\n"
                    "       changing. tau_sub scales as eps^2, so set -dtmax in the\n"
                    "       GEOMETRY file (which owns eps), not in a shared experiment\n"
                    "       file. Suggested: -dtmax %.2e  (tau_sub/10)\n\n",
                    (double)dtmax, (double)ratio, (double)tau_sub,
                    (double)(tau_sub / 10.0));
        }
        {   /* Realised vs requested kinetics. beta_bare = tau_sub*d0_sub0/eps^2
             * is the coefficient the sharp-interface limit actually delivers; it
             * equals -beta_sub0 only when the thin-interface counter-terms are
             * absent, which is the default. See docs/gt_deficit/. */
            PetscReal b_bare = tau_sub * user.d0_sub0 / (user.eps * user.eps);
            PetscPrintf(PETSC_COMM_WORLD,
                "   tau_sub terms:  kinetic %.4e s", user.tau_kin);
            if (user.thin_iface_corr)
                PetscPrintf(PETSC_COMM_WORLD, " + thermal %.4e + vapor %.4e",
                            user.tau_therm, user.tau_vap);
            PetscPrintf(PETSC_COMM_WORLD, "   (-thin_iface_corr %d)\n",
                        (int)user.thin_iface_corr);
            PetscPrintf(PETSC_COMM_WORLD,
                "   beta realised  =  %.4e s/m  =  %.4f x -beta_sub0%s\n",
                b_bare, b_bare / user.beta_sub0,
                user.thin_iface_corr ? "   <-- NOT the beta you requested" : "");
        }

        PetscPrintf(PETSC_COMM_WORLD, "   mob_sub  =  %.4e m/s   [M&F: 4.33e-7]\n", user.mob_sub);
        PetscPrintf(PETSC_COMM_WORLD, "   alph_sub =  %.4e 1/s%s\n", user.alph_sub,
                    (user.alph_sub == 0.0) ? "   (phase-change DECOUPLED)" : "");
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "   [temperature-dependent kinetics active]\n");
    }
    if (user.decouple_phase_change)
        PetscPrintf(PETSC_COMM_WORLD,
                    "   decouple_phase_change: ON  (pure AC dynamics, no latent-heat/mass source)\n");
    PetscPrintf(PETSC_COMM_WORLD, "   xi_T     =  %.4e   (thermal conduction + latent heat scale)\n",
                user.xi_T);
    PetscPrintf(PETSC_COMM_WORLD, "   xi_v     =  %.4e   (vapor diffusion + rho_ice source scale)\n",
                user.xi_v);

    /* --- Solver ------------------------------------------------------------ */
    PetscPrintf(PETSC_COMM_WORLD, "\n SOLVER\n");
    if (vi_bounds)
        PetscPrintf(PETSC_COMM_WORLD, "   VI bounds:  ON  (ice in [%.4f, %.4f])\n", vi_lo, vi_hi);
    else
        PetscPrintf(PETSC_COMM_WORLD,
                    "   VI bounds:  OFF (unbounded Newton — pair with -snes_type newtonls)\n");

    /* --- Boundary conditions ------------------------------------------------
     * One row per wall, one column per equation, printed from the bc_dirichlet
     * record filled in where the conditions were actually applied. The old
     * three-line summary reported only the FLAGS, which is not the same thing:
     * -flag_BC_Tfix with a transverse gradient pins two walls and leaves the
     * other two insulating, and -rhovfix_axis pins only one axis's pair, but
     * both cases printed a flat "T: Dirichlet" / "rho_v: Dirichlet" that
     * implied all four. */
    {
        static const char *const wall_name[3][2] = {
            {"x = 0   (left)",   "x = Lx  (right)"},
            {"y = 0   (bottom)", "y = Ly  (top)"},
            {"z = 0   (back)",   "z = Lz  (front)"},
        };
        const char *unit[3] = {"", "C", "kg/m3"};
        /* rho_v runs ~1e-4 kg/m3, which %g renders as an unreadable string of
         * leading zeros; T is a plain few-digit number. Format each to suit. */
        const char *vfmt[3] = {"", "Dirichlet %.4f %s", "Dirichlet %.4e %s"};

        PetscPrintf(PETSC_COMM_WORLD, "\n BOUNDARY CONDITIONS\n");
        if (user.periodic == 1) {
            PetscPrintf(PETSC_COMM_WORLD,
                "   PERIODIC on every axis — the domain has no walls, and all\n"
                "   Dirichlet flags are forced off (-flag_BC_Tfix/-flag_BC_rhovfix\n"
                "   are ignored under -periodic 1).\n");
        } else {
            PetscPrintf(PETSC_COMM_WORLD,
                "   Per wall, per equation. \"Neumann\" is the NATURAL condition —\n"
                "   zero normal flux, imposed by omission rather than explicitly:\n"
                "   dphi/dn = 0 for the phase field (a 90° contact angle where ice\n"
                "   meets the wall), zero heat flux for T (insulating), and zero\n"
                "   vapor flux for rho_v (sealed — no mass enters or leaves there).\n\n");
            PetscPrintf(PETSC_COMM_WORLD,
                "   %-18s %-22s %-24s %s\n",
                "Wall", "phi_i (ice)", "T (temperature)", "rho_v (vapor)");
            PetscPrintf(PETSC_COMM_WORLD,
                "   %-18s %-22s %-24s %s\n",
                "------------------", "----------------------",
                "------------------------", "------------------------");
            for (PetscInt l = 0; l < dim; l++) {
                for (PetscInt m = 0; m < 2; m++) {
                    char cell[3][40];
                    const char *name = wall_name[l][m];
                    /* In axisymmetric r-z mode the y=0 face is the symmetry
                     * axis, not a physical wall — say so, since "Neumann"
                     * there is the exact axis condition, not a modelling
                     * choice about a boundary. */
                    if (user.axisym && l == 1 && m == 0) name = "r = 0   (axis)";
                    if (user.wall_face[l][m] && user.costhet != 0.0)
                        PetscSNPrintf(cell[0], sizeof(cell[0]),
                                      "regolith  theta=%.1f°",
                                      (double)(PetscAcosReal(user.costhet)
                                               * 180.0 / PETSC_PI));
                    else
                        PetscSNPrintf(cell[0], sizeof(cell[0]), "Neumann  dphi/dn=0");
                    for (PetscInt d = 1; d < 3; d++) {
                        if (bc_dirichlet[l][m][d])
                            PetscSNPrintf(cell[d], sizeof(cell[d]), vfmt[d],
                                          (double)bc_value[l][m][d], unit[d]);
                        else
                            PetscSNPrintf(cell[d], sizeof(cell[d]), "Neumann  (%s)",
                                          (d == 1) ? "insulating" : "sealed");
                    }
                    PetscPrintf(PETSC_COMM_WORLD, "   %-18s %-22s %-24s %s\n",
                                name, cell[0], cell[1], cell[2]);
                }
            }
            /* phi_i has no Dirichlet path anywhere in the code, so the column
             * above is constant by construction; call that out rather than
             * leaving the reader to wonder whether a flag could change it. */
            /* phi_i still has no Dirichlet path anywhere in the code. What it
             * DOES have, since -wall_faces, is a natural condition that is no
             * longer always zero: a regolith face carries the wall free-energy
             * term and enforces dphi/dn = cos(theta)*phi(1-phi)/eps. */
            if (user.wall_any && user.costhet != 0.0) {
                PetscPrintf(PETSC_COMM_WORLD,
                    "\n   phi_i is never pinned (no Dirichlet path exists for it). On a\n"
                    "   face marked \"regolith\" above, its NATURAL condition is the wall\n"
                    "   free-energy term  dphi/dn = cos(theta)*phi(1-phi)/eps  rather than\n"
                    "   zero; every other face keeps dphi/dn = 0 (a 90° contact angle).\n");
                PetscPrintf(PETSC_COMM_WORLD,
                    "   gamma_ia = %.4e   gamma_is = %.4e   gamma_as = %.4e  J/m²\n",
                    (double)user.gamma_ia, (double)user.gamma_is,
                    (double)user.gamma_as);
                if (user.costhet_direct)
                    PetscPrintf(PETSC_COMM_WORLD,
                        "   cos(theta) = %.6f  ->  theta = %.2f°   "
                        "[-contact_angle_deg: Young's equation BYPASSED, debug only]\n",
                        (double)user.costhet,
                        (double)(PetscAcosReal(user.costhet) * 180.0 / PETSC_PI));
                else
                    PetscPrintf(PETSC_COMM_WORLD,
                        "   cos(theta) = (gamma_as - gamma_is)/gamma_ia = %.6f"
                        "  ->  theta = %.2f°\n",
                        (double)user.costhet,
                        (double)(PetscAcosReal(user.costhet) * 180.0 / PETSC_PI));
            } else {
                PetscPrintf(PETSC_COMM_WORLD,
                    "\n   phi_i is natural Neumann on every wall — dphi/dn = 0, i.e. a 90°\n"
                    "   contact angle. Set -wall_faces (with -gamma_is/-gamma_as) to\n"
                    "   prescribe a different angle where ice meets regolith.\n");
            }
            if (flag_BC_Tfix && !bc_dirichlet[0][0][1] && !bc_dirichlet[1][0][1])
                PetscPrintf(PETSC_COMM_WORLD,
                    "   NOTE: -flag_BC_Tfix is set but no wall was pinned.\n");
            if (flag_BC_rhovfix && user.rhovfix_axis >= 0)
                PetscPrintf(PETSC_COMM_WORLD,
                    "   Vapor reservoir restricted to axis %d by -rhovfix_axis "
                    "(lo x%.4g, hi x%.4g of rho_vs(T0)).\n",
                    (int)user.rhovfix_axis, (double)user.rhovfix_lo,
                    (double)user.rhovfix_hi);
        }
    }

    PetscPrintf(PETSC_COMM_WORLD,
        "\n================================================================================\n\n");

    /* Create solution vector (ice, temperature, vapor) */
    Vec U;
    ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);
    ierr = VecZeroEntries(U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Setting up initial conditions... \n");

    /* ---- Restart from a previous run's snapshot -------------------------
     * -initial_cond <sol_NNNNN.dat> loads a solution vector written by
     * IGAWriteVec instead of building an initial condition. The option had
     * existed since before the fork but was dead: it was parsed into
     * user.initial_cond and never read, and nothing in the tree called
     * IGAReadVec, so passing it silently did nothing and the run started from
     * the -ic_type geometry as usual.
     *
     * The vector carries no mesh of its own, so the IGA here must match the one
     * that wrote it -- same -Nx/-Ny, -p, -C, -dof and domain. IGAReadVec checks
     * the length and errors on a mismatch, which catches the common case of
     * restarting against the wrong geometry file.
     *
     * -t_start sets the clock so a continuation reports absolute time rather
     * than starting again from zero; without it the second leg's t would
     * overlap the first and any velocity measured across the join would be
     * wrong. */
    if (user.initial_cond[0] != '\0') {
        PetscPrintf(PETSC_COMM_WORLD,
            "  RESTART: loading solution from %s\n"
            "           clock resumes at t = %.6e s (%.2f days)\n",
            user.initial_cond, (double)t_start, (double)(t_start / 86400.0));
        ierr = IGAReadVec(iga, U, user.initial_cond); CHKERRQ(ierr);
        user.readFlag = PETSC_TRUE;
        {   /* report what was actually loaded, so a wrong file is obvious */
            PetscReal pmin, pmax;
            ierr = VecStrideMin(U, 0, NULL, &pmin); CHKERRQ(ierr);
            ierr = VecStrideMax(U, 0, NULL, &pmax); CHKERRQ(ierr);
            PetscPrintf(PETSC_COMM_WORLD,
                "           phi range [%.6f, %.6f]\n", (double)pmin, (double)pmax);
        }
    } else if (dim == 1) {
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
        } else if (strcmp(ic_type, "multi_grains_file") == 0) {
            /* Same builder, but the grain list comes from -grains_file rather
             * than the -ice_grain_* option arrays (packings have ~400-500
             * grains, far past what is reasonable to inline in a .opts). */
            ierr = ReadGrainsFromFile(&user);               CHKERRQ(ierr);
            ierr = FormInitialMultiGrains2D(iga, U, &user); CHKERRQ(ierr);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Unknown -ic_type. Valid: two_ice_grains_boundary ice_slab "
                    "single_ice multi_grains multi_grains_file");
        }
    }

    /* ---- -test_wall_jacobian: analytic Jacobian vs finite differences ------
     * -snes_test_jacobian would do this, but only from inside a TS step, where
     * TSALPHA's restart solve re-enters the test repeatedly and a debug PETSc
     * build takes minutes on even a tiny mesh. This checks the same thing
     * directly and in seconds: for random directions v,
     *     J*v  ==  [F(U + h v) - F(U - h v)] / (2h)
     * to second order in h. Run it with and without -wall_faces to isolate the
     * boundary block -- the interior blocks are identical between the two, so
     * any discrepancy that appears only with -wall_faces is the wall term's.
     *
     * shift = 0 and V = 0, so this tests dR/dU alone. The wall term has no
     * phi_t dependence, which is exactly the part of J it contributes to. */
    if (test_wall_jacobian) {
        Mat Jfull, Jint;
        Vec Vz, Fp, Fm, Gp, Gm, v, Jv, Up;
        PetscReal h = 1.0e-7, worst_full = 0.0, worst_bnd = 0.0, scale_bnd = 0.0;
        PetscRandom rnd;
        const PetscReal cos_save = user.costhet;

        ierr = IGACreateMat(iga, &Jfull); CHKERRQ(ierr);
        ierr = IGACreateMat(iga, &Jint);  CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Vz); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Fp); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Fm); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Gp); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Gm); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &v);  CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Jv); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Up); CHKERRQ(ierr);
        ierr = VecZeroEntries(Vz); CHKERRQ(ierr);

        /* J with the wall term, and J with it switched off. Setting costhet = 0
         * makes WallPointActive() return false, so Jint is the interior form
         * alone -- every other coefficient is untouched. Their difference is
         * exactly the wall block. */
        ierr = IGAComputeIJacobian(iga, 0.0, Vz, 0.0, U, Jfull); CHKERRQ(ierr);
        user.costhet = 0.0;
        ierr = IGAComputeIJacobian(iga, 0.0, Vz, 0.0, U, Jint); CHKERRQ(ierr);
        user.costhet = cos_save;

        ierr = PetscRandomCreate(PETSC_COMM_WORLD, &rnd); CHKERRQ(ierr);
        ierr = PetscRandomSetType(rnd, PETSCRAND48); CHKERRQ(ierr);
        ierr = PetscRandomSetInterval(rnd, -1.0, 1.0); CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD,
            "\n WALL JACOBIAN CHECK (-test_wall_jacobian)   h = %.1e,  "
            "cos(theta) = %.6f\n"
            "   The full-system column is dominated by the interior form, so it\n"
            "   is insensitive to the wall block; the isolated column differences\n"
            "   both J and F between costhet on and off and is the real gate.\n\n"
            "   dir   full system      wall block       ||Jbnd*v||\n",
            (double)h, (double)user.costhet);

        for (PetscInt k = 0; k < 5; k++) {
            PetscReal nd_f, nj_f, nd_b, nj_b;

            ierr = PetscRandomSetSeed(rnd, (unsigned long)(12345 + k)); CHKERRQ(ierr);
            ierr = PetscRandomSeed(rnd); CHKERRQ(ierr);
            ierr = VecSetRandom(v, rnd); CHKERRQ(ierr);

            /* F(U +- h v) with the wall on (F) and off (G). */
            ierr = VecWAXPY(Up,  h, v, U); CHKERRQ(ierr);
            ierr = IGAComputeIFunction(iga, 0.0, Vz, 0.0, Up, Fp); CHKERRQ(ierr);
            user.costhet = 0.0;
            ierr = IGAComputeIFunction(iga, 0.0, Vz, 0.0, Up, Gp); CHKERRQ(ierr);
            user.costhet = cos_save;

            ierr = VecWAXPY(Up, -h, v, U); CHKERRQ(ierr);
            ierr = IGAComputeIFunction(iga, 0.0, Vz, 0.0, Up, Fm); CHKERRQ(ierr);
            user.costhet = 0.0;
            ierr = IGAComputeIFunction(iga, 0.0, Vz, 0.0, Up, Gm); CHKERRQ(ierr);
            user.costhet = cos_save;

            /* --- full system --- */
            ierr = VecCopy(Fp, Up); CHKERRQ(ierr);
            ierr = VecAXPY(Up, -1.0, Fm); CHKERRQ(ierr);
            ierr = VecScale(Up, 1.0 / (2.0 * h)); CHKERRQ(ierr);
            ierr = MatMult(Jfull, v, Jv); CHKERRQ(ierr);
            ierr = VecNorm(Jv, NORM_2, &nj_f); CHKERRQ(ierr);
            ierr = VecAXPY(Up, -1.0, Jv); CHKERRQ(ierr);
            ierr = VecNorm(Up, NORM_2, &nd_f); CHKERRQ(ierr);

            /* --- wall block alone: (F - G) differenced, vs (Jfull - Jint)*v --- */
            ierr = VecAXPY(Fp, -1.0, Gp); CHKERRQ(ierr);   /* wall part at +h */
            ierr = VecAXPY(Fm, -1.0, Gm); CHKERRQ(ierr);   /* wall part at -h */
            ierr = VecAXPY(Fp, -1.0, Fm); CHKERRQ(ierr);
            ierr = VecScale(Fp, 1.0 / (2.0 * h)); CHKERRQ(ierr);
            ierr = MatMult(Jfull, v, Jv); CHKERRQ(ierr);
            ierr = MatMult(Jint,  v, Gp); CHKERRQ(ierr);
            ierr = VecAXPY(Jv, -1.0, Gp); CHKERRQ(ierr);   /* Jbnd * v */
            ierr = VecNorm(Jv, NORM_2, &nj_b); CHKERRQ(ierr);
            ierr = VecAXPY(Fp, -1.0, Jv); CHKERRQ(ierr);
            ierr = VecNorm(Fp, NORM_2, &nd_b); CHKERRQ(ierr);

            if (nj_f > 0.0 && nd_f / nj_f > worst_full) worst_full = nd_f / nj_f;
            if (nj_b > 0.0 && nd_b / nj_b > worst_bnd)  worst_bnd  = nd_b / nj_b;
            if (nj_b > scale_bnd) scale_bnd = nj_b;

            PetscPrintf(PETSC_COMM_WORLD, "   %3d   %.6e    %.6e    %.6e\n",
                        (int)k,
                        (double)(nj_f > 0.0 ? nd_f / nj_f : nd_f),
                        (double)(nj_b > 0.0 ? nd_b / nj_b : nd_b),
                        (double)nj_b);
        }
        PetscPrintf(PETSC_COMM_WORLD,
            "   worst %.6e    %.6e\n"
            "   (||Jbnd*v|| = 0 would mean the wall block is never exercised --\n"
            "    the interface must actually cross a face named by -wall_faces.)\n",
            (double)worst_full, (double)worst_bnd);

        ierr = PetscRandomDestroy(&rnd); CHKERRQ(ierr);
        ierr = MatDestroy(&Jfull); CHKERRQ(ierr);
        ierr = MatDestroy(&Jint);  CHKERRQ(ierr);
        ierr = VecDestroy(&Vz); CHKERRQ(ierr);
        ierr = VecDestroy(&Fp); CHKERRQ(ierr);
        ierr = VecDestroy(&Fm); CHKERRQ(ierr);
        ierr = VecDestroy(&Gp); CHKERRQ(ierr);
        ierr = VecDestroy(&Gm); CHKERRQ(ierr);
        ierr = VecDestroy(&v);  CHKERRQ(ierr);
        ierr = VecDestroy(&Jv); CHKERRQ(ierr);
        ierr = VecDestroy(&Up); CHKERRQ(ierr);
        goto cleanup;
    }

    /* ---- -test_wall_measure: verify the boundary surface measure ----------
     * The one part of the wall term that cannot be checked by reading the code
     * is whether PetIGA integrates the boundary form with the right dS. On a
     * plain Cartesian patch (no -geom_file) petigaelem.c takes the detS = 1.0
     * branch rather than computing a geometric surface Jacobian, so the face
     * measure comes entirely from the surviving axes' quadrature weights.
     *
     * Set phi = 1/2 uniformly, T = temp0, rho_v = rho_vs(temp0), and phi_t = 0.
     * Every interior contribution to R[.][0] then vanishes IDENTICALLY:
     *   phi_t = 0; grad phi = 0; f1(1/2) = 0; and rho_v - rho_vs = 0 kills the
     *   sublimation source. Whatever is left in the assembled vector is the
     *   boundary term alone. With cos(theta) = 1 it integrates to
     *       sum_a F[a][0] = -3*M*(1/4)*|Gamma_wall| = -0.75*M*|Gamma_wall|
     * because the shape functions are a partition of unity. */
    if (test_wall_measure) {
        Vec Uw, Vw, Fw;
        PetscReal rho_vs_w, area = 0.0, got, want;
        const PetscReal LL[3] = {user.Lx, user.Ly, user.Lz};

        if (!user.wall_any)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                    "-test_wall_measure needs -wall_faces (and a non-zero "
                    "contact angle) or there is no boundary term to measure.");

        RhoVS_I(&user, user.temp0, &rho_vs_w, NULL);

        ierr = IGACreateVec(iga, &Uw); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Vw); CHKERRQ(ierr);
        ierr = IGACreateVec(iga, &Fw); CHKERRQ(ierr);
        ierr = VecZeroEntries(Vw); CHKERRQ(ierr);
        ierr = VecStrideSet(Uw, 0, 0.5); CHKERRQ(ierr);
        ierr = VecStrideSet(Uw, 1, user.temp0); CHKERRQ(ierr);
        ierr = VecStrideSet(Uw, 2, rho_vs_w); CHKERRQ(ierr);

        ierr = IGAComputeIFunction(iga, 0.0, Vw, 0.0, Uw, Fw); CHKERRQ(ierr);
        ierr = VecStrideSum(Fw, 0, &got); CHKERRQ(ierr);

        /* |Gamma_wall|: each flagged face contributes the measure of the
         * domain face perpendicular to its axis. */
        for (PetscInt l = 0; l < dim; l++) {
            PetscReal face = 1.0;
            for (PetscInt k = 0; k < dim; k++) if (k != l) face *= LL[k];
            for (PetscInt m = 0; m < 2; m++) if (user.wall_face[l][m]) area += face;
        }
        want = -3.0 * user.mob_sub * 0.25 * user.costhet * area;

        PetscPrintf(PETSC_COMM_WORLD,
            "\n WALL SURFACE-MEASURE CHECK (-test_wall_measure)\n"
            "   |Gamma_wall|   = %.12e m^%d\n"
            "   cos(theta)     = %.12f\n"
            "   sum F[.][0]    = %.12e   (assembled)\n"
            "   expected       = %.12e   (-3*M*phi(1-phi)*cos(theta)*|Gamma|)\n"
            "   rel. error     = %.3e\n",
            (double)area, (int)(dim - 1), (double)user.costhet,
            (double)got, (double)want,
            (double)(want != 0.0 ? PetscAbsReal((got - want) / want)
                                 : PetscAbsReal(got)));

        ierr = VecDestroy(&Uw); CHKERRQ(ierr);
        ierr = VecDestroy(&Vw); CHKERRQ(ierr);
        ierr = VecDestroy(&Fw); CHKERRQ(ierr);
        goto cleanup;
    }

    /* Solve the system */
    ierr = TSSolve(ts, U); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Solution completed. \n");

cleanup:
    /* Cleanup Resources */
    if (user.ssa_view) { ierr = PetscViewerDestroy(&user.ssa_view); CHKERRQ(ierr); }
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = VecDestroy(&user.cfl_U_prev); CHKERRQ(ierr);
    ierr = TSDestroy(&ts); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    for (PetscInt d = 0; d < 3; d++) { ierr = PetscFree(user.cent[d]); CHKERRQ(ierr); }
    ierr = PetscFree(user.radius);       CHKERRQ(ierr);
    ierr = PetscFree(user.ice_grain_ax); CHKERRQ(ierr);
    ierr = PetscFree(user.ice_grain_ay); CHKERRQ(ierr);
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
