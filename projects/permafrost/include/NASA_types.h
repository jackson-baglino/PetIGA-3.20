#ifndef NASA_TYPES_H
#define NASA_TYPES_H

#include <petsc.h>
#include "petiga.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/* Field definitions for node data */
typedef struct {
  PetscScalar soil;
} FieldS;
typedef struct {
  PetscScalar ice, tem, rhov, sed;
} Field;
/* Application context structure */
typedef struct {
  IGA       iga;  // Isogeometric analysis (IGA) structure for managing geometry and basis functions

  // Physical parameters related to phase field and thermodynamics
  PetscReal eps;  // Interface width parameter for phase field method
  PetscReal mob_sub;  // Mobility for ice phase evolution
  PetscReal mob_sed;  // Mobility for sediment phase evolution equation
  PetscReal mob_air; // Mobility for air phase evolution equation (if applicable)
  PetscReal Etai, Etam, Etaa;  // Activation energy terms for different phases (ice, sediment, air)
  PetscReal alph_sub;  // Substrate interaction coefficient
  PetscReal Lambd;  // Parameter related to thermal conductivity or latent heat (context-dependent)
  PetscReal beta_sub0, d0_sub0;  // Parameters related to phase change at the substrate

  // Thermophysical properties of different phases
  PetscReal thcond_ice, thcond_sed, thcond_air;  // Thermal conductivities of ice, sediment, and air
  PetscReal cp_ice, cp_sed, cp_air;  // Specific heat capacities of ice, sediment, and air
  PetscReal rho_ice, rho_sed, rho_air;  // Densities of ice, sediment, and air
  PetscReal dif_vap;  // Vapor diffusivity in air
  PetscReal lat_sub;  // Latent heat of sublimation
  PetscReal diff_sub;  // Diffusivity related to sublimation

  // Environmental conditions and threshold parameters
  PetscReal air_lim;  // Air phase fraction limit (threshold to distinguish between ice/air)
  PetscReal xi_v, xi_T;  // Characteristic non-dimensional parameters for vapor and temperature

  // Initial and boundary condition parameters
  PetscReal T_melt;  // Melting temperature of ice
  PetscReal temp0, hum0;  // Initial temperature and humidity
  PetscReal grad_temp0[3];  // Initial temperature gradient in x, y, and z directions

  // Domain size and resolution
  PetscReal Lx, Ly, Lz;  // Physical domain dimensions in x, y, and z
  PetscInt  Nx, Ny, Nz;  // Number of elements in x, y, and z directions

  // Radius of curvature parameters (possibly for computing capillary effects)
  PetscReal RCice, RCsed;  // Mean radius of curvature for ice and sediment grains
  PetscReal RCice_dev, RCsed_dev;  // Standard deviation of radius of curvature for ice and sediment

  // Per-grain radii and separation for enclosed grain pair IC
  PetscReal RCice0, RCice1;   /* Outer ice radius of grain 0 / grain 1 (default: RCice) */
  PetscReal RCsed0, RCsed1;   /* Sediment core radius of grain 0 / grain 1 (default: RCsed) */
  PetscReal grain_sep;         /* Air gap between outer ice surfaces (m); 0 = tangent */
  PetscReal x_slab_frac;      /* Fraction of Lx occupied by the right-side ice slab (slab_and_grains IC) */

  // Arrays storing geometry information for ice and sediment grains
  PetscReal cent[3][200];  // Coordinates of ice grain centers (3D array for x, y, z positions)
  PetscReal radius[200];  // Radii of individual ice grains
  PetscReal centsed[3][200];  // Coordinates of sediment grain centers (3D array for x, y, z positions)
  PetscReal radiussed[200];  // Radii of individual sediment grains

  // Initial normal vector components (possibly for a structured interface)
  PetscReal norm0[4];  // Per-DOF initial residual norms for SNES convergence check

  // Flags for controlling different simulation options
  PetscInt  flag_tIC;        // IC geometry variant: 0=centered slab, 2=flat interface
  PetscInt  outp;            // output control flag
  PetscInt  n_relax;         // AC-only relaxation steps before full physics (0 = none)
  PetscBool flag_relax;      // PETSC_TRUE while relaxation steps are active
  PetscBool flag_Tdep;       // temperature-dependent material properties
  PetscBool flag_sed_frozen; // PETSC_FALSE = 3-phase active; PETSC_TRUE = 2-phase (sediment frozen)

  /* Residual avenue selector:
   *   1 = Allen-Cahn, penalty vapor, freeze sed → zero RHS after t_sed_freeze
   *   2 = Allen-Cahn, penalty vapor, freeze sed → k_sed penalty (default)
   *   3 = Cahn-Hilliard, no penalties (requires p ≥ 2, C ≥ 1) */
  PetscInt  flag_avenue;

  PetscReal t_sed_freeze;    // duration of 3-phase period (s); 0 = start immediately in 2-phase

  // Numerical method and discretization parameters
  PetscInt p;  // Polynomial degree of basis functions (for IGA)
  PetscInt C;  // Continuity of basis functions
  PetscInt dim;  // Spatial dimension of the problem (2D or 3D)
  PetscInt dof;  // Degrees of freedom per node (ice, temperature, vapor, sediment)
  PetscInt periodic;  // Periodicity flag (0 = non-periodic, 1 = periodic boundaries)

  // Time stepping parameters
  PetscReal t_out;  // Output time interval
  PetscReal t_interv;  // Intermediate time step interval
  PetscReal t_IC;  // Total duration for initial condition phase

  // Counters for active ice and sediment grains
  PetscInt NCice, NCsed;  // Number of ice and sediment grains
  PetscInt n_act, n_actsed;  // Number of currently active grains (ice and sediment)

  // Arrays for field variables
  PetscScalar *Phi_sed0;  // Sediment phase field variable
  PetscInt npoints;

  /* Penalty parameters (tunable via CLI) */
  PetscReal difvap_pen;   // multiplicative factor: D_pen = difvap_pen * difvap (dimensionless)
  PetscReal k_pen;        // vapour interface equilibrium stiffness
  PetscReal k_sed_pen;    // sediment shape-restoring stiffness

  /* Phase-field bounds: simulation aborts if any phi leaves [phase_lo, phase_hi] */
  PetscReal phase_lo;     // lower bound for phi_ice, phi_sed, phi_air (default -0.25)
  PetscReal phase_hi;     // upper bound for phi_ice, phi_sed, phi_air (default  1.25)
  PetscReal *alph;     // Alpha field, possibly phase fraction or related property
  PetscReal *mob;      // Ice mobility field, spatially varying (T-dependent)

  // Flag for reading input files
  PetscBool readFlag; // read initial field data from file

  // Output file path
  char output_path[PETSC_MAX_PATH_LEN];  // Path for output files
  char initial_cond[PETSC_MAX_PATH_LEN];  // Path for initial condition file
  char initial_PFgeom[PETSC_MAX_PATH_LEN];  // Path for initial geometry file

  // Capillary neck parameters
  PetscReal R1;  // Radius of capillary neck

  // Initial domain integrals (set at step 0, used for percentage reporting)
  PetscReal tot_ice_0;     // initial ∫ φ_i dΩ
  PetscReal tot_air_0;     // initial ∫ φ_a dΩ
  PetscReal tot_sed_0;     // initial ∫ φ_s dΩ
  PetscReal tot_rhov_0;    // initial ∫ ρ_v φ_a dΩ (vapor mass in air phase)
  PetscReal tot_mass_0;    // initial ρ_ice·∫φ_i + ρ_sed·∫φ_s + ∫ρ_v·φ_a (total system mass)

} AppCtx;/* Field definitions for node data */

#endif // NASA_TYPES_H