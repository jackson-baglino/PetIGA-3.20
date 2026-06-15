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
  PetscScalar ice, tem, rhov;
} Field;
/* Application context structure */
typedef struct {
  IGA       iga;  // Isogeometric analysis (IGA) structure for managing geometry and basis functions
  SNES      snes; // Nonlinear solver handle (cached so Residual can call SNESSetFunctionDomainError)

  // Physical parameters related to phase field and thermodynamics
  PetscReal eps;  // Interface width parameter for phase field method
  PetscReal mob_sub;  // Mobility for ice phase evolution
  PetscReal Etai, Etam, Etaa;  // Surface energy terms: Sigma_i (ice-vapor), Etam (unused), Sigma_a (air-vapor side)
  PetscReal alph_sub;  // Substrate interaction coefficient
  PetscReal Lambd;  // Higher-order penalty coefficient Lambda in the double-well free energy
  PetscReal beta_sub0, d0_sub0;  // Parameters related to phase change at the substrate

  // Thermophysical properties of different phases
  PetscReal thcond_ice, thcond_air;  // Thermal conductivities of ice and air
  PetscReal cp_ice, cp_air;  // Specific heat capacities of ice and air
  PetscReal rho_ice, rho_air;  // Densities of ice and air
  PetscReal dif_vap;  // Vapor diffusivity in air
  PetscReal lat_sub;  // Latent heat of sublimation
  PetscReal diff_sub;  // Diffusivity related to sublimation

  // Environmental conditions and threshold parameters
  PetscReal air_lim;  // Air phase fraction floor used in the vapor equation (avoids division/degeneracy as phi_a -> 0)

  // Initial and boundary condition parameters
  PetscReal T_melt;  // Melting temperature of ice
  PetscReal temp0, hum0;  // Initial temperature and humidity
  PetscReal grad_temp0[3];  // Initial temperature gradient in x, y, and z directions

  // Domain size and resolution
  PetscReal Lx, Ly, Lz;  // Physical domain dimensions in x, y, and z
  PetscInt  Nx, Ny, Nz;  // Number of elements in x, y, and z directions

  // Radius of curvature parameters (possibly for computing capillary effects)
  PetscReal RCice;  // Mean radius of curvature for ice grains
  PetscReal RCice_dev;  // Standard deviation of radius of curvature for ice grains

  // Per-grain radii for boundary-centered two-grain IC
  PetscReal RCice0, RCice1;   /* Outer ice radius of grain 0 / grain 1 (default: RCice) */

  // Sediment-grain "bump" geometry parameter (must match -geom_file's
  // build_geometry_sediment_grain.py R_sed; 0 => flat domain, no distortion)
  PetscReal geom_bump_R;

  // Arrays storing geometry information for ice grains
  PetscReal cent[3][200];  // Coordinates of ice grain centers (3D array for x, y, z positions)
  PetscReal radius[200];  // Radii of individual ice grains

  // Initial normal vector components (possibly for a structured interface)
  PetscReal norm0[3];  // Per-DOF initial residual norms for SNES convergence check

  // Flags for controlling different simulation options
  PetscInt  flag_tIC;        // IC geometry variant: 0=centered slab, 2=flat interface
  PetscInt  outp;            // output control flag
  PetscBool flag_Tdep;       // temperature-dependent material properties

  // Numerical method and discretization parameters
  PetscInt p;  // Polynomial degree of basis functions (for IGA)
  PetscInt C;  // Continuity of basis functions
  PetscInt dim;  // Spatial dimension of the problem (2D or 3D)
  PetscInt dof;  // Degrees of freedom per node (ice, temperature, vapor)
  PetscInt periodic;  // Periodicity flag (0 = non-periodic, 1 = periodic boundaries)

  // Time stepping parameters
  PetscReal t_out;  // Output time interval
  PetscReal t_interv;  // Intermediate time step interval
  PetscReal t_IC;  // Total duration for initial condition phase

  // Counters for active ice grains
  PetscInt NCice;  // Number of ice grains
  PetscInt n_act;  // Number of currently active ice grains

  PetscInt npoints;

  /* Phase-field bounds: simulation aborts if any phi leaves [phase_lo, phase_hi] */
  PetscReal phase_lo;     // lower bound for phi_ice, phi_air (default -0.25)
  PetscReal phase_hi;     // upper bound for phi_ice, phi_air (default  1.25)
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
  PetscReal tot_rhov_0;    // initial ∫ ρ_v φ_a dΩ (vapor mass in air phase)
  PetscReal tot_mass_0;    // initial ρ_ice·∫φ_i + ∫ρ_v·φ_a (total system mass)

  // Deferred bounds-rollback request — set by Monitor() when phase fields go
  // out of bounds, consumed by BoundsRollbackPreStep() before the next TSStep.
  // We can't call TSRollBack inside Monitor because ts->vec_sol is read-locked
  // there; the PreStep callback runs when the vector is writable.
  PetscBool bounds_violated;
  PetscReal bounds_new_dt;

} AppCtx;/* Field definitions for node data */

#endif // NASA_TYPES_H