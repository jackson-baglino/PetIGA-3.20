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
  PetscReal Lambd;  // Parameter related to thermal conductivity or latent heat (context-dependent)
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

  // Multi-grain sediment "bump" geometry: the bottom edge of a -geom_file
  // bump geometry is the sum of n_sed_grains individual SedimentBump()
  // humps. Each bump has a center (sed_grain_x), half-width (sed_grain_R),
  // and peak height (sed_grain_h, defaults to sed_grain_R if not provided).
  // Must match build_geometry_multi_grain.py's SEDIMENT_GRAINS list.
  // n_sed_grains == 0 falls back to the single-bump geom_bump_R behavior.
#define MAX_SED_GRAINS 24
  PetscInt  n_sed_grains;
  PetscReal sed_grain_x[MAX_SED_GRAINS];
  PetscReal sed_grain_R[MAX_SED_GRAINS];
  PetscReal sed_grain_h[MAX_SED_GRAINS];  /* peak height; default = R */

  // Top-wall (ceiling) bumps — same C∞ shape, push DOWN from Ly.
  // Must match TOP_GRAINS in build_geometry_multi_grain.py.
  PetscInt  n_top_grains;
  PetscReal top_grain_x[MAX_SED_GRAINS];
  PetscReal top_grain_R[MAX_SED_GRAINS];
  PetscReal top_grain_h[MAX_SED_GRAINS];

  // Arrays storing geometry information for ice grains
  PetscReal cent[3][200];  // Coordinates of ice grain centers (3D array for x, y, z positions)
  PetscReal radius[200];   // Radii of individual ice grains (isotropic; used as default ax/ay)
  PetscReal ice_grain_ax[200]; /* ellipse semi-axis in x; defaults to radius[k] if -ice_grain_ax not set */
  PetscReal ice_grain_ay[200]; /* ellipse semi-axis in y; defaults to radius[k] if -ice_grain_ay not set */

  // Ice "shell" capping a floor bump at constant thickness, conformal to
  // the bump's own surface (true distance to the SedimentBumpField(x)
  // curve, not just a vertical offset -- see SedimentBumpFieldDeriv() and
  // the -ice_shell_x loop in FormInitialMultiGrains2D), windowed laterally
  // to [ice_shell_x[k]-ice_shell_R[k], ice_shell_x[k]+ice_shell_R[k]] so it
  // only covers the bump itself, not the whole floor. Added on top of the
  // ice_grain_* ellipses, not a replacement.
  PetscInt  n_ice_shells;
  PetscReal ice_shell_x[MAX_SED_GRAINS];
  PetscReal ice_shell_R[MAX_SED_GRAINS];
  PetscReal ice_shell_thickness[MAX_SED_GRAINS];

  // Flat ice layer encapsulating a floor bump: ice fills everything below
  // the ABSOLUTE height ice_flat_height[k] (not relative to the bump's own
  // surface, unlike ice_shell_*), windowed laterally to
  // [ice_flat_x[k]-ice_flat_R[k], ice_flat_x[k]+ice_flat_R[k]]. Gives a flat
  // (non-rounded) ice-air interface burying the bump, instead of a domed
  // cap or a conformal coating. Added on top of ice_grain_*/ice_shell_*.
  PetscInt  n_ice_flats;
  PetscReal ice_flat_x[MAX_SED_GRAINS];
  PetscReal ice_flat_R[MAX_SED_GRAINS];
  PetscReal ice_flat_height[MAX_SED_GRAINS];

  // Initial normal vector components (possibly for a structured interface)
  PetscReal norm0[3];  // Per-DOF initial residual norms for SNES convergence check

  // Flags for controlling different simulation options
  PetscInt  flag_tIC;        // IC geometry variant: 0=centered slab, 2=flat interface
  PetscInt  outp;            // output control flag
  PetscBool flag_Tdep;       // temperature-dependent material properties
  PetscBool decouple_phase_change;  // -decouple_phase_change 1: zero the ice_t-driven
                                     // source terms in R_tem (latent heat) and R_vap
                                     // (mass-balance) too, not just S_sub in R_ice --
                                     // isolates pure AC curvature relaxation from any
                                     // leakage into temperature/vapor (see assembly.c)

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

  /* Gibbs-Thomson capillary length: rhoI_vs_eff = rhoI_vs*(1 + d0_GT*kappa).
   * Defaults to the same microscopic capillary length used to derive the
   * kinetic coefficients (d0_sub, see permafrost2.c) -- 0 disables GT
   * (flat-interface rhoI_vs, no curvature-driven Ostwald ripening). */
  PetscReal d0_GT;
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