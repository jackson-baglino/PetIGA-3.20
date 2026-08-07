#ifndef NASA_TYPES_H
#define NASA_TYPES_H

#include <petsc.h>
#include "petiga.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/* Field definitions for node data */
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
  // With -flag_BC_rhovfix: which axis's two faces get the Dirichlet vapor
  // reservoir. -1 (default) = every face, the legacy behaviour. A wall-bounded
  // pore channel wants 0 (x-faces only): its top/bottom are solid regolith, and
  // pinning vapor there feeds the ice at its contact line.
  PetscInt  rhovfix_axis;
  // Per-face vapor reservoir strength, as a multiple of rho_vs(temp0): lo is
  // the m=0 face (x=0, "left"), hi the m=1 face (x=Lx, "right"). Both default
  // to hum0, which reproduces the previous single-value behaviour exactly.
  //
  // These must be specified at the 1e-6 level, NOT as a humidity. A grain's own
  // Gibbs-Thomson equilibrium sits within ~5e-6 of rho_vs (d0/r), so a
  // humidity-style 0.99 is ~2000x the grain scale and simply sublimates
  // everything uniformly, hiding the curvature physics the run is measuring.
  PetscReal rhovfix_lo, rhovfix_hi;

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

  // Affine baseline of each wall, added UNDER the bumps above (see
  // WallBottom()/WallTop() in src/initial_conditions.c):
  //     y_bot(x) = wall_bot_y0 + wall_bot_slope*x + SedimentBumpField(x)
  //     y_top(x) = wall_top_y0 + wall_top_slope*x - TopBumpField(x)
  // Bumps have compact support and cannot express a linear ramp, so a
  // tapered (wedge) domain needs these. Must match the --bot-y0/--bot-slope/
  // --top-y0/--top-slope arguments of build_geometry_multi_grain.py, which
  // cuts the mesh from the same two curves.
  // Defaults (0, 0, Ly, 0) reproduce the previous flat-baseline behaviour.
  PetscReal wall_bot_y0, wall_bot_slope;
  PetscReal wall_top_y0, wall_top_slope;

  // Wedge-bridging ice band: the annulus wedge_band_r1 <= |X - apex| <=
  // wedge_band_r2 about (wedge_apex_x, wedge_apex_y). An apex-centred arc is
  // perpendicular to every ray from the apex, i.e. to both wedge walls, so
  // this is the shape that meets both at the natural 90-degree contact angle.
  // Inactive unless wedge_band_r2 > wedge_band_r1 (both default 0).
  // n_wedge_bands annuli sharing one apex; band k spans wedge_band_r1[k] to
  // wedge_band_r2[k]. Two bands let a wedge hold a pair of grains at different
  // distances from the apex, i.e. at different confinement.
  PetscReal wedge_apex_x, wedge_apex_y;
  PetscInt  n_wedge_bands;
  PetscReal wedge_band_r1[MAX_SED_GRAINS], wedge_band_r2[MAX_SED_GRAINS];

  // Arrays storing geometry information for ice grains.
  //
  // Heap-allocated with capacity n_grain_max (-n_grain_max, default 2000).
  // These were fixed [200] arrays until 2026-07-27; granular packings from
  // preprocess/generate_packing.py routinely carry 400-500 grains at
  // 2 mm / 45 um mean radius, which silently overran them. Capacity is
  // allocated once in main() and freed at the end; indexing is unchanged.
  PetscInt  n_grain_max;   // Allocated capacity of the four arrays below
  PetscReal *cent[3];      // Ice grain centre coordinates: cent[dim][k]
  PetscReal *radius;       // Radii of individual ice grains (isotropic; default ax/ay)
  PetscReal *ice_grain_ax; /* ellipse semi-axis in x; defaults to radius[k] if -ice_grain_ax not set */
  PetscReal *ice_grain_ay; /* ellipse semi-axis in y; defaults to radius[k] if -ice_grain_ay not set */

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

  // Ice plug spanning a pore throat. Each plug is the INTERSECTION of two
  // meniscus discs, both centred ON the channel axis y=bridge_cy:
  //     left  meniscus: circle (bridge_cxL, cy) radius bridge_rL
  //     right meniscus: circle (bridge_cxR, cy) radius bridge_rR
  // Ice lies INSIDE both discs, so each meniscus is CONVEX and the plug is
  // BARREL-shaped -- wider at the axis than at its wall contacts. That is what
  // a 90-degree contact angle forces in a converging-diverging channel; the
  // hourglass shape of a wetting capillary bridge needs theta < 90, which this
  // model's Neumann wall BC cannot produce.
  //
  // 90 degrees requires the arc centre to lie along the wall's TANGENT at the
  // contact (NOT its normal, which gives tangency, i.e. a 0-degree contact):
  //     xc = xp + (cy - y_wall(xp)) / y_wall'(xp)
  // For a straight wall this returns the wedge apex, which is the cross-check
  // that the sign is right. The radii must be solved against the actual wall
  // shape -- generate them with preprocess/build_geometry_two_throat.py, which
  // asserts 90 degrees at every contact, and never hand-write them.
  //
  // Inactive when n_bridges == 0.
  PetscInt  n_bridges;
  PetscReal bridge_cxL[MAX_SED_GRAINS], bridge_rL[MAX_SED_GRAINS];
  PetscReal bridge_cxR[MAX_SED_GRAINS], bridge_rR[MAX_SED_GRAINS];
  PetscReal bridge_cy[MAX_SED_GRAINS];

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

  /* Anti-trapping / stabilization multipliers (Moure & Fu 2024, K&P 2009).
   * xi_T scales thermal conduction and latent heat in R_tem.
   * xi_v scales vapor diffusion and the rho_ice part of the mass exchange
   *      source in R_vap — reduces spurious interface fluxes from the
   *      finite-width diffuse interface to the physical sharp-interface limit.
   * Defaults: xi_T = 1.0 (no effect), xi_v = 1e-3 (M&F sublimation value). */
  PetscReal xi_T;
  PetscReal xi_v;

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
  char grains_file[PETSC_MAX_PATH_LEN];  // Grain list for -ic_type multi_grains_file

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

  // Interface-CFL timestep limiter (InterfaceCFLMonitor in monitoring.c):
  // clamps the NEXT dt so max pointwise |d(phi)| per step stays below
  // cfl_dphimax, measured from the last accepted step's rate
  // ||phi^n - phi^{n-1}||_inf / dt. Lets dtmax be large (quiet phases
  // cruise) while fast interface events (grain collapse, neck formation)
  // throttle dt automatically instead of Gibbs-rippling the B-spline front.
  PetscBool flag_dtCFL;      // -dtCFL (default 1)
  PetscReal cfl_dphimax;     // -dtCFL_dphimax (default 0.2)
  Vec       cfl_U_prev;      // previous accepted solution (lazy-created)
  PetscReal cfl_t_prev;      // time of previous accepted step

  // Axisymmetric (r-z) mode: x = z (symmetry axis direction), y = r (radial),
  // axis on the y = 0 boundary. Every residual/Jacobian/monitor integrand is
  // weighted by the quadrature point's r-coordinate (2*pi*r dr dz measure);
  // the azimuthal curvature mode is generated by that weight automatically
  // (see docs/axisymmetric_plan.md section 1b). Grains must be centered on
  // the axis (ice_grain_cy = 0).
  PetscBool axisym;          // -axisym (default 0 = planar)

  // -ic_grain_union: build the multi_grains IC from the signed distance to the
  // UNION of the grains (sdf = min_k sdf_k, phi = 0.5-0.5*tanh(0.5*sdf/eps))
  // instead of SUMMING each grain's tanh profile.
  //
  // Why it matters for OVERLAPPING grains: the additive form's phi=0.5 contour
  // depends on eps at the neck, because both grains' tails contribute there.
  // (Away from the neck a lone grain crosses 0.5 at r=R for any eps, which is
  // why only the neck misbehaves.) With the union form phi=0.5 <=> sdf=0, i.e.
  // exactly the sharp union surface, for ANY eps -- so an eps series shares one
  // initial contour instead of one per eps (38.19/35.60/32.91 um measured in
  // the 2026-07-15 series).
  //
  // Caveat: the sharp union has a zero-radius concave crease at the neck, which
  // no diffuse interface can hold; AC dynamics rounds it to O(eps) within a few
  // steps. Contours coincide at t=0 and then diverge by O(eps) at the neck --
  // that residual IS the eps error being measured.
  //
  // Default 0 (additive) preserves every pre-2026-07-15 run. The two forms
  // agree for well-separated grains, where the far grain's tanh is negligible.
  PetscBool ic_grain_union;  // -ic_grain_union (default 0 = additive)

  // Persistent viewer for SSA_evo.dat. Opened ONCE on the first Monitor() call
  // and flushed every step, instead of re-opening a viewer per step (which
  // silently exhausted file descriptors ~step 15 -- the calls had no CHKERRQ --
  // truncating the file while outp.txt kept going). NULL until first use.
  PetscViewer ssa_view;

} AppCtx;/* Field definitions for node data */

#endif // NASA_TYPES_H