#include "material_properties.h"
#include "NASA_types.h"

/**
 * @brief Computes the effective thermal conductivity based on phase fractions.
 *
 * This function calculates the thermal conductivity of a material as a weighted
 * sum of the thermal conductivities of ice and air phases (air = 1 - ice). It
 * also computes the derivative of conductivity with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase (0 to 1).
 * @param cond Pointer to store the computed thermal conductivity.
 * @param dcond_ice Pointer to store the derivative of conductivity with respect to ice fraction.
 */
void ThermalCond(AppCtx *user, PetscScalar ice,
                 PetscScalar *cond, PetscScalar *dcond_ice)
{
    // Define derivatives (initialized to 1.0, modified if phase fractions are negative)
    PetscReal dice = 1.0, dair = 1.0;

    // Compute air fraction (air = 1 - ice)
    PetscReal air = 1.0 - ice;

    // Ensure phase fractions are non-negative (corrects numerical issues)
    if (ice < 0.0) {
        ice = 0.0;
        dice = 0.0; // If ice fraction is zero, its derivative should also be zero
    }
    if (air < 0.0) {
        air = 0.0;
        dair = 0.0; // If air fraction is zero, its derivative should also be zero
    }

    // Retrieve material thermal conductivities from the user context
    PetscReal cond_ice = user->thcond_ice;
    PetscReal cond_air = user->thcond_air;

    // Compute the effective thermal conductivity as a weighted sum
    if (cond)
        (*cond) = ice * cond_ice + air * cond_air;

    // Compute the derivative of thermal conductivity with respect to ice fraction
    if (dcond_ice)
        (*dcond_ice) = cond_ice * dice - cond_air * dair;

    return; // Explicit return statement for clarity
}

/**
 * @brief Computes the effective heat capacity and its derivative with respect to ice.
 *
 * This function calculates the heat capacity as a weighted sum of contributions
 * from ice and air (air = 1 - ice). It also computes the derivative of heat
 * capacity with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param cp Pointer to store computed heat capacity.
 * @param dcp_ice Pointer to store derivative of heat capacity with respect to ice.
 */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp,
  PetscScalar *dcp_ice)
{
PetscReal dice = 1.0, dair = 1.0;
PetscReal air = 1.0 - ice;

// Ensure phase fractions are non-negative
if (ice < 0.0) { ice = 0.0; dice = 0.0; }
if (air < 0.0) { air = 0.0; dair = 0.0; }

// Retrieve heat capacities
PetscReal cp_ice = user->cp_ice;
PetscReal cp_air = user->cp_air;

// Compute effective heat capacity
if (cp)
(*cp) = ice * cp_ice + air * cp_air;

// Compute derivative with respect to ice
if (dcp_ice)
(*dcp_ice) = cp_ice * dice - cp_air * dair;
return;
}

/**
 * @brief Computes the effective density and its derivative with respect to ice.
 *
 * This function calculates the density as a weighted sum of contributions
 * from ice and air (air = 1 - ice). It also computes the derivative of
 * density with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param rho Pointer to store computed density.
 * @param drho_ice Pointer to store derivative of density with respect to ice.
 */
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho,
  PetscScalar *drho_ice)
{
PetscReal dice = 1.0, dair = 1.0;
PetscReal air = 1.0 - ice;

// Ensure phase fractions are non-negative
if (ice < 0.0) { ice = 0.0; dice = 0.0; }
if (air < 0.0) { air = 0.0; dair = 0.0; }

// Retrieve densities
PetscReal rho_ice = user->rho_ice;
PetscReal rho_air = user->rho_air;

// Compute effective density
if (rho)
(*rho) = ice * rho_ice + air * rho_air;

// Compute derivative with respect to ice
if (drho_ice)
(*drho_ice) = rho_ice * dice - rho_air * dair;
return;
}

/**
 * @brief Computes vapor diffusivity and its derivative with respect to temperature.
 *
 * @param user Pointer to the application context containing material properties.
 * @param tem Temperature in Celsius.
 * @param difvap Pointer to store computed vapor diffusivity.
 * @param d_difvap Pointer to store derivative of vapor diffusivity with respect to temperature.
 */
void VaporDiffus(AppCtx *user, PetscScalar tem, PetscScalar *difvap,
  PetscScalar *d_difvap)
{
PetscReal dif_vap = user->dif_vap;
PetscReal Kratio = (tem + 273.15) / 273.15; // Convert temperature to Kelvin ratio
// PetscReal meanT = user->temp0, Kratio = (meanT + 273.15) / 273.15;
PetscReal aa = 1.81;

// Compute vapor diffusivity
if (difvap)
(*difvap) = dif_vap * pow(Kratio, aa);

// Compute derivative with respect to temperature
if (d_difvap)
(*d_difvap) = dif_vap * aa * pow(Kratio, aa - 1.0) / 273.15;

return;
}

/**
 * @brief Computes the saturation vapor density and its derivative with respect to temperature.
 *
 * @param user Pointer to the application context containing material properties.
 * @param tem Temperature in Celsius.
 * @param rho_vs Pointer to store computed saturation vapor density.
 * @param d_rhovs Pointer to store derivative of saturation vapor density with respect to temperature.
 */
void RhoVS_I(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs,
  PetscScalar *d_rhovs)
{
PetscReal rho_air = user->rho_air;
PetscReal Patm = 101325.0;
PetscReal bb = 0.62;
PetscReal temK = tem + 273.15; // Convert temperature to Kelvin

// Empirical coefficients for saturation pressure
PetscReal K0 = -0.5865e4, K1 = 0.2224e2, K2 = 0.1375e-1;
PetscReal K3 = -0.3403e-4, K4 = 0.2697e-7, K5 = 0.6918;

// Compute saturation vapor pressure
PetscReal Pvs = exp(K0 * pow(temK, -1.0) + K1 + K2 * temK +
             K3 * pow(temK, 2.0) + K4 * pow(temK, 3.0) +
             K5 * log(temK));

// Compute derivative of vapor pressure with respect to temperature
PetscReal Pvs_T = Pvs * (-K0 * pow(temK, -2.0) + K2 +
                  2.0 * K3 * temK + 3.0 * K4 * pow(temK, 2.0) +
                  K5 / temK);

// Compute saturation vapor density
if (rho_vs)
(*rho_vs) = rho_air * bb * Pvs / (Patm - Pvs);

// Compute derivative with respect to temperature
if (d_rhovs)
(*d_rhovs) = rho_air * bb * (Pvs_T * (Patm - Pvs) + Pvs * Pvs_T) /
          ((Patm - Pvs) * (Patm - Pvs));

return;
}

void Mobility(AppCtx *user, PetscScalar ice, PetscScalar *mob)
{
    /* Constant, phase-independent mobility. */
    (void)ice;  /* unused — kept in signature for callers */
    if (mob) (*mob) = user->mob_sub;
}


/**
 * @brief Computes σ₀ (sigma zero) as a function of temperature using lookup tables and
 * interpolation.
 *
 * This function determines the value of σ₀ based on predefined temperature values
 * using a logarithmic interpolation technique. If the temperature is out of range,
 * a warning is printed, and the function assigns the closest available value.
 *
 * @param temp Temperature in Celsius.
 * @param sigm0 Pointer to store the computed σ₀ value.
 */
void Sigma0(PetscScalar temp, PetscScalar *sigm0)
{
    // Lookup table for σ₀ values at different temperatures
    PetscReal sig[10], tem[10];

    // Predefined values for σ₀ corresponding to specific temperatures
    sig[0] = 3.0e-3;  sig[1] = 4.1e-3;  sig[2] = 5.5e-3;  sig[3] = 8.0e-3;  sig[4] = 4.0e-3;
    sig[5] = 6.0e-3;  sig[6] = 3.5e-2;  sig[7] = 7.0e-2;  sig[8] = 1.1e-1;  sig[9] = 0.75;

    // Corresponding temperature values (in Celsius)
    tem[0] = -0.0001;  tem[1] = -2.0;   tem[2] = -4.0;   tem[3] = -6.0;   tem[4] = -7.0;
    tem[5] = -10.0;    tem[6] = -20.0;  tem[7] = -30.0;  tem[8] = -40.0;  tem[9] = -100.0;

    PetscInt ii, interv = 0;
    PetscReal t0, t1, s0, s1;

    // Warn user if temperature is outside the valid range
    if (temp > tem[0] || temp < tem[9])
        PetscPrintf(PETSC_COMM_WORLD, "Warning: Temperature (%g) out of range in Sigma0 function.\n", temp);

    // Find the interval in which `temp` falls
    for (ii = 0; ii < 10; ii++) {
        if (temp <= tem[ii])
            interv = ii;
    }

    // If temperature is above the highest defined value, use the first table value
    if (temp > tem[0])
        interv = -1;

    // Assign σ₀ based on the identified interval
    if (interv == -1) {
        *sigm0 = sig[0];  // Assign σ₀ at highest temperature
    }
    else if (interv == 9) {
        *sigm0 = sig[9];  // Assign σ₀ at lowest temperature
    }
    else {
        // Get the bounding temperature and σ₀ values
        t0 = fabs(tem[interv]);
        t1 = fabs(tem[interv + 1]);
        s0 = sig[interv];
        s1 = sig[interv + 1];

        // Compute logarithmic interpolation for σ₀
        *sigm0 = pow(10.0, log10(s0) +
                          (log10(s1) - log10(s0)) / (log10(t1) - log10(t0)) *
                          (log10(fabs(temp)) - log10(t0)));
    }

    return;
}

/* ==========================================================================
 * Three-phase (ice / sediment / air) material properties and triple well.
 * Used only by the 3-phase model (dof=4, assembly_sed.c). air = 1 - ice - sed.
 * Callers pass clamped fractions; these are exact linear mixtures with
 * constant analytic derivatives (d(air)/d(ice) = d(air)/d(sed) = -1).
 * ========================================================================== */

/* Effective thermal conductivity K(phi) = ice*k_i + sed*k_s + air*k_a. */
void ThermalCond3(AppCtx *user, PetscScalar ice, PetscScalar sed,
                  PetscScalar *cond, PetscScalar *dcond_ice, PetscScalar *dcond_sed)
{
    PetscReal ki = user->thcond_ice, ka = user->thcond_air, ks = user->thcond_sed;
    PetscScalar air = 1.0 - ice - sed;
    if (cond)      (*cond)      = ice*ki + sed*ks + air*ka;
    if (dcond_ice) (*dcond_ice) = ki - ka;
    if (dcond_sed) (*dcond_sed) = ks - ka;
}

/* Effective heat capacity cp(phi) = ice*cp_i + sed*cp_s + air*cp_a. */
void HeatCap3(AppCtx *user, PetscScalar ice, PetscScalar sed,
              PetscScalar *cp, PetscScalar *dcp_ice, PetscScalar *dcp_sed)
{
    PetscReal ci = user->cp_ice, ca = user->cp_air, cs = user->cp_sed;
    PetscScalar air = 1.0 - ice - sed;
    if (cp)      (*cp)      = ice*ci + sed*cs + air*ca;
    if (dcp_ice) (*dcp_ice) = ci - ca;
    if (dcp_sed) (*dcp_sed) = cs - ca;
}

/* Effective density rho(phi) = ice*rho_i + sed*rho_s + air*rho_a. */
void Density3(AppCtx *user, PetscScalar ice, PetscScalar sed,
              PetscScalar *rho, PetscScalar *drho_ice, PetscScalar *drho_sed)
{
    PetscReal ri = user->rho_ice, ra = user->rho_air, rs = user->rho_sed;
    PetscScalar air = 1.0 - ice - sed;
    if (rho)      (*rho)      = ice*ri + sed*rs + air*ra;
    if (drho_ice) (*drho_ice) = ri - ra;
    if (drho_sed) (*drho_sed) = rs - ra;
}

/* Triple-well driving force for the ice equation:
 *   g = dF^tri/dphi_i - dF^tri/dphi_a   (the beta-eliminated difference)
 * with F^tri = 1/2 * sum_k Sigma_k phi_k^2 (1-phi_k)^2 + Lambda phi_i^2 phi_s^2 phi_a^2,
 * phi_a = 1 - phi_i - phi_s. Returns g and its derivatives w.r.t. the two
 * independent DOFs (phi_i, phi_s). Reduces to the 2-phase double well
 * (Sigma_i+Sigma_a) phi_i(1-phi_i)(1-2phi_i) when sed = 0. See the
 * derivation in studies/icy_regolith/explicit_sediment_phase/PLAN.md.
 *
 *   w(p)  = p(1-p)(1-2p),   w'(p) = 1 - 6p + 6p^2
 *   g        = Si w(i) - Sa w(a) + 2 L s^2 ( i a^2 - i^2 a )
 *   dg/di    = Si w'(i) + Sa w'(a) + 2 L s^2 ( a^2 - 4 i a + i^2 )
 *   dg/ds    = Sa w'(a) + 2 L [ 2 s i a (a - i) + s^2 ( i^2 - 2 i a ) ]
 */
void TripleWell(AppCtx *user, PetscScalar ice, PetscScalar sed,
                PetscScalar *g, PetscScalar *dg_ice, PetscScalar *dg_sed)
{
    PetscReal Si = user->Sigma_i, Sa = user->Sigma_a, L = user->Lambd;
    PetscScalar i = ice, s = sed, a = 1.0 - ice - sed;
#define WELL(p)  ((p)*(1.0-(p))*(1.0-2.0*(p)))
#define DWELL(p) (1.0 - 6.0*(p) + 6.0*(p)*(p))
    if (g)
        (*g)      = Si*WELL(i) - Sa*WELL(a) + 2.0*L*s*s*(i*a*a - i*i*a);
    if (dg_ice)
        (*dg_ice) = Si*DWELL(i) + Sa*DWELL(a) + 2.0*L*s*s*(a*a - 4.0*i*a + i*i);
    if (dg_sed)
        (*dg_sed) = Sa*DWELL(a)
                  + 2.0*L*( 2.0*s*i*a*(a - i) + s*s*(i*i - 2.0*i*a) );
#undef WELL
#undef DWELL
}
