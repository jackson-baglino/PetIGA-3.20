#include "material_properties.h"
#include "NASA_types.h"

/**
 * @brief Computes the effective thermal conductivity based on phase fractions.
 * 
 * This function calculates the thermal conductivity of a material as a weighted 
 * sum of the thermal conductivities of ice, sediment, and air phases. It also computes 
 * the derivative of conductivity with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase (0 to 1).
 * @param sed Fraction of the sediment phase (0 to 1).
 * @param cond Pointer to store the computed thermal conductivity.
 * @param dcond_ice Pointer to store the derivative of conductivity with respect to ice fraction.
 */
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar sed, 
                 PetscScalar *cond, PetscScalar *dcond_ice)
{
    // Define derivatives (initialized to 1.0, modified if phase fractions are negative)
    PetscReal dice = 1.0, dair = 1.0;

    // Compute air fraction (air = 1 - ice - sediment)
    PetscReal air = 1.0 - sed - ice;

    // Ensure phase fractions are non-negative (corrects numerical issues)
    if (sed < 0.0) sed = 0.0;
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
    PetscReal cond_sed = user->thcond_sed;
    PetscReal cond_air = user->thcond_air;

    // Compute the effective thermal conductivity as a weighted sum
    if (cond) 
        (*cond) = ice * cond_ice + sed * cond_sed + air * cond_air;

    // Compute the derivative of thermal conductivity with respect to ice fraction
    if (dcond_ice) 
        (*dcond_ice) = cond_ice * dice - cond_air * dair;

    return; // Explicit return statement for clarity
}

/**
 * @brief Computes the effective heat capacity and its derivative with respect to ice.
 * 
 * This function calculates the heat capacity as a weighted sum of contributions 
 * from ice, sediment, and air. It also computes the derivative of heat capacity 
 * with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param sed Fraction of the sediment phase.
 * @param cp Pointer to store computed heat capacity.
 * @param dcp_ice Pointer to store derivative of heat capacity with respect to ice.
 */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *cp, 
  PetscScalar *dcp_ice)
{
PetscReal dice = 1.0, dair = 1.0;
PetscReal air = 1.0 - sed - ice;

// Ensure phase fractions are non-negative
if (sed < 0.0) sed = 0.0;
if (ice < 0.0) { ice = 0.0; dice = 0.0; }
if (air < 0.0) { air = 0.0; dair = 0.0; }

// Retrieve heat capacities
PetscReal cp_ice = user->cp_ice;
PetscReal cp_sed = user->cp_sed;
PetscReal cp_air = user->cp_air;

// Compute effective heat capacity
if (cp) 
(*cp) = ice * cp_ice + sed * cp_sed + air * cp_air;

// Compute derivative with respect to ice
if (dcp_ice) 
(*dcp_ice) = cp_ice * dice - cp_air * dair;
return;
}

/**
 * @brief Computes the effective density and its derivative with respect to ice.
 *
 * This function calculates the density as a weighted sum of contributions 
 * from ice, sediment, and air. It also computes the derivative of density 
 * with respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param sed Fraction of the sediment phase.
 * @param rho Pointer to store computed density.
 * @param drho_ice Pointer to store derivative of density with respect to ice.
 */
void Density(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *rho, 
  PetscScalar *drho_ice)
{
PetscReal dice = 1.0, dair = 1.0;
PetscReal air = 1.0 - sed - ice;

// Ensure phase fractions are non-negative
if (sed < 0.0) sed = 0.0;
if (ice < 0.0) { ice = 0.0; dice = 0.0; }
if (air < 0.0) { air = 0.0; dair = 0.0; }

// Retrieve densities
PetscReal rho_ice = user->rho_ice;
PetscReal rho_sed = user->rho_sed;
PetscReal rho_air = user->rho_air;

// Compute effective density
if (rho) 
(*rho) = ice * rho_ice + sed * rho_sed + air * rho_air;

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
 * @brief Computes a smooth approximation of the Heaviside step function.
 *
 * This function evaluates a smooth polynomial approximation of the Heaviside
 * step function, commonly used in phase-field and level-set methods to provide
 * a differentiable transition instead of a sharp discontinuity.
 *
 * The smooth Heaviside function is defined as:
 * H(φ) = φ³(3 - 2φ)
 *
 * This provides a smooth transition from 0 to 1 over the interval [0, 1].
 *
 * @param[in]  phi  Input scalar value (typically a phase field or level-set value)
 * @param[out] g    Pointer to output scalar where the result is stored.
 *                  If NULL, computation is skipped.
 * @param[out] dg_dphi Pointer to output scalar where the derivative of g with respect to phi is 
 *                     stored. If NULL, computation is skipped.
 *
 * @return void
 *
 * @note The function assumes phi is bounded within [0, 1] for meaningful results.
 *       Values outside this range will produce results outside [0, 1].
 */
void SmoothHeavisidePoly(PetscScalar phi, PetscScalar *g, PetscScalar *dg_dphi)
{
    if (g)
        (*g) = phi*phi*phi * (3 - 2*phi); // Smooth approximation of Heaviside function
    if (dg_dphi)
        (*dg_dphi) = 6*phi*phi * (1 - phi); // Derivative of the smooth Heaviside function with respect to phi
    return;
}

/* Penalty weight: SmoothHeavisidePoly applied to a shifted variable so the
 * penalty is concentrated near phi = 1 (deep solid) rather than spread across
 * the diffuse interface. Zero for phi <= PENALTY_PHI_LO, unity for
 * phi >= PENALTY_PHI_HI, smooth in between. Width is small enough that the
 * mid-thickness of a typical ice-air interface (ice+sed ≈ 0.5) falls fully in
 * the zero region, so Gibbs-Thomson can emerge at the interface, while the
 * solid interior (ice+sed → 1) still has its rhov pinned to rhov_eq. */
#define PENALTY_PHI_LO 0.90
#define PENALTY_PHI_HI 1.00
void PenaltyWeight(PetscScalar phi, PetscScalar *g, PetscScalar *dg_dphi)
{
    const PetscReal lo = PENALTY_PHI_LO;
    const PetscReal hi = PENALTY_PHI_HI;
    if (phi <= lo) {
        if (g)       (*g)       = 0.0;
        if (dg_dphi) (*dg_dphi) = 0.0;
        return;
    }
    if (phi >= hi) {
        if (g)       (*g)       = 1.0;
        if (dg_dphi) (*dg_dphi) = 0.0;
        return;
    }
    const PetscReal inv_w = 1.0 / (hi - lo);
    PetscScalar u = (phi - lo) * inv_w;
    PetscScalar g_u, dg_du;
    SmoothHeavisidePoly(u, &g_u, &dg_du);
    if (g)       (*g)       = g_u;
    if (dg_dphi) (*dg_dphi) = dg_du * inv_w;
}

/* Curvature of an isosurface of phi from its gradient and Hessian. See header
 * for full description. */
void Curvature(PetscInt dim,
               const PetscScalar grad_phi[],
               const PetscScalar hess_phi[],
               PetscReal eps_reg,
               PetscScalar *kappa,
               PetscScalar dkappa_dg[],
               PetscScalar dkappa_dH[])
{
    /* 1D: curvature is identically zero physically; the regularized formula
     * would otherwise produce a small artifact. Return zeros and exit. */
    if (dim == 1) {
        if (kappa) (*kappa) = 0.0;
        if (dkappa_dg) dkappa_dg[0] = 0.0;
        if (dkappa_dH) dkappa_dH[0] = 0.0;
        return;
    }

    PetscInt k, l;

    /* G^2 = |g|^2 + eps_reg^2  (Tikhonov-style regularization)
     * Keeps kappa bounded where |grad_phi| -> 0 in bulk regions. */
    PetscScalar G2 = (PetscScalar)(eps_reg * eps_reg);
    for (l = 0; l < dim; l++) G2 += grad_phi[l] * grad_phi[l];
    PetscScalar G  = PetscSqrtScalar(G2);
    PetscScalar G3 = G2 * G;
    PetscScalar G5 = G3 * G2;

    /* Laplacian L = trace(H) */
    PetscScalar L = 0.0;
    for (k = 0; k < dim; k++) L += hess_phi[k * dim + k];

    /* g . H . g */
    PetscScalar gHg = 0.0;
    for (k = 0; k < dim; k++)
        for (l = 0; l < dim; l++)
            gHg += grad_phi[k] * hess_phi[k * dim + l] * grad_phi[l];

    if (kappa) {
        (*kappa) = -L / G + gHg / G3;
    }

    if (dkappa_dg) {
        /* (H g)_l = sum_k H_{lk} g_k */
        PetscScalar Hg[3] = {0.0, 0.0, 0.0};
        for (l = 0; l < dim; l++) {
            PetscScalar s = 0.0;
            for (k = 0; k < dim; k++) s += hess_phi[l * dim + k] * grad_phi[k];
            Hg[l] = s;
        }
        for (l = 0; l < dim; l++) {
            dkappa_dg[l] = (L * grad_phi[l] + 2.0 * Hg[l]) / G3
                         - 3.0 * gHg * grad_phi[l] / G5;
        }
    }

    if (dkappa_dH) {
        /* d kappa / d H_{kl} = -delta_{kl} / G + g_k g_l / G^3 */
        for (k = 0; k < dim; k++)
            for (l = 0; l < dim; l++) {
                PetscScalar delta = (k == l) ? 1.0 : 0.0;
                dkappa_dH[k * dim + l] = -delta / G
                                       + grad_phi[k] * grad_phi[l] / G3;
            }
    }
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

/**
 * @brief Computes the free energy function for the ice phase.
 * 
 * This function calculates the free energy contribution of the ice phase in the 
 * phase-field model. It also computes the derivative of this energy with 
 * respect to the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param sed Fraction of the sediment phase.
 * @param fice Pointer to store the computed free energy for the ice phase.
 * @param dfice_ice Pointer to store the derivative of fice with respect to ice.
 */
void Fice(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fice, 
  PetscScalar *dfice_ice)
{
// Retrieve material parameters from user context
PetscReal Lambd = user->Lambd;
PetscReal etai  = user->Etai;

// Compute the air fraction
PetscReal air = 1.0 - sed - ice;

// Compute the free energy function for the ice phase
if (fice) 
(*fice) = etai * ice * (1.0 - ice) * (1.0 - 2.0 * ice) + 
          2.0 * Lambd * ice * sed * sed * air * air;

// Compute the derivative with respect to the ice fraction
if (dfice_ice) 
(*dfice_ice) = etai * (1.0 - 6.0 * ice + 6.0 * ice * ice) + 
              2.0 * Lambd * sed * sed * air * air - 
              2.0 * Lambd * ice * sed * sed * 2.0 * air;

return;
}

/**
 * @brief Computes the phase evolution function for the water phase (sediment fraction).
 * 
 * This function models the free energy contribution of the water (sediment) phase 
 * in a phase-field simulation. It also computes the derivative with respect to 
 * the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param sed Fraction of the sediment phase.
 * @param fsed Pointer to store the computed phase function for water.
 * @param dfsed_ice Pointer to store the derivative of fsed with respect to ice.
 */
void Fsed(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fsed, 
          PetscScalar *dfsed_ice)
{
    // Retrieve material parameters from user context
    PetscReal Lambd = user->Lambd;
    PetscReal etam  = user->Etam;

    // Compute the air fraction
    PetscReal air = 1.0 - sed - ice;

    // Compute the phase evolution function for the water phase
    if (fsed) 
        (*fsed) = etam * sed * (1.0 - sed) * (1.0 - 2.0 * sed) + 
                  2.0 * Lambd * ice * ice * sed * air * air;

    // Compute the derivative with respect to the ice fraction
    if (dfsed_ice) {
        (*dfsed_ice) = 2.0 * Lambd * 2.0 * ice * sed * air * air 
                     - 2.0 * Lambd * ice * ice * sed * 2.0 * air;
    }

    return;
}

/**
 * @brief Computes the phase evolution function for the air phase.
 * 
 * This function models the free energy contribution of the air phase in a 
 * phase-field simulation. It also computes the derivative with respect to 
 * the ice fraction.
 *
 * @param user Pointer to the application context containing material properties.
 * @param ice Fraction of the ice phase.
 * @param sed Fraction of the sediment phase.
 * @param fair Pointer to store the computed phase function for air.
 * @param dfair_ice Pointer to store the derivative of fair with respect to ice.
 */
void Fair(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fair, 
          PetscScalar *dfair_ice)
{
    // Retrieve material parameters from user context
    PetscReal Lambd = user->Lambd;
    PetscReal etaa  = user->Etaa;

    // Compute the air fraction
    PetscReal air = 1.0 - sed - ice;

    // Compute the phase evolution function for the air phase
    if (fair) 
        (*fair) = etaa * air * (1.0 - air) * (1.0 - 2.0 * air) + 
                  2.0 * Lambd * ice * ice * sed * sed * air;

    // Compute the derivative with respect to the ice fraction
    if (dfair_ice) {
        (*dfair_ice)  = -etaa * (1.0 - air) * (1.0 - 2.0 * air) 
                      + etaa * air * (1.0 - 2.0 * air) 
                      + etaa * air * (1.0 - air) * 2.0;
        (*dfair_ice) += 2.0 * Lambd * 2.0 * ice * sed * sed * air 
                      - 2.0 * Lambd * ice * ice * sed * sed;
    }

    return;
}

void Mobility(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *mob)
{
    /* Constant, phase-independent mobility.
     *
     * The previous phase-weighted form
     *      mob = mob_sub*ice + mob_sed*sed + mob_air*air
     * gave a 0.5x slowdown at ice-sed interfaces. With mob_sed=0 (sediment
     * inert by default) and mob_air=mob_sub, this collapsed to mob_sub*(1-sed)
     * — i.e. full mobility in air/ice regions but half at ice-sed contacts.
     *
     * The slowdown was originally protective: it prevented AC from dissolving
     * the sed interior during a 3-phase pre-pin window (t < t_sed_freeze).
     * With t_sed_freeze=0 (the new default), sed is pinned by R_sed = sed_t = 0
     * from t=0 and never has a 3-phase AC window, so the slowdown serves no
     * purpose — it only unphysically reduced the ice AC rate at ice-sed
     * contacts (slowing sintering near sediment grains).
     *
     * Now mob_sub everywhere. If t_sed_freeze > 0 is ever reintroduced for a
     * legitimate AC pre-pin period, restore the phase-weighted form here.
     */
    (void)ice;  /* unused — kept in signature for callers */
    (void)sed;  /* unused */
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