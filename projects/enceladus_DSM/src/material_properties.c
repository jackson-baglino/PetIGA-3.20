#include "material_properties.h"

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


/**
 * @brief Condensation coefficient alpha_c(T, rho_v) and its derivatives.
 *
 * alpha_c is the model's one genuinely free parameter, and it should be a
 * FUNCTION OF THE LOCAL STATE rather than a number chosen per run: it is set
 * by the temperature and by the local supersaturation, both of which the
 * solver already carries at every quadrature point. Fixing it as a scalar
 * turns a determined quantity into a knob.
 *
 * Three models, selected by user->alpha_model:
 *
 *   ALPHA_MODEL_CONST (0)  alpha_c = user->alpha_c0.  Reproduces the previous
 *                          behaviour exactly and is the default, so nothing
 *                          changes until a run opts in.
 *
 *   ALPHA_MODEL_ARRH (1)   alpha_c = A*exp(-Ea/(R*T)).  Temperature only.
 *                          Anchored on the alpha_c range Braun et al. (2024)
 *                          report; note THEY give no temperature dependence
 *                          (they hold alpha constant and say so), so the
 *                          T-dependence here is our assumption, not theirs.
 *
 *   ALPHA_MODEL_LIBB2 (2)  alpha_c = A*exp(-f*sigma0(T)/sigma), the Libbrecht
 *                          form with BOTH free parameters fitted rather than
 *                          A hardwired to 1. sigma = |rho_v - rho_vs|/rho_vs
 *                          is the local supersaturation, so this is the full
 *                          alpha_c(T, rho_v).
 *
 * CLAMPED to [user->alpha_lo, user->alpha_hi] in every model. The floor is
 * not cosmetic: alpha_c -> 0 as sigma -> 0 sends beta_sub ~ 1/alpha_c to
 * infinity, which is both unphysical and unsolvable. The clamp is applied
 * BEFORE the derivatives, and the derivatives are zeroed where it binds --
 * otherwise the Jacobian would advertise a sensitivity the residual does not
 * have, and Newton would chase it.
 *
 * Mirrors preprocess/comp_eps.py (alpha_arrhenius / alpha_libbrecht2). Keep
 * the two in step: comp_eps sizes the mesh from the same alpha_c the solver
 * evaluates, and a divergence between them silently invalidates the mesh.
 */
void AlphaCondensation(AppCtx *user, PetscScalar tem, PetscScalar rhov,
                       PetscScalar *alpha, PetscScalar *dalpha_dtem,
                       PetscScalar *dalpha_drhov)
{
    PetscReal a = user->alpha_c0, da_dT = 0.0, da_drv = 0.0;
    const PetscReal T_K = PetscRealPart(tem) + 273.15;

    if (user->alpha_model == ALPHA_MODEL_ARRH) {
        /* a = A exp(-Ea/(R T));  da/dT = a * Ea/(R T^2). */
        a     = user->alpha_A * exp(-user->alpha_Ea / (ALPHA_R_GAS * T_K));
        da_dT = a * user->alpha_Ea / (ALPHA_R_GAS * T_K * T_K);
    }
    else if (user->alpha_model == ALPHA_MODEL_LIBB2) {
        PetscScalar rho_vs, d_rho_vs, s0;
        RhoVS_I(user, tem, &rho_vs, &d_rho_vs);
        Sigma0(tem, &s0);

        const PetscReal rvs   = PetscRealPart(rho_vs);
        const PetscReal drvs  = PetscRealPart(d_rho_vs);
        const PetscReal dev   = PetscRealPart(rhov) - rvs;      /* signed */
        const PetscReal sgn   = (dev >= 0.0) ? 1.0 : -1.0;
        const PetscReal sig   = fabs(dev) / rvs;

        /* Below this the exponent overflows and alpha underflows to zero; the
         * clamp would catch it anyway, but evaluating exp(-1e30) first is a
         * needless way to generate a denormal. */
        if (sig <= 1.0e-30) {
            a = user->alpha_lo;
        } else {
            const PetscReal expo = user->alpha_f * PetscRealPart(s0) / sig;
            if (expo > 700.0) {
                a = 0.0;
            } else {
                a = user->alpha_A * exp(-expo);
                /* d sigma/d rho_v = sgn/rvs
                 * d sigma/d T     = -sgn*dev*drvs/rvs^2 - sgn*drvs/rvs ... via
                 *   sigma = |rho_v - rvs(T)|/rvs(T):
                 *   dsig/dT = [-sgn*drvs*rvs - |dev|*drvs] / rvs^2
                 *           = -drvs*(sgn*rvs + fabs(dev))/rvs^2
                 * and da/dx = a * (f*s0/sigma^2) * dsigma/dx  (s0's own T
                 * dependence is a table interpolant; its derivative is
                 * piecewise-constant and is deliberately omitted -- see NOTE). */
                const PetscReal ds_drv = sgn / rvs;
                const PetscReal ds_dT  = -drvs * (sgn * rvs + fabs(dev)) / (rvs * rvs);
                const PetscReal pref   = a * user->alpha_f * PetscRealPart(s0) / (sig * sig);
                da_drv = pref * ds_drv;
                da_dT  = pref * ds_dT;
            }
        }
        /* NOTE: d(sigma0)/dT is dropped. Sigma0() is a 10-point log-log table
         * interpolant with a genuine kink at -6/-7 C, so its derivative is
         * discontinuous and, at the kink, meaningless. Including it would make
         * the Jacobian inconsistent with the residual exactly where the table
         * misbehaves. The omission makes the Jacobian slightly inexact in T;
         * check -snes_test_jacobian after enabling this model. */
    }

    /* Clamp, and kill the derivatives where the clamp binds. */
    if (a <= user->alpha_lo) { a = user->alpha_lo; da_dT = 0.0; da_drv = 0.0; }
    if (a >= user->alpha_hi) { a = user->alpha_hi; da_dT = 0.0; da_drv = 0.0; }

    if (alpha)        (*alpha)        = a;
    if (dalpha_dtem)  (*dalpha_dtem)  = da_dT;
    if (dalpha_drhov) (*dalpha_drhov) = da_drv;
}


/**
 * @brief Pointwise sublimation kinetics: mob_sub and alph_sub from alpha_c(T, rho_v).
 *
 * main() computes these once as scalars from a single alpha_c. That is only
 * right if alpha_c is a constant. Once alpha_c depends on the local state,
 * everything downstream of it does too, via the M&F SI Eq. 9 chain:
 *
 *   chi        = rho_vs(T)/rho_ice
 *   d0_sub     = d0_sub0 * chi
 *   lambda_sub = a1*eps/d0_sub
 *   beta_s     = (1/alpha_c) * sqrt(2*pi*m/(kB*T))      [scaled, K&P beta']
 *   B          = beta_s/a1 + a2*eps/diff_sub + a2*eps/D_v(T)
 *   tau_sub    = eps*lambda_sub*B
 *   mob_sub    = eps/(3*tau_sub) = 1/(3*lambda_sub*B)
 *   alph_sub   = lambda_sub/tau_sub = 1/(eps*B)
 *
 * Note alph_sub depends on T and rho_v ONLY through B, while mob_sub also
 * picks up lambda_sub's temperature dependence through d0_sub.
 *
 * With -alpha_model 0 (the default) alpha_c is constant, and the only residual
 * state dependence is D_v(T) and rho_vs(T). So even in CONST mode this is a
 * slight refinement, not a no-op. To fall back to main()'s scalars exactly,
 * set -alpha_pointwise 0.
 *
 * Both outputs are then multiplied by user->mob_scale / user->alph_scale
 * (-mob_scale / -alph_scale, 1.0 by default), which main() applies to its own
 * scalars as well -- so those two knobs bite on either path, unlike the
 * absolute -mob_sub / -alph_sub overrides, which this function ignores.
 *
 * Derivatives are analytic and exact for the CONST and ARRH models. For LIBB2
 * they inherit the dropped d(sigma0)/dT term from AlphaCondensation (~1% in
 * dalpha/dT) -- check -snes_test_jacobian before trusting that model.
 */
void SubKinetics(AppCtx *user, PetscScalar tem, PetscScalar rhov,
                 PetscScalar *mob,  PetscScalar *alph,
                 PetscScalar *dmob_dT,  PetscScalar *dmob_drv,
                 PetscScalar *dalph_dT, PetscScalar *dalph_drv)
{
    const PetscReal a1 = 5.0, a2 = 0.1581;      /* M&F SI footnote 1 */
    const PetscReal kB = K_BOLTZ, m_h2o = M_H2O_KG;
    const PetscReal eps = user->eps;
    const PetscReal T_K = PetscRealPart(tem) + 273.15;

    PetscScalar alpha, da_dT, da_drv;
    AlphaCondensation(user, tem, rhov, &alpha, &da_dT, &da_drv);

    PetscScalar rho_vs, d_rho_vs, dv, d_dv;
    RhoVS_I(user, tem, &rho_vs, &d_rho_vs);
    VaporDiffus(user, tem, &dv, &d_dv);

    const PetscReal rvs = PetscRealPart(rho_vs), drvs = PetscRealPart(d_rho_vs);
    const PetscReal Dv = PetscRealPart(dv),      dDv = PetscRealPart(d_dv);

    /* lambda_sub = a1*eps/(d0_sub0*chi),  chi = rho_vs/rho_ice */
    const PetscReal chi   = rvs / user->rho_ice;
    const PetscReal dchi  = drvs / user->rho_ice;
    const PetscReal d0s   = user->d0_sub0 * chi;
    const PetscReal lam   = a1 * eps / d0s;
    const PetscReal dlam  = -lam * dchi / chi;          /* d lambda / dT */

    /* beta_s = C(T)/alpha,  C = sqrt(2*pi*m/(kB*T)),  dC/dT = -C/(2T) */
    const PetscReal Cth   = sqrt(2.0 * PETSC_PI * m_h2o / (kB * T_K));
    const PetscReal beta_s = Cth / PetscRealPart(alpha);
    const PetscReal dbeta_dT  = (-Cth / (2.0 * T_K)) / PetscRealPart(alpha)
                              - beta_s * PetscRealPart(da_dT)  / PetscRealPart(alpha);
    const PetscReal dbeta_drv = - beta_s * PetscRealPart(da_drv) / PetscRealPart(alpha);

    /* B and its derivatives (diff_sub is a constant; D_v carries T) */
    const PetscReal B    = beta_s / a1 + a2 * eps / user->diff_sub + a2 * eps / Dv;
    const PetscReal dB_dT  = dbeta_dT  / a1 - a2 * eps * dDv / (Dv * Dv);
    const PetscReal dB_drv = dbeta_drv / a1;

    if (alph)      (*alph)      = 1.0 / (eps * B);
    if (dalph_dT)  (*dalph_dT)  = -dB_dT  / (eps * B * B);
    if (dalph_drv) (*dalph_drv) = -dB_drv / (eps * B * B);

    if (mob)       (*mob)       = 1.0 / (3.0 * lam * B);
    if (dmob_dT)   (*dmob_dT)   = -(dlam * B + lam * dB_dT) / (3.0 * lam * lam * B * B);
    if (dmob_drv)  (*dmob_drv)  = -dB_drv / (3.0 * lam * B * B);

    /* Empirical multipliers (-mob_scale / -alph_scale, both 1.0 by default;
     * main() applies the same two to its scalars). The scaling is linear, so
     * each derivative carries its own factor and the analytic Jacobian stays
     * exact -- verify with -snes_test_jacobian after touching this. */
    if (alph)      (*alph)      *= user->alph_scale;
    if (dalph_dT)  (*dalph_dT)  *= user->alph_scale;
    if (dalph_drv) (*dalph_drv) *= user->alph_scale;

    if (mob)       (*mob)       *= user->mob_scale;
    if (dmob_dT)   (*dmob_dT)   *= user->mob_scale;
    if (dmob_drv)  (*dmob_drv)  *= user->mob_scale;
}
