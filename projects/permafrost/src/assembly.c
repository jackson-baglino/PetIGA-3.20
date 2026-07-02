#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * 2-phase (ice / temperature / vapor) system. phi_a = 1 - phi_i is
 * computed algebraically; there is no sediment DOF.
 *
 * Free energy: single double-well in phi_i,
 *   f1(phi_i) = phi_i*(1-phi_i)*(1-2*phi_i)
 *
 * Ice (Allen-Cahn — pure curvature relaxation; sublimation source off):
 *   R_ice = N0[a]*ice_t + 3*mob_sub*eps*grad_N.grad_ice
 *         + (3*mob_sub/eps)*f1(phi_i)*N0[a]
 *
 * This is the met=0 reduction of dry_snow_metamorphism's 3-phase ice
 * equation: R_ice += mob*3/eps/ETA * ((Etam+Etaa)*fice - Etaa*fmet - Etam*fair).
 * With met=0, ETA=Etaa*Etai, fmet=0, fice=Etai*f1(ice), so the expression
 * collapses to mob*3/eps*f1(ice) -- Sigma_i/Sigma_a/Lambda all cancel and
 * play no role in a true 2-phase (no third phase) system.
 *
 * Temperature (row-scaled by S_T = 1/(rho_ice*lat_sub), numerical
 * preconditioning only — no latent-heat source while sublimation is off):
 *   R_tem = rho*cp*N0[a]*tem_t + thcond*grad_N.grad_tem
 *
 * Vapor (air_lim floor as in dry_snow_metamorphism):
 *   R_vap = N0[a]*air_eff*rhov_t + difvap*air_eff*grad_N.grad_rhov
 *   air_eff = (air > air_lim) ? air : air_lim,   air = 1 - ice
 * ========================================================================= */

/* Double-well derivative f1(phi_i) = phi_i(1-phi_i)(1-2phi_i) and its
 * derivative df1/dphi_i = 1 - 6 phi_i + 6 phi_i^2. */
static void DoubleWellDeriv(PetscReal ice, PetscReal *f1, PetscReal *df1)
{
    if (f1)  *f1  = ice * (1.0 - ice) * (1.0 - 2.0 * ice);
    if (df1) *df1 = 1.0 - 6.0 * ice + 6.0 * ice * ice;
}


/* =========================================================================
 * Residual
 * ========================================================================= */
PetscErrorCode Residual_A1(IGAPoint pnt,
                           PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U,
                           PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    /* Quick-and-dirty: model eps < IC eps (75%) to counter observed
     * post-IC interface growth/coarsening. TODO: split into a proper
     * user->eps_model field if this sticks. */
    PetscReal eps     = 0.75 * user->eps;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice   = sol[0],  ice_t  = sol_t[0];
    PetscScalar tem   = sol[1],  tem_t  = sol_t[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    (void)rhov; /* used only in the commented-out sublimation source below */
    PetscScalar grad_ice [dim];
    PetscScalar grad_tem [dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_ice [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar air = 1.0 - ice;

    /* Sublimation/GT terms are commented out; hessian no longer needed.
    PetscScalar hess_sol[3][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);
    PetscScalar hess_ice[dim * dim];
    for (l = 0; l < dim * dim; l++) hess_ice[l] = hess_sol[0][l / dim][l % dim];
    */

    /* SNES domain-error catch — if a trial Newton iterate has phi out of the
     * configured bounds, signal an invalid state so the line search backs off
     * the step. Cheaper than committing a bad state and rolling back later. */
    {
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (PetscRealPart(ice) < lo || PetscRealPart(ice) > hi ||
            PetscRealPart(air) < lo || PetscRealPart(air) > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    /* Clamp clipping copies of phi for material-property evaluation only —
     * keeps Density/HeatCap/etc. from getting negative inputs near the bounds. */
    PetscScalar ice_c = ice, air_c = air;
    if (PetscRealPart(ice_c) < 0.0) ice_c = 0.0;
    if (PetscRealPart(ice_c) > 1.0) ice_c = 1.0;
    if (PetscRealPart(air_c) < 0.0) air_c = 0.0;
    if (PetscRealPart(air_c) > 1.0) air_c = 1.0;

    /* Material properties */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, ice_c, &thcond,  NULL);
    HeatCap    (user, ice_c, &cp,      NULL);
    Density    (user, ice_c, &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, ice_c, &mob_sub);

    /* Gibbs-Thomson curvature correction (sublimation off; commented out).
    PetscScalar kappa = 0.0;
    if (user->d0_GT != 0.0)
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps, &kappa, NULL, NULL);
    PetscReal rhoI_vs_eff = rho_vs * (1.0 + user->d0_GT * PetscRealPart(kappa));
    */

    PetscReal f1;
    DoubleWellDeriv(PetscRealPart(ice_c), &f1, NULL);

    /* air_eff = max(air, air_lim) */
    PetscReal air_r = PetscRealPart(air_c);
    PetscReal air_eff = (air_r > air_lim) ? air_r : air_lim;

    /* Phase-change source (sublimation off; commented out).
    PetscReal loc   = PetscRealPart(ice_c) * PetscRealPart(ice_c)
                    * air_r * air_r;
    PetscReal S_sub = alph_sub * loc * (PetscRealPart(rhov) - rhoI_vs_eff) / rho_ice;
    */

    /* T-equation row-scale (numerical preconditioning; no physics change). */
    PetscReal S_T = 1.0 / (rho_ice * lat_sub);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar (*R)[3] = (PetscScalar (*)[3])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal grad_N_dot_grad_ice  = 0.0;
        PetscReal grad_N_dot_grad_tem  = 0.0;
        PetscReal grad_N_dot_grad_rhov = 0.0;
        for (l = 0; l < dim; l++) {
            grad_N_dot_grad_ice  += N1[a][l] * PetscRealPart(grad_ice [l]);
            grad_N_dot_grad_tem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            grad_N_dot_grad_rhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }

        PetscScalar R_ice = N0[a] * ice_t
                          + 3.0 * mob_sub * eps * grad_N_dot_grad_ice
                          + (3.0 * mob_sub / eps) * f1 * N0[a];

        PetscScalar R_tem = rho * cp * N0[a] * tem_t
                          + thcond * grad_N_dot_grad_tem;

        PetscScalar R_vap = N0[a] * air_eff * rhov_t
                          + dif_vap * air_eff * grad_N_dot_grad_rhov;

        R[a][0] = R_ice;
        R[a][1] = S_T * R_tem;
        R[a][2] = R_vap;
    }
    return 0;
}


/* =========================================================================
 * Residual dispatcher (single avenue)
 * ========================================================================= */
PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    return Residual_A1(pnt, shift, V, t, U, Re, ctx);
}


/* =========================================================================
 * Jacobian.  J[a][i][b][j] = d R[a][i] / d u[b][j] + shift * d R[a][i] / d u_t[b][j]
 *
 * DOF layout: 0=ice, 1=tem, 2=rhov.
 * ========================================================================= */
static PetscErrorCode Jacobian_A1(IGAPoint pnt,
                                  PetscReal shift, const PetscScalar *V,
                                  PetscReal t, const PetscScalar *U,
                                  PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    /* Quick-and-dirty: model eps < IC eps (75%), matches Residual_A1. */
    PetscReal eps     = 0.75 * user->eps;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice    = sol[0];
    PetscScalar tem    = sol[1];
    PetscScalar rhov   = sol[2],  rhov_t = sol_t[2];
    (void)rhov; /* used only in the commented-out sublimation Jacobian terms below */
    PetscScalar grad_ice[dim];
    PetscScalar grad_tem[dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_ice [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar air = 1.0 - ice;

    /* Hessian of ice for GT-curvature derivatives (sublimation off; commented out).
    PetscScalar hess_sol[3][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);
    PetscScalar hess_ice[dim * dim];
    for (l = 0; l < dim * dim; l++) hess_ice[l] = hess_sol[0][l / dim][l % dim];
    */

    /* Clamp copies for material-property evaluation */
    PetscScalar ice_c = ice, air_c = air;
    if (PetscRealPart(ice_c) < 0.0) ice_c = 0.0;
    if (PetscRealPart(ice_c) > 1.0) ice_c = 1.0;
    if (PetscRealPart(air_c) < 0.0) air_c = 0.0;
    if (PetscRealPart(air_c) > 1.0) air_c = 1.0;

    /* Material properties + derivatives */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, ice_c, &thcond,  NULL);
    HeatCap    (user, ice_c, &cp,      NULL);
    Density    (user, ice_c, &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, ice_c, &mob_sub);

    PetscReal dthcond_dice;
    ThermalCond(user, ice_c, NULL, &dthcond_dice);

    PetscReal d_dif_vap;
    VaporDiffus(user, tem, NULL, &d_dif_vap);

    /* Gibbs-Thomson curvature correction (sublimation off; commented out).
    PetscScalar kappa = 0.0;
    PetscScalar dkappa_dg[3] = {0.0, 0.0, 0.0};
    PetscScalar dkappa_dH[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (user->d0_GT != 0.0)
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps, &kappa, dkappa_dg, dkappa_dH);
    PetscReal rhoI_vs_eff        = rho_vs   * (1.0 + user->d0_GT * PetscRealPart(kappa));
    PetscReal d_rhoI_vs_eff_dtem = d_rho_vs * (1.0 + user->d0_GT * PetscRealPart(kappa));
    */

    PetscReal df1;
    DoubleWellDeriv(PetscRealPart(ice_c), NULL, &df1);

    /* air_eff = max(air, air_lim) */
    PetscReal air_r = PetscRealPart(air_c);
    PetscBool air_above_lim = (air_r > air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal air_eff = air_above_lim ? air_r : air_lim;

    /* loc and its ice-derivative (sublimation off; commented out).
    PetscReal loc       = PetscRealPart(ice_c) * PetscRealPart(ice_c) * air_r * air_r;
    PetscReal dloc_dice = 2.0 * PetscRealPart(ice_c) * air_r * (air_r - PetscRealPart(ice_c));
    PetscReal rho_v_minus_rho_vs = PetscRealPart(rhov) - rhoI_vs_eff;
    */

    PetscReal S_T = 1.0 / (rho_ice * lat_sub);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[3][nen][3] = (PetscScalar (*)[3][nen][3])Je;

    for (a = 0; a < nen; a++) {
        PetscReal grad_Na_dot_grad_tem  = 0.0;
        PetscReal grad_Na_dot_grad_rhov = 0.0;
        for (l = 0; l < dim; l++) {
            grad_Na_dot_grad_tem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            grad_Na_dot_grad_rhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }

        for (b = 0; b < nen; b++) {
            PetscReal Na_Nb   = N0[a] * N0[b];
            PetscReal N1a_N1b = 0.0;
            for (l = 0; l < dim; l++) N1a_N1b += N1[a][l] * N1[b][l];

            /* ====================== [ ice , * ] ============================= */
            J[a][0][b][0] += shift * Na_Nb
                           + 3.0 * mob_sub * eps * N1a_N1b
                           + (3.0 * mob_sub / eps) * df1 * Na_Nb;
            /* Sublimation source Jacobian terms (commented out):
            J[a][0][b][0] -= alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb;
            J[a][0][b][1] += alph_sub * loc * d_rhoI_vs_eff_dtem / rho_ice * Na_Nb;
            J[a][0][b][2] -= alph_sub * loc / rho_ice * Na_Nb;
            if (user->d0_GT != 0.0) { ... GT curvature chain-rule block ... }
            */

            /* ====================== [ tem , * ] (row-scaled by S_T) ========= */
            J[a][1][b][0] += S_T * dthcond_dice * grad_Na_dot_grad_tem * N0[b];
            /* Latent-heat coupling (sublimation off; commented out):
            J[a][1][b][0] -= S_T * couple * rho_ice * lat_sub * shift * Na_Nb;
            */
            J[a][1][b][1] += S_T * (
                  shift * rho * cp * Na_Nb
                + thcond * N1a_N1b
            );

            /* ====================== [ vap , * ] ============================= */
            /* The air_above_lim block is d(air_eff)/d(ice) acting on the
             * storage and diffusion terms — geometric coupling, always active. */
            /* Mass-exchange coupling (sublimation off; commented out):
            J[a][2][b][0] += couple * (rho_ice - PetscRealPart(rhov)) * shift * Na_Nb;
            */
            if (air_above_lim) {
                J[a][2][b][0] += -Na_Nb * PetscRealPart(rhov_t)
                               - dif_vap * N0[b] * grad_Na_dot_grad_rhov;
            }
            J[a][2][b][1] += d_dif_vap * air_eff * grad_Na_dot_grad_rhov * N0[b];
            J[a][2][b][2] += air_eff * shift * Na_Nb
                           + dif_vap * air_eff * N1a_N1b;
            /* Mass-exchange Jacobian term (sublimation off; commented out):
            J[a][2][b][2] -= couple * Na_Nb * PetscRealPart(ice_t);
            */
        }
    }
    return 0;
}


PetscErrorCode Jacobian(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Je, void *ctx)
{
    return Jacobian_A1(pnt, shift, V, t, U, Je, ctx);
}


/* =========================================================================
 * Per-element scalar integrals for the monitor table.
 *   S[0] = ice
 *   S[1] = ice^2 air^2   (ice-air interface measure)
 *   S[2] = air
 *   S[3] = tem
 *   S[4] = rhov * air
 * ========================================================================= */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                           PetscScalar *S, void *ctx)
{
    PetscFunctionBegin;
    (void)n; (void)ctx;
    PetscScalar sol[3];
    IGAPointFormValue(pnt, U, &sol[0]);

    PetscReal ice  = PetscRealPart(sol[0]);
    PetscReal tem  = PetscRealPart(sol[1]);
    PetscReal rhov = PetscRealPart(sol[2]);
    PetscReal air  = 1.0 - ice;

    S[0] = ice;
    S[1] = ice*ice * air*air;
    S[2] = air;
    S[3] = tem;
    S[4] = rhov * air;

    PetscFunctionReturn(0);
}
