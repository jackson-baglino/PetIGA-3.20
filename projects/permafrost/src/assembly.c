#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * 2-phase (ice / temperature / vapor) system. phi_a = 1 - phi_i is
 * computed algebraically; there is no sediment DOF.
 *
 * Free energy: single double-well in phi_i,
 *   f1(phi_i) = phi_i*(1-phi_i)*(1-2*phi_i)
 *
 * Ice (Allen-Cahn):
 *   R_ice = N0[a]*ice_t + 3*mob_sub*eps*grad_N.grad_ice
 *         + (3*mob_sub/eps)*f1(phi_i)*N0[a] - S_sub*N0[a]
 *
 * This is the met=0 reduction of dry_snow_metamorphism's 3-phase ice
 * equation: R_ice += mob*3/eps/ETA * ((Etam+Etaa)*fice - Etaa*fmet - Etam*fair).
 * With met=0, ETA=Etaa*Etai, fmet=0, fice=Etai*f1(ice), so the expression
 * collapses to mob*3/eps*f1(ice) -- Sigma_i/Sigma_a/Lambda all cancel and
 * play no role in a true 2-phase (no third phase) system.
 *
 * Temperature and vapor are sourced by the FULL d(phi_i)/dt = ice_t (not
 * just the kinetic S_sub term) per Moure & Fu (2024), SI Eqs. (6)-(7) --
 * the sublimation-equivalence reduction of the wet-snow model. ice_t
 * already includes the AC curvature/double-well relaxation, so using it
 * (rather than S_sub alone) makes the ice<->vapor mass exchange and the
 * latent-heat release exact, regardless of mob_sub. This assumes ice_t is
 * always a real ice<->vapor exchange; -decouple_phase_change 1 zeroes this
 * coupling (multiplies it by `couple`, 0 or 1) for diagnostics that ablate
 * -alph_sub but leave mob_sub active, where ice_t is pure curvature
 * relaxation, not phase change -- without decoupling, that motion still
 * gets dumped into vapor at full rho_ice, blowing up the vapor field:
 *
 * Temperature (row-scaled by S_T = 1/(rho_ice*lat_sub), numerical
 * preconditioning only). The latent-heat source uses the constant rho_ice
 * (mass of the converting ice phase, matching R_vap's rho_SE=rho_ice), NOT
 * the local mixture density rho(phi_i) -- otherwise the energy exchanged
 * is inconsistent with the mass exchanged in R_vap at the diffuse interface:
 *   R_tem = rho*cp*N0[a]*tem_t + thcond*grad_N.grad_tem - rho_ice*lat_sub*ice_t*N0[a]
 *
 * Vapor (air_lim floor as in dry_snow_metamorphism):
 *   R_vap = N0[a]*air_eff*rhov_t + difvap*air_eff*grad_N.grad_rhov
 *         + N0[a]*(rho_ice - rhov)*ice_t
 *   air_eff = (air > air_lim) ? air : air_lim,   air = 1 - ice
 *
 * Phase-change source (S_sub > 0 = deposition, vapor -> ice; appears only
 * in R_ice, which defines ice_t):
 *   loc = phi_i^2 * phi_a^2
 *   S_sub = alph_sub * loc * (rho_v - rho_vs(T)) / rho_ice
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
    PetscReal alph_sub= user->alph_sub;
    PetscReal air_lim = user->air_lim;
    /* -decouple_phase_change 1: zero the ice_t-driven source terms in R_tem
     * and R_vap (see file header). Lets -alph_sub 0 truly isolate pure AC
     * curvature relaxation -- without this, ice motion from curvature
     * relaxation alone (not real sublimation) still gets dumped into vapor
     * at full rho_ice, producing an unphysical vapor blowup amplified by
     * rho_ice/rho_v (~1e5-1e6x). */
    PetscReal couple = user->decouple_phase_change ? 0.0 : 1.0;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice   = sol[0],  ice_t  = sol_t[0];
    PetscScalar tem   = sol[1],  tem_t  = sol_t[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    PetscScalar grad_ice [dim];
    PetscScalar grad_tem [dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_ice [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar air = 1.0 - ice;

    /* Hessian of ice, needed only for the Gibbs-Thomson curvature term
     * below (kappa = -div(grad_ice/|grad_ice|)). Read unconditionally; the
     * GT correction itself is skipped when d0_GT == 0, but the read is cheap. */
    PetscScalar hess_sol[3][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);
    PetscScalar hess_ice[dim * dim];
    for (l = 0; l < dim * dim; l++) hess_ice[l] = hess_sol[0][l / dim][l % dim];

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
    PetscReal thcond, cp, rho, dif_vap, rho_vs, mob_sub;
    ThermalCond(user, ice_c, &thcond,  NULL);
    HeatCap    (user, ice_c, &cp,      NULL);
    Density    (user, ice_c, &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    RhoVS_I    (user, tem,    &rho_vs,  NULL);
    Mobility   (user, ice_c, &mob_sub);

    /* Gibbs-Thomson curvature correction to the saturation vapor density:
     * rhoI_vs_eff = rho_vs * (1 + d0_GT * kappa), kappa = -div(grad_ice/|grad_ice|).
     * Without this, every interface shares one flat rho_vs(T) and there is no
     * driving force for Ostwald ripening: larger grains/flatter necks (smaller
     * kappa) get a lower rhoI_vs_eff than smaller/sharper ones, so vapor flows
     * small -> large. d0_GT == 0 disables it (rhoI_vs_eff = rho_vs). */
    PetscScalar kappa = 0.0;
    if (user->d0_GT != 0.0)
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps, &kappa, NULL, NULL);
    PetscReal rhoI_vs_eff = rho_vs * (1.0 + user->d0_GT * PetscRealPart(kappa));

    PetscReal f1;
    DoubleWellDeriv(PetscRealPart(ice_c), &f1, NULL);

    /* air_eff = max(air, air_lim) */
    PetscReal air_r = PetscRealPart(air_c);
    PetscReal air_eff = (air_r > air_lim) ? air_r : air_lim;

    /* Phase-change source. S_sub > 0 = deposition (vapor -> ice). */
    PetscReal loc   = PetscRealPart(ice_c) * PetscRealPart(ice_c)
                    * air_r * air_r;
    PetscReal S_sub = alph_sub * loc * (PetscRealPart(rhov) - rhoI_vs_eff) / rho_ice;

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
                          + (3.0 * mob_sub / eps) * f1 * N0[a]
                          - S_sub * N0[a];

        PetscScalar R_tem = rho * cp * N0[a] * tem_t
                          + thcond * grad_N_dot_grad_tem
                          - couple * rho_ice * lat_sub * ice_t * N0[a];

        PetscScalar R_vap = N0[a] * air_eff * rhov_t
                          + dif_vap * air_eff * grad_N_dot_grad_rhov
                          + couple * N0[a] * (rho_ice - rhov) * ice_t;

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
 *
 * loc(phi_i, phi_a) = phi_i^2 phi_a^2,  phi_a = 1-phi_i
 *   dloc/dphi_i = 2 phi_i phi_a (phi_a - phi_i)
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
    PetscReal alph_sub= user->alph_sub;
    PetscReal air_lim = user->air_lim;
    PetscReal couple  = user->decouple_phase_change ? 0.0 : 1.0;  /* see Residual_A1 */

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice    = sol[0],  ice_t  = sol_t[0];
    PetscScalar tem    = sol[1];
    PetscScalar rhov   = sol[2],  rhov_t = sol_t[2];
    PetscScalar grad_ice[dim];
    PetscScalar grad_tem[dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_ice [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar air = 1.0 - ice;

    /* Hessian of ice for the GT-curvature derivatives (dkappa/dgrad, dkappa/dHess). */
    PetscScalar hess_sol[3][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);
    PetscScalar hess_ice[dim * dim];
    for (l = 0; l < dim * dim; l++) hess_ice[l] = hess_sol[0][l / dim][l % dim];

    /* Clamp copies for material-property evaluation */
    PetscScalar ice_c = ice, air_c = air;
    if (PetscRealPart(ice_c) < 0.0) ice_c = 0.0;
    if (PetscRealPart(ice_c) > 1.0) ice_c = 1.0;
    if (PetscRealPart(air_c) < 0.0) air_c = 0.0;
    if (PetscRealPart(air_c) > 1.0) air_c = 1.0;

    /* Material properties + derivatives */
    PetscReal thcond, cp, rho, dif_vap, rho_vs, mob_sub;
    ThermalCond(user, ice_c, &thcond,  NULL);
    HeatCap    (user, ice_c, &cp,      NULL);
    Density    (user, ice_c, &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    RhoVS_I    (user, tem,    &rho_vs,  NULL);
    Mobility   (user, ice_c, &mob_sub);

    PetscReal dthcond_dice;
    ThermalCond(user, ice_c, NULL, &dthcond_dice);

    PetscReal d_dif_vap, d_rho_vs;
    VaporDiffus(user, tem, NULL, &d_dif_vap);
    RhoVS_I    (user, tem, NULL, &d_rho_vs);

    /* Gibbs-Thomson curvature correction (see Residual_A1 for rationale).
     * kappa doesn't depend on tem, so d(rhoI_vs_eff)/dtem = d_rho_vs*(1+d0_GT*kappa). */
    PetscScalar kappa = 0.0;
    PetscScalar dkappa_dg[3] = {0.0, 0.0, 0.0};
    PetscScalar dkappa_dH[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (user->d0_GT != 0.0)
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps, &kappa, dkappa_dg, dkappa_dH);
    PetscReal rhoI_vs_eff      = rho_vs   * (1.0 + user->d0_GT * PetscRealPart(kappa));
    PetscReal d_rhoI_vs_eff_dtem = d_rho_vs * (1.0 + user->d0_GT * PetscRealPart(kappa));

    PetscReal df1;
    DoubleWellDeriv(PetscRealPart(ice_c), NULL, &df1);

    /* air_eff = max(air, air_lim) */
    PetscReal air_r = PetscRealPart(air_c);
    PetscBool air_above_lim = (air_r > air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal air_eff = air_above_lim ? air_r : air_lim;

    /* loc and its ice-derivative */
    PetscReal loc       = PetscRealPart(ice_c) * PetscRealPart(ice_c) * air_r * air_r;
    PetscReal dloc_dice = 2.0 * PetscRealPart(ice_c) * air_r * (air_r - PetscRealPart(ice_c));
    PetscReal rho_v_minus_rho_vs = PetscRealPart(rhov) - rhoI_vs_eff;

    PetscReal S_T = 1.0 / (rho_ice * lat_sub);

    const PetscReal *N0, (*N1)[dim], (*N2)[dim][dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    /* N2 (second-derivative shape funs) -- needed only for the GT curvature
     * chain-rule term below, but cheap to fetch unconditionally. */
    IGAPointGetShapeFuns(pnt, 2, (const PetscReal**)&N2);

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
                           + (3.0 * mob_sub / eps) * df1 * Na_Nb
                           - alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb;
            J[a][0][b][1] += alph_sub * loc * d_rhoI_vs_eff_dtem / rho_ice * Na_Nb;
            J[a][0][b][2] += -alph_sub * loc / rho_ice * Na_Nb;

            /* Gibbs-Thomson curvature contribution to [ice, ice]. S_sub
             * depends on grad(phi_ice) and Hess(phi_ice) through
             * rhoI_vs_eff = rhoI_vs*(1 + d0_GT*kappa). Chain rule:
             *   d S_sub / d phi_ice^b |_GT
             *     = -alph_sub*loc*rhoI_vs*d0_GT/rho_ice * dkappa/dphi_ice^b,
             *   dkappa/dphi_ice^b = sum_l dkappa_dg[l]*N1[b][l]
             *                     + sum_kl dkappa_dH[kl]*N2[b][k][l].
             * R_ice -= S_sub*N0[a], so J[ice,ice] += -d(S_sub)/dphi*N0[a]
             *   = +alph_sub*loc*rhoI_vs*d0_GT/rho_ice * N0[a]*(grad+hess). */
            if (user->d0_GT != 0.0) {
                PetscReal grad_part = 0.0;
                for (l = 0; l < dim; l++)
                    grad_part += (PetscReal)dkappa_dg[l] * N1[b][l];
                PetscReal hess_part = 0.0;
                PetscInt _k, _l;
                for (_k = 0; _k < dim; _k++)
                    for (_l = 0; _l < dim; _l++)
                        hess_part += (PetscReal)dkappa_dH[_k * dim + _l] * N2[b][_k][_l];
                J[a][0][b][0] += alph_sub * loc * rho_vs * user->d0_GT
                              / rho_ice * N0[a] * (grad_part + hess_part);
            }

            /* ====================== [ tem , * ] (row-scaled by S_T) ========= */
            J[a][1][b][0] += S_T * (
                  dthcond_dice * grad_Na_dot_grad_tem * N0[b]
                - couple * rho_ice * lat_sub * shift * Na_Nb
            );
            J[a][1][b][1] += S_T * (
                  shift * rho * cp * Na_Nb
                + thcond * N1a_N1b
            );

            /* ====================== [ vap , * ] ============================= */
            /* Only the d/d(ice) of the couple*(rho_ice-rhov)*ice_t source term
             * is gated here -- the air_above_lim block below is d(air_eff)/d(ice)
             * acting on the (unconditional) mass-matrix/diffusion terms, a
             * separate, legitimate geometric coupling that stays on regardless
             * of -decouple_phase_change. */
            J[a][2][b][0] += couple * (rho_ice - PetscRealPart(rhov)) * shift * Na_Nb;
            if (air_above_lim) {
                J[a][2][b][0] += -Na_Nb * PetscRealPart(rhov_t)
                               - dif_vap * N0[b] * grad_Na_dot_grad_rhov;
            }
            J[a][2][b][1] += d_dif_vap * air_eff * grad_Na_dot_grad_rhov * N0[b];
            J[a][2][b][2] += -couple * Na_Nb * PetscRealPart(ice_t)
                           + air_eff * shift * Na_Nb
                           + dif_vap * air_eff * N1a_N1b;
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
