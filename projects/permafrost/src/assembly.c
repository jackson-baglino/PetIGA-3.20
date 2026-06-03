#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * Clean 2-phase formulation.
 *
 * Strong form (sediment is identically frozen — dphi_s/dt = 0 — so dphi_a = -dphi_i):
 *
 *   dphi_i/dt = -3 M_0 / (eta_i + eta_a) [
 *                 (dF_tri/dphi_i - dF_tri/dphi_a) / eps
 *                 - (eta_i + eta_a) eps grad^2 phi_i
 *                 - eta_a            grad^2 phi_s ]
 *             + S_sub
 *
 *   dphi_s/dt = 0
 *
 *   rho cp dT/dt = div(K grad T) + rho L_sub S_sub
 *
 *   drho_v/dt = div(D_eff grad rho_v)
 *             - k_pen g(phi_i + phi_s) (rho_v - rho_v_eq)
 *             - rho_ice S_sub
 *
 *   D_eff(phi) = D_v(T) g(phi_a) + D_pen (1 - g(phi_a))
 *   g(x)       = x^3 (3 - 2x)         (plain SmoothHeavisidePoly)
 *   rho_v_eq   = (phi_i + phi_s) rho_vs(T) + phi_a rho_v
 *   S_sub      = alpha_sub * (phi_i*phi_a)^2 * (rho_v - rho_vs(T)) / rho_ice
 *
 *   phi_a = 1 - phi_i - phi_s
 *
 * Sign convention:
 *   S_sub > 0 = deposition (vapor -> ice). Ice gains (+S_sub),
 *   vapor loses (-rho_ice S_sub), heat is released (+rho L_sub S_sub).
 *
 * The temperature residual is row-scaled by S_T = 1 / (rho_ice L_sub) — purely
 * numerical preconditioning. The physics is unchanged (dividing both sides of
 * the T equation by a constant); the residual the solver sees is O(sub_src)
 * instead of O(rho_ice L_sub sub_src ~ 1e9 sub_src).
 * ========================================================================= */


/* -------------------------------------------------------------------------
 * Sed-direction derivatives of the chemical-potential factors fi, fa, fs.
 *   fi = eta_i phi_i(1-phi_i)(1-2 phi_i) + 2 Lambda phi_i phi_s^2 phi_a^2
 *   fa = eta_a phi_a(1-phi_a)(1-2 phi_a) + 2 Lambda phi_i^2 phi_s^2 phi_a
 *   fs = eta_s phi_s(1-phi_s)(1-2 phi_s) + 2 Lambda phi_i^2 phi_s phi_a^2
 *
 * Fice / Fair / Fsed in material_properties.c return the d/dphi_i derivatives;
 * the d/dphi_s derivatives are computed here (treating phi_a = 1-phi_i-phi_s,
 * so dphi_a/dphi_s = -1).
 * ------------------------------------------------------------------------- */
static void ChemPot_dsed(AppCtx *user,
                         PetscReal ice, PetscReal sed,
                         PetscReal *dfi_dsed_out,
                         PetscReal *dfs_dsed_out,
                         PetscReal *dfa_dsed_out)
{
    PetscReal Etaa   = user->Etaa;
    PetscReal Etased = user->Etam;
    PetscReal Lambd  = user->Lambd;
    PetscReal air    = 1.0 - ice - sed;

    /* dfi/ds = 4 Lambda phi_i phi_s phi_a (phi_a - phi_s) */
    if (dfi_dsed_out)
        *dfi_dsed_out = 4.0 * Lambd * ice * sed * air * (air - sed);

    /* dfs/ds = eta_s (1 - 6 phi_s + 6 phi_s^2) + 2 Lambda phi_i^2 phi_a (phi_a - 2 phi_s) */
    if (dfs_dsed_out)
        *dfs_dsed_out = Etased * (1.0 - 6.0*sed + 6.0*sed*sed)
                      + 2.0 * Lambd * ice*ice * air * (air - 2.0*sed);

    /* dfa/ds = -eta_a (1 - 6 phi_a + 6 phi_a^2) + 2 Lambda phi_i^2 phi_s (2 phi_a - phi_s) */
    if (dfa_dsed_out)
        *dfa_dsed_out = -Etaa * (1.0 - 6.0*air + 6.0*air*air)
                      + 2.0 * Lambd * ice*ice * sed * (2.0*air - sed);
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
    PetscReal eps     = user->eps;
    PetscReal Etai    = user->Etai;
    PetscReal Etaa    = user->Etaa;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal alph_sub= user->alph_sub;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice   = sol[0],  ice_t  = sol_t[0];
    PetscScalar tem   = sol[1],  tem_t  = sol_t[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    PetscScalar sed   = sol[3],  sed_t  = sol_t[3];
    PetscScalar grad_ice [dim], grad_sed[dim];
    PetscScalar grad_tem [dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_ice [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
        grad_sed [l] = grad_sol[3][l];
    }
    PetscScalar air = 1.0 - ice - sed;

    /* SNES domain-error catch — if a trial Newton iterate has phi out of the
     * configured bounds, signal an invalid state so the line search backs off
     * the step. Cheaper than committing a bad state and rolling back later. */
    {
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (PetscRealPart(ice) < lo || PetscRealPart(ice) > hi ||
            PetscRealPart(sed) < lo || PetscRealPart(sed) > hi ||
            PetscRealPart(air) < lo || PetscRealPart(air) > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    /* Clamp clipping copies of phi for material-property evaluation only —
     * keeps Density/HeatCap/etc. from getting negative inputs near the bounds. */
    PetscScalar ice_c = ice, sed_c = sed, air_c = air;
    if (PetscRealPart(ice_c) < 0.0) ice_c = 0.0;
    if (PetscRealPart(ice_c) > 1.0) ice_c = 1.0;
    if (PetscRealPart(sed_c) < 0.0) sed_c = 0.0;
    if (PetscRealPart(sed_c) > 1.0) sed_c = 1.0;
    if (PetscRealPart(air_c) < 0.0) air_c = 0.0;
    if (PetscRealPart(air_c) > 1.0) air_c = 1.0;

    /* Material properties */
    PetscReal thcond, cp, rho, dif_vap, rho_vs, fi, fa, mob_sub;
    ThermalCond(user, ice_c, sed_c, &thcond,  NULL);
    HeatCap    (user, ice_c, sed_c, &cp,      NULL);
    Density    (user, ice_c, sed_c, &rho,     NULL);
    VaporDiffus(user, tem,           &dif_vap, NULL);
    RhoVS_I    (user, tem,           &rho_vs,  NULL);
    Fice       (user, ice_c, sed_c, &fi,      NULL);
    Fair       (user, ice_c, sed_c, &fa,      NULL);
    Mobility   (user, ice_c, sed_c, &mob_sub);

    /* Effective vapor diffusivity: D_eff = D_v g(phi_a) + D_pen (1 - g(phi_a)) */
    PetscScalar g_a;
    SmoothHeavisidePoly(air_c, &g_a, NULL);
    PetscReal D_pen = user->difvap_pen * dif_vap;
    PetscReal D_eff = (PetscReal)g_a * dif_vap + (1.0 - (PetscReal)g_a) * D_pen;

    /* Penalty: -k_pen g(phi_i + phi_s) (rho_v - rho_v_eq)
     *          rho_v_eq = (phi_i + phi_s) rho_vs + phi_a rho_v   (self-referential in air) */
    PetscScalar g_solid;
    SmoothHeavisidePoly(ice_c + sed_c, &g_solid, NULL);
    PetscScalar rhov_eq = (ice_c + sed_c) * rho_vs + air_c * rhov;
    PetscScalar vap_pen = user->k_pen * g_solid * (rhov - rhov_eq);

    /* Phase-change source. S_sub > 0 = deposition (vapor -> ice). */
    PetscReal loc   = PetscRealPart(ice_c) * PetscRealPart(ice_c)
                    * PetscRealPart(air_c) * PetscRealPart(air_c);
    PetscReal S_sub = alph_sub * loc * (PetscRealPart(rhov) - rho_vs) / rho_ice;

    /* Ice driving-force prefactor: K2P = 3 M_0 / (eta_i + eta_a). */
    PetscReal K2P = 3.0 * mob_sub / (Etai + Etaa);

    /* T-equation row-scale (numerical preconditioning; no physics change). */
    PetscReal S_T = 1.0 / (rho_ice * lat_sub);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal grad_N_dot_grad_ice  = 0.0;
        PetscReal grad_N_dot_grad_sed  = 0.0;
        PetscReal grad_N_dot_grad_tem  = 0.0;
        PetscReal grad_N_dot_grad_rhov = 0.0;
        for (l = 0; l < dim; l++) {
            grad_N_dot_grad_ice  += N1[a][l] * PetscRealPart(grad_ice [l]);
            grad_N_dot_grad_sed  += N1[a][l] * PetscRealPart(grad_sed [l]);
            grad_N_dot_grad_tem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            grad_N_dot_grad_rhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }

        /* --- Ice ----------------------------------------------------- */
        PetscScalar R_ice = N0[a] * ice_t
                          + 3.0 * mob_sub * eps * grad_N_dot_grad_ice
                          + K2P * Etaa      * grad_N_dot_grad_sed
                          + (K2P / eps) * (fi - fa) * N0[a]
                          - S_sub * N0[a];

        /* --- Sediment: identically frozen (R_sed = N0 * sed_t forces sed_t = 0) */
        PetscScalar R_sed = N0[a] * sed_t;

        /* --- Temperature (row-scaled by S_T) --------------------------
         * Strong form: rho cp dT/dt = div(K grad T) + rho L_sub S_sub.
         * Residual = LHS - RHS, so the latent-heat term enters with a MINUS sign.
         * After IBP on the divergence: R_tem = N rho cp tem_t + K (grad_N . grad_T)
         *                                   - rho L_sub S_sub N. */
        PetscScalar R_tem = rho * cp * N0[a] * tem_t
                          + thcond * grad_N_dot_grad_tem
                          - rho * lat_sub * S_sub * N0[a];

        /* --- Vapor ---------------------------------------------------- */
        PetscScalar R_vap = N0[a] * rhov_t
                          + D_eff * grad_N_dot_grad_rhov
                          + vap_pen * N0[a]
                          + rho_ice * S_sub * N0[a];

        R[a][0] = R_ice;
        R[a][1] = S_T * R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;
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
 * DOF layout: 0=ice, 1=tem, 2=rhov, 3=sed.
 *
 * Sub_src localiser:  loc(phi_i, phi_a) = phi_i^2 phi_a^2
 *   dloc/dphi_i = 2 phi_i phi_a (phi_a - phi_i)         [chain through phi_a = 1-phi_i-phi_s]
 *   dloc/dphi_s = -2 phi_i^2 phi_a
 *
 * rho_v_eq derivatives (self-referential form: rho_v_eq = (phi_i+phi_s) rho_vs + phi_a rho_v):
 *   d(rho_v_eq)/dphi_i = rho_vs - rho_v    (because dphi_a/dphi_i = -1)
 *   d(rho_v_eq)/dphi_s = rho_vs - rho_v
 *   d(rho_v_eq)/drho_v = phi_a
 *   d(rho_v_eq)/dT     = (phi_i + phi_s) d_rho_vs/dT
 * ========================================================================= */
static PetscErrorCode Jacobian_A1(IGAPoint pnt,
                                  PetscReal shift, const PetscScalar *V,
                                  PetscReal t, const PetscScalar *U,
                                  PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    PetscReal eps     = user->eps;
    PetscReal Etai    = user->Etai;
    PetscReal Etaa    = user->Etaa;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal alph_sub= user->alph_sub;

    (void)V;  /* time derivatives enter the Jacobian only through the `shift` factor */

    if (pnt->atboundary) return 0;

    PetscScalar sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice  = sol[0];
    PetscScalar tem  = sol[1];
    PetscScalar rhov = sol[2];
    PetscScalar sed  = sol[3];
    PetscScalar grad_tem[dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar air = 1.0 - ice - sed;

    /* Clamp copies for material-property evaluation */
    PetscScalar ice_c = ice, sed_c = sed, air_c = air;
    if (PetscRealPart(ice_c) < 0.0) ice_c = 0.0;
    if (PetscRealPart(ice_c) > 1.0) ice_c = 1.0;
    if (PetscRealPart(sed_c) < 0.0) sed_c = 0.0;
    if (PetscRealPart(sed_c) > 1.0) sed_c = 1.0;
    if (PetscRealPart(air_c) < 0.0) air_c = 0.0;
    if (PetscRealPart(air_c) > 1.0) air_c = 1.0;

    /* Material properties + derivatives */
    PetscReal thcond, cp, rho, dif_vap, rho_vs, fi, fa, mob_sub;
    ThermalCond(user, ice_c, sed_c, &thcond,  NULL);
    HeatCap    (user, ice_c, sed_c, &cp,      NULL);
    Density    (user, ice_c, sed_c, &rho,     NULL);
    VaporDiffus(user, tem,           &dif_vap, NULL);
    RhoVS_I    (user, tem,           &rho_vs,  NULL);
    Fice       (user, ice_c, sed_c, &fi,      NULL);
    Fair       (user, ice_c, sed_c, &fa,      NULL);
    Mobility   (user, ice_c, sed_c, &mob_sub);

    PetscReal dthcond_dice, dthcond_dsed;
    ThermalCond(user, ice_c, sed_c, NULL, &dthcond_dice);
    /* HeatCap, Density currently report d/d(ice) only via second arg; that
     * matches the prior treatment of rho, cp as "frozen" for the latent-heat
     * derivative. Conduction's d/d(sed) is implicit in dthcond_dice via the
     * weighted-average form: dthcond/dsed has a separate algebraic form, but
     * the existing API only exposes the ice-direction derivative, so we follow
     * the existing modified-Newton convention there. */
    dthcond_dsed = 0.0;

    PetscReal d_dif_vap, d_rho_vs;
    VaporDiffus(user, tem, NULL, &d_dif_vap);
    RhoVS_I    (user, tem, NULL, &d_rho_vs);

    PetscReal dfi_dice, dfa_dice;
    Fice(user, ice_c, sed_c, NULL, &dfi_dice);
    Fair(user, ice_c, sed_c, NULL, &dfa_dice);

    PetscReal dfi_dsed, dfs_dsed, dfa_dsed;
    ChemPot_dsed(user, PetscRealPart(ice_c), PetscRealPart(sed_c),
                 &dfi_dsed, &dfs_dsed, &dfa_dsed);

    /* Smooth Heaviside on air and on (ice + sed) */
    PetscScalar g_a, dg_a;
    SmoothHeavisidePoly(air_c, &g_a, &dg_a);
    PetscScalar g_solid, dg_solid;
    SmoothHeavisidePoly(ice_c + sed_c, &g_solid, &dg_solid);

    PetscReal D_pen   = user->difvap_pen * dif_vap;
    PetscReal D_eff   = (PetscReal)g_a * dif_vap + (1.0 - (PetscReal)g_a) * D_pen;
    /* d(D_eff)/dphi_a = (D_v - D_pen) * g_a';   dphi_a/dphi_i = -1, dphi_a/dphi_s = -1
     *   => d(D_eff)/dphi_i = d(D_eff)/dphi_s = -(D_v - D_pen) g_a'              */
    PetscReal dDeff_dair  = (dif_vap - D_pen) * (PetscReal)dg_a;
    PetscReal dDeff_dice  = -dDeff_dair;
    PetscReal dDeff_dsed  = -dDeff_dair;
    /* d(D_eff)/dT: D_v depends on T; the linear combination scales: */
    PetscReal dDeff_dtem  = d_dif_vap * ((PetscReal)g_a + (1.0 - (PetscReal)g_a) * user->difvap_pen);

    /* rho_v_eq = (phi_i + phi_s) rho_vs + phi_a rho_v */
    PetscReal rhov_eq = (PetscRealPart(ice_c) + PetscRealPart(sed_c)) * rho_vs
                      + PetscRealPart(air_c) * PetscRealPart(rhov);
    PetscReal d_rhoveq_dice = rho_vs - PetscRealPart(rhov);
    PetscReal d_rhoveq_dsed = rho_vs - PetscRealPart(rhov);
    PetscReal d_rhoveq_drv  = PetscRealPart(air_c);
    PetscReal d_rhoveq_dtem = (PetscRealPart(ice_c) + PetscRealPart(sed_c)) * d_rho_vs;

    /* sub_src and its derivatives */
    PetscReal loc       = PetscRealPart(ice_c) * PetscRealPart(ice_c)
                        * PetscRealPart(air_c) * PetscRealPart(air_c);
    PetscReal dloc_dice = 2.0 * PetscRealPart(ice_c) * PetscRealPart(air_c)
                        * (PetscRealPart(air_c) - PetscRealPart(ice_c));
    PetscReal dloc_dsed = -2.0 * PetscRealPart(ice_c) * PetscRealPart(ice_c)
                        * PetscRealPart(air_c);
    PetscReal rho_v_minus_rho_vs = PetscRealPart(rhov) - rho_vs;

    PetscReal K2P = 3.0 * mob_sub / (Etai + Etaa);
    PetscReal S_T = 1.0 / (rho_ice * lat_sub);
    PetscReal k_pen = user->k_pen;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;

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

            /* ====================== [ ice , * ] ====================== */
            /* Time deriv + AC gradient + bulk drive + sub_src.
             * dR_ice/dphi_i^b: shift mass + 3 M_0 eps grad_grad + (K2P/eps) (dfi-dfa) Na_Nb
             *                  - alpha_sub dloc_dice (rho_v - rho_vs)/rho_ice Na_Nb */
            J[a][0][b][0] += shift * Na_Nb
                           + 3.0 * mob_sub * eps * N1a_N1b
                           + (K2P / eps) * (dfi_dice - dfa_dice) * Na_Nb
                           - alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb;

            /* dR_ice/dT^b: through rho_vs(T) in sub_src (note sign — sub_src appears
             *               in R_ice as -S_sub so d(-S_sub)/dT = +alpha_sub loc d_rho_vs/rho_ice) */
            J[a][0][b][1] += alph_sub * loc * d_rho_vs / rho_ice * Na_Nb;

            /* dR_ice/drho_v^b */
            J[a][0][b][2] += -alph_sub * loc / rho_ice * Na_Nb;

            /* dR_ice/dphi_s^b: cross-coupling through (1) the grad^2(phi_s) term in
             *                  the ice equation, (2) the (fi-fa) bulk via ChemPot_dsed,
             *                  (3) sub_src through dloc_dsed. */
            J[a][0][b][3] += K2P * Etaa * N1a_N1b
                           + (K2P / eps) * (dfi_dsed - dfa_dsed) * Na_Nb
                           - alph_sub * dloc_dsed * rho_v_minus_rho_vs / rho_ice * Na_Nb;

            /* ====================== [ sed , * ] ====================== */
            /* R_sed = N0 sed_t   =>   only [sed,sed] gets a mass-matrix contribution */
            J[a][3][b][3] += shift * Na_Nb;

            /* ====================== [ tem , * ] (row-scaled by S_T) =====
             * R_tem has -rho*L_sub*S_sub * N, so d(R_tem)/d(state) for the latent
             * heat enters with the same minus sign relative to the corresponding
             * sub_src derivatives. */
            /* dR_tem/dphi_i^b: conduction through dK/dphi_i + latent through dloc_dice */
            J[a][1][b][0] += S_T * (
                  dthcond_dice * grad_Na_dot_grad_tem * N0[b]
                - rho * lat_sub * alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb
            );

            /* dR_tem/dT^b: shift mass + conduction stiffness + latent through d_rho_vs */
            J[a][1][b][1] += S_T * (
                  shift * rho * cp * Na_Nb
                + thcond * N1a_N1b
                + rho * lat_sub * alph_sub * loc * d_rho_vs / rho_ice * Na_Nb
            );

            /* dR_tem/drho_v^b: latent through d(sub_src)/d(rhov) */
            J[a][1][b][2] += S_T * (
                - rho * lat_sub * alph_sub * loc / rho_ice * Na_Nb
            );

            /* dR_tem/dphi_s^b: conduction through dK/dphi_s + latent through dloc_dsed */
            J[a][1][b][3] += S_T * (
                  dthcond_dsed * grad_Na_dot_grad_tem * N0[b]
                - rho * lat_sub * alph_sub * dloc_dsed * rho_v_minus_rho_vs / rho_ice * Na_Nb
            );

            /* ====================== [ vap , * ] ====================== */
            /* Diffusion through D_eff(phi_a, T), penalty through g_solid and rho_v_eq,
             * sub_src source. */
            /* dR_vap/dphi_i^b:
             *   - dDeff/dphi_i changes the flux stiffness
             *   - dg_solid/dphi_i affects penalty weight
             *   - d(rho_v_eq)/dphi_i affects penalty value (= rho_vs - rho_v)
             *   - sub_src term has +rho_ice S_sub, derivative is +alpha_sub*dloc_dice*(rho_v-rho_vs) Na_Nb
             *     (the rho_ice cancels with the /rho_ice in S_sub).                       */
            J[a][2][b][0] += dDeff_dice * grad_Na_dot_grad_rhov * N0[b]
                           + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                           - k_pen * (PetscReal)g_solid  * d_rhoveq_dice * Na_Nb
                           + alph_sub * dloc_dice * rho_v_minus_rho_vs * Na_Nb;

            /* dR_vap/dT^b: D_v(T) in diffusion; rho_vs(T) in penalty (via rho_v_eq);
             *               rho_vs(T) in sub_src. */
            J[a][2][b][1] += dDeff_dtem * grad_Na_dot_grad_rhov * N0[b]
                           - k_pen * (PetscReal)g_solid * d_rhoveq_dtem * Na_Nb
                           - alph_sub * loc * d_rho_vs * Na_Nb;

            /* dR_vap/drho_v^b: shift mass + diffusion stiffness + penalty (1 - phi_a) + sub_src linearity */
            J[a][2][b][2] += shift * Na_Nb
                           + D_eff * N1a_N1b
                           + k_pen * (PetscReal)g_solid * (1.0 - d_rhoveq_drv) * Na_Nb
                           + alph_sub * loc * Na_Nb;

            /* dR_vap/dphi_s^b: same structural pattern as [vap, ice] (g_solid depends on
             *                   phi_i+phi_s, D_eff on phi_a = 1-phi_i-phi_s). */
            J[a][2][b][3] += dDeff_dsed * grad_Na_dot_grad_rhov * N0[b]
                           + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                           - k_pen * (PetscReal)g_solid  * d_rhoveq_dsed * Na_Nb
                           + alph_sub * dloc_dsed * rho_v_minus_rho_vs * Na_Nb;
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
 * (Unchanged from prior; columns: ice, ice^2 sed^2 air^2, air, T, rho_v*air,
 *  ice^2 air^2, sed, sed^2 air^2, sed^2 ice^2.)
 * ========================================================================= */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                           PetscScalar *S, void *ctx)
{
    PetscFunctionBegin;
    PetscScalar sol[4];
    IGAPointFormValue(pnt, U, &sol[0]);

    PetscReal ice  = PetscRealPart(sol[0]);
    PetscReal tem  = PetscRealPart(sol[1]);
    PetscReal rhov = PetscRealPart(sol[2]);
    PetscReal sed  = PetscRealPart(sol[3]);
    PetscReal air  = 1.0 - sed - ice;

    S[0] = ice;
    S[1] = SQ(air) * SQ(sed) * SQ(ice);
    S[2] = air;
    S[3] = tem;
    S[4] = rhov * air;
    S[5] = air*air * ice*ice;
    S[6] = sed;
    S[7] = sed*sed * air*air;
    S[8] = sed*sed * ice*ice;

    PetscFunctionReturn(0);
}
