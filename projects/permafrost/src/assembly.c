#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * Two formulations, switched by user->flag_relax.
 *
 * ---------- Three-phase relaxation window (flag_relax == TRUE) -------------
 * Full beta-eliminated Kim-Steinbach AC: dphi_i + dphi_s + dphi_a = 0.
 *
 *   dphi_i/dt = -3 M_0 / (eps Sigma_T) [
 *                  (Sigma_s + Sigma_a) dF/dphi_i
 *                  - Sigma_s            dF/dphi_a
 *                  - Sigma_a            dF/dphi_s ]
 *             + 3 M_0 eps grad^2 phi_i  + S_sub
 *
 *   dphi_s/dt = -3 M_0 / (eps Sigma_T) [
 *                  -Sigma_a            dF/dphi_i
 *                  -Sigma_i            dF/dphi_a
 *                  -(Sigma_i+Sigma_a)  dF/dphi_s ]
 *             + 3 M_0 eps grad^2 phi_s
 *
 *   rho cp dT/dt = div(K grad T) - rho L_sub dphi_a/dt
 *                = div(K grad T) + rho L_sub (dphi_i/dt + dphi_s/dt)
 *
 *   drho_v/dt = div(D_eff grad rho_v)
 *             - k_pen g(phi_i + phi_s) (rho_v - rho_v_eq)
 *             - rho_i dphi_i/dt
 *
 * Note: the 3-phase equations have NO grad^2(other phase) cross-coupling —
 * that only emerges after beta-elimination under the 2-phase constraint.
 *
 * ---------- Two-phase post-relax window (flag_relax == FALSE) --------------
 * Sediment is identically frozen (dphi_s/dt = 0), so dphi_a = -dphi_i.
 *
 *   dphi_i/dt = -3 M_0 / (eta_i + eta_a) [
 *                 (dF/dphi_i - dF/dphi_a) / eps
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
 * The grad^2(phi_s) cross-coupling in the post-relax ice equation is the
 * algebraic trace of the 3-phase derivation under the dphi_s/dt = 0
 * constraint.
 *
 * ---------- Common pieces ---------------------------------------------------
 *   D_eff(phi) = D_v(T) g(phi_a) + D_pen (1 - g(phi_a))
 *   g(x)       = x^3 (3 - 2x)         (plain SmoothHeavisidePoly)
 *   rho_v_eq   = (phi_i + phi_s) rho_vs(T) + phi_a rho_v
 *   S_sub      = alpha_sub * (phi_i*phi_a)^2 * (rho_v - rho_vs(T)) / rho_ice
 *   phi_a      = 1 - phi_i - phi_s
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
    PetscReal Etas    = user->Etam;   /* sediment eta (a.k.a. Etam) */
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
    PetscReal thcond, cp, rho, dif_vap, rho_vs, fi, fa, fs, mob_sub;
    ThermalCond(user, ice_c, sed_c, &thcond,  NULL);
    HeatCap    (user, ice_c, sed_c, &cp,      NULL);
    Density    (user, ice_c, sed_c, &rho,     NULL);
    VaporDiffus(user, tem,           &dif_vap, NULL);
    RhoVS_I    (user, tem,           &rho_vs,  NULL);
    Fice       (user, ice_c, sed_c, &fi,      NULL);
    Fair       (user, ice_c, sed_c, &fa,      NULL);
    Fsed       (user, ice_c, sed_c, &fs,      NULL);
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

    /* Two-phase ice prefactor: K2P = 3 M_0 / (eta_i + eta_a). */
    PetscReal K2P    = 3.0 * mob_sub / (Etai + Etaa);
    /* Three-phase prefactor (used only during relax): C3 = 3 M_0 / (eps * Sigma_T). */
    PetscReal SigmaT = Etai + Etas + Etaa;
    PetscReal C3     = 3.0 * mob_sub / (eps * SigmaT);

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

        PetscScalar R_ice, R_sed, R_tem, R_vap;

        if (user->flag_relax) {
            /* ============ Three-phase relaxation (beta-eliminated) ============
             * dphi_i/dt = -C3 [(Sigma_s+Sigma_a) fi - Sigma_s fa - Sigma_a fs]
             *             + 3 M_0 eps grad^2(phi_i) + S_sub
             * dphi_s/dt = -C3 [-Sigma_a fi - Sigma_i fa - (Sigma_i+Sigma_a) fs]
             *             + 3 M_0 eps grad^2(phi_s)
             * rho cp dT/dt = div(K grad T) - rho L_sub d phi_a/dt
             *              = div(K grad T) + rho L_sub (ice_t + sed_t)
             * drho_v/dt = div(D_eff grad rho_v)
             *           - k_pen g(phi_i+phi_s) (rho_v - rho_v_eq)
             *           - rho_ice d phi_i/dt
             * No grad^2(other phase) cross-coupling in the 3-phase form — that
             * only appears after beta-elimination under the 2-phase constraint. */
            R_ice = N0[a] * ice_t
                  + 3.0 * mob_sub * eps * grad_N_dot_grad_ice
                  + C3 * ((Etas + Etaa) * fi - Etas * fa - Etaa * fs) * N0[a]
                  - S_sub * N0[a];

            R_sed = N0[a] * sed_t
                  + 3.0 * mob_sub * eps * grad_N_dot_grad_sed
                  + C3 * (-Etaa * fi - Etai * fa - (Etai + Etaa) * fs) * N0[a];

            R_tem = rho * cp * N0[a] * tem_t
                  + thcond * grad_N_dot_grad_tem
                  - rho * lat_sub * (ice_t + sed_t) * N0[a];

            R_vap = N0[a] * rhov_t
                  + D_eff * grad_N_dot_grad_rhov
                  + vap_pen * N0[a]
                  + rho_ice * ice_t * N0[a];
        } else {
            /* ============ Two-phase (sediment frozen) =========================
             * dphi_i/dt = -K2P [(fi-fa)/eps - (Sigma_i+Sigma_a) eps grad^2(phi_i)
             *                                - Sigma_a grad^2(phi_s)] + S_sub
             * dphi_s/dt = 0
             * rho cp dT/dt = div(K grad T) + rho L_sub S_sub
             * drho_v/dt   = div(D_eff grad rho_v) - k_pen g(...)(rho_v - rho_v_eq)
             *                                     - rho_ice S_sub
             * The grad^2(phi_s) cross-coupling in R_ice is the algebraic trace of
             * the 3-phase derivation under the dphi_s/dt = 0 constraint. */
            R_ice = N0[a] * ice_t
                  + 3.0 * mob_sub * eps * grad_N_dot_grad_ice
                  + K2P * Etaa      * grad_N_dot_grad_sed
                  + (K2P / eps) * (fi - fa) * N0[a]
                  - S_sub * N0[a];

            R_sed = N0[a] * sed_t;

            R_tem = rho * cp * N0[a] * tem_t
                  + thcond * grad_N_dot_grad_tem
                  - rho * lat_sub * S_sub * N0[a];

            R_vap = N0[a] * rhov_t
                  + D_eff * grad_N_dot_grad_rhov
                  + vap_pen * N0[a]
                  + rho_ice * S_sub * N0[a];
        }

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
    PetscReal Etas    = user->Etam;  /* sediment eta */
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

    PetscReal dfi_dice, dfa_dice, dfs_dice;
    Fice(user, ice_c, sed_c, NULL, &dfi_dice);
    Fair(user, ice_c, sed_c, NULL, &dfa_dice);
    Fsed(user, ice_c, sed_c, NULL, &dfs_dice);  /* Lambda cross-term only */

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

    PetscReal K2P    = 3.0 * mob_sub / (Etai + Etaa);    /* used post-relax */
    PetscReal SigmaT = Etai + Etas + Etaa;
    PetscReal C3     = 3.0 * mob_sub / (eps * SigmaT);   /* used during relax */
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

            /* Common ice-equation derivatives that don't change between branches:
             *   S_sub appears in R_ice as -S_sub, and S_sub = alpha_sub*loc*(rhov-rhovs)/rho_ice.
             *   So [ice,T] gets +alpha_sub*loc*d_rho_vs/rho_ice, [ice,rhov] gets -alpha_sub*loc/rho_ice. */
            J[a][0][b][1] += alph_sub * loc * d_rho_vs / rho_ice * Na_Nb;
            J[a][0][b][2] += -alph_sub * loc / rho_ice * Na_Nb;

            if (user->flag_relax) {
                /* ============ Three-phase Jacobian =============================
                 * No grad^2(other phase) cross-coupling in the AC blocks.
                 * R_tem couples to ice_t + sed_t (time-derivative), giving
                 * -rho*lat_sub*shift*N*N in [tem,ice] and [tem,sed].
                 * R_vap couples to ice_t, giving +rho_ice*shift*N*N in [vap,ice]. */
                J[a][0][b][0] += shift * Na_Nb
                               + 3.0 * mob_sub * eps * N1a_N1b
                               + C3 * ((Etas + Etaa) * dfi_dice
                                       - Etas * dfa_dice
                                       - Etaa * dfs_dice) * Na_Nb
                               - alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb;

                J[a][0][b][3] += C3 * ((Etas + Etaa) * dfi_dsed
                                       - Etas * dfa_dsed
                                       - Etaa * dfs_dsed) * Na_Nb
                               - alph_sub * dloc_dsed * rho_v_minus_rho_vs / rho_ice * Na_Nb;

                J[a][3][b][0] += C3 * (-Etaa * dfi_dice
                                       - Etai * dfa_dice
                                       - (Etai + Etaa) * dfs_dice) * Na_Nb;
                J[a][3][b][3] += shift * Na_Nb
                               + 3.0 * mob_sub * eps * N1a_N1b
                               + C3 * (-Etaa * dfi_dsed
                                       - Etai * dfa_dsed
                                       - (Etai + Etaa) * dfs_dsed) * Na_Nb;

                J[a][1][b][0] += S_T * (
                      dthcond_dice * grad_Na_dot_grad_tem * N0[b]
                    - rho * lat_sub * shift * Na_Nb
                );
                J[a][1][b][1] += S_T * (
                      shift * rho * cp * Na_Nb
                    + thcond * N1a_N1b
                );
                /* J[a][1][b][2] = 0 — latent in 3-phase is rho*L_sub*(ice_t+sed_t), no rhov dependence */
                J[a][1][b][3] += S_T * (
                      dthcond_dsed * grad_Na_dot_grad_tem * N0[b]
                    - rho * lat_sub * shift * Na_Nb
                );

                J[a][2][b][0] += dDeff_dice * grad_Na_dot_grad_rhov * N0[b]
                               + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                               - k_pen * (PetscReal)g_solid  * d_rhoveq_dice * Na_Nb
                               + rho_ice * shift * Na_Nb;
                J[a][2][b][1] += dDeff_dtem * grad_Na_dot_grad_rhov * N0[b]
                               - k_pen * (PetscReal)g_solid * d_rhoveq_dtem * Na_Nb;
                J[a][2][b][2] += shift * Na_Nb
                               + D_eff * N1a_N1b
                               + k_pen * (PetscReal)g_solid * (1.0 - d_rhoveq_drv) * Na_Nb;
                J[a][2][b][3] += dDeff_dsed * grad_Na_dot_grad_rhov * N0[b]
                               + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                               - k_pen * (PetscReal)g_solid  * d_rhoveq_dsed * Na_Nb;
            } else {
                /* ============ Two-phase Jacobian (sediment frozen) ============= */
                J[a][0][b][0] += shift * Na_Nb
                               + 3.0 * mob_sub * eps * N1a_N1b
                               + (K2P / eps) * (dfi_dice - dfa_dice) * Na_Nb
                               - alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb;

                J[a][0][b][3] += K2P * Etaa * N1a_N1b
                               + (K2P / eps) * (dfi_dsed - dfa_dsed) * Na_Nb
                               - alph_sub * dloc_dsed * rho_v_minus_rho_vs / rho_ice * Na_Nb;

                /* R_sed = N0*sed_t — only [sed,sed] has shift mass */
                J[a][3][b][3] += shift * Na_Nb;

                J[a][1][b][0] += S_T * (
                      dthcond_dice * grad_Na_dot_grad_tem * N0[b]
                    - rho * lat_sub * alph_sub * dloc_dice * rho_v_minus_rho_vs / rho_ice * Na_Nb
                );
                J[a][1][b][1] += S_T * (
                      shift * rho * cp * Na_Nb
                    + thcond * N1a_N1b
                    + rho * lat_sub * alph_sub * loc * d_rho_vs / rho_ice * Na_Nb
                );
                J[a][1][b][2] += S_T * (
                    - rho * lat_sub * alph_sub * loc / rho_ice * Na_Nb
                );
                J[a][1][b][3] += S_T * (
                      dthcond_dsed * grad_Na_dot_grad_tem * N0[b]
                    - rho * lat_sub * alph_sub * dloc_dsed * rho_v_minus_rho_vs / rho_ice * Na_Nb
                );

                J[a][2][b][0] += dDeff_dice * grad_Na_dot_grad_rhov * N0[b]
                               + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                               - k_pen * (PetscReal)g_solid  * d_rhoveq_dice * Na_Nb
                               + alph_sub * dloc_dice * rho_v_minus_rho_vs * Na_Nb;
                J[a][2][b][1] += dDeff_dtem * grad_Na_dot_grad_rhov * N0[b]
                               - k_pen * (PetscReal)g_solid * d_rhoveq_dtem * Na_Nb
                               - alph_sub * loc * d_rho_vs * Na_Nb;
                J[a][2][b][2] += shift * Na_Nb
                               + D_eff * N1a_N1b
                               + k_pen * (PetscReal)g_solid * (1.0 - d_rhoveq_drv) * Na_Nb
                               + alph_sub * loc * Na_Nb;
                J[a][2][b][3] += dDeff_dsed * grad_Na_dot_grad_rhov * N0[b]
                               + k_pen * (PetscReal)dg_solid * (PetscRealPart(rhov) - rhov_eq) * Na_Nb
                               - k_pen * (PetscReal)g_solid  * d_rhoveq_dsed * Na_Nb
                               + alph_sub * dloc_dsed * rho_v_minus_rho_vs * Na_Nb;
            }
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
