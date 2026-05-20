#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * Helper: derivatives of fi, fs, fa w.r.t. sediment
 * (the material_properties.c functions only provide derivatives w.r.t. ice;
 * these are needed for the Cahn-Hilliard gradient term in Avenue 3)
 * ========================================================================= */
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

    /* fi = Etai*i*(1-i)*(1-2i) + 2L*i*s²*a²
     * ∂fi/∂s = 4L*i*s*a*(a - s) */
    if (dfi_dsed_out)
        *dfi_dsed_out = 4.0 * Lambd * ice * sed * air * (air - sed);

    /* fs = Etam*s*(1-s)*(1-2s) + 2L*i²*s*a²
     * ∂fs/∂s = Etam*(1-6s+6s²) + 2L*i²*a*(a - 2s) */
    if (dfs_dsed_out)
        *dfs_dsed_out = Etased * (1.0 - 6.0*sed + 6.0*sed*sed)
                      + 2.0 * Lambd * ice*ice * air * (air - 2.0*sed);

    /* fa = Etaa*a*(1-a)*(1-2a) + 2L*i²*s²*a
     * ∂fa/∂s = -Etaa*(1-6a+6a²) + 2L*i²*s*(2a - s)   [∂air/∂s = -1] */
    if (dfa_dsed_out)
        *dfa_dsed_out = -Etaa * (1.0 - 6.0*air + 6.0*air*air)
                      + 2.0 * Lambd * ice*ice * sed * (2.0*air - sed);
}


/* =========================================================================
 * Avenue 1: Allen-Cahn + penalty vapor; after t_sed_freeze: sediment RHS = 0
 *
 * THREE-PHASE (flag_sed_frozen = PETSC_FALSE):
 *   Full 3-phase Allen-Cahn for both ice and sediment.
 *   Vapor uses penalised diffusivity (D_pen = difvap_pen * difvap) and
 *   interface equilibrium stiffness (k_pen).
 *
 * TWO-PHASE (flag_sed_frozen = PETSC_TRUE, t >= t_sed_freeze):
 *   Sediment RHS is zeroed (phi_sed is stationary).
 *   Ice uses the 2-phase Allen-Cahn form pinned to the frozen sediment boundary.
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
    PetscReal Etased  = user->Etam;
    PetscReal EtaT    = Etai*Etased + Etai*Etaa + Etased*Etaa;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;
    PetscReal xi_v    = user->xi_v;
    PetscReal xi_T    = user->xi_T;
    PetscReal alph_sub = user->alph_sub;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    /* Hessian of all DOFs (needed for Gibbs-Thomson curvature of phi_ice).
     * IGAPointFormHess fills hess_sol[dof][dim][dim] from the second
     * derivatives of the basis at this quadrature point. Layout matches
     * grad_sol: dof-major, then dim-major. With p=2 C^1 B-splines the
     * Hessian exists in the interior of each element; quadrature points
     * are interior so this is well-defined. */
    PetscScalar hess_sol[4][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);

    PetscScalar ice = sol[0], ice_t = sol_t[0];
    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];
    PetscScalar hess_ice[dim * dim];
    {
        PetscInt _k, _l;
        for (_k = 0; _k < dim; _k++)
            for (_l = 0; _l < dim; _l++)
                hess_ice[_k * dim + _l] = hess_sol[0][_k][_l];
    }

    PetscScalar sed = sol[3], sed_t = sol_t[3];
    PetscScalar grad_sed[dim];
    for (l = 0; l < dim; l++) grad_sed[l] = grad_sol[3][l];

    PetscScalar air   = 1.0 - sed - ice;
    PetscScalar air_t = -ice_t - sed_t;

    /* SNES domain check: if the *trial* Newton iterate has any phi outside
     * [phase_lo, phase_hi], tell SNES the current line-search trial is in an
     * invalid region. Line search will treat the residual norm as infinite,
     * back off the step, and try a smaller update. Catches the dt-induced AC
     * instability while Newton is still iterating — earlier and cheaper than
     * the deferred-rollback safety net in Monitor(), which only fires after
     * a step has been accepted.
     *
     * The check uses the *unclamped* values read straight from the trial
     * iterate; the clamps below are only for downstream material-property
     * arithmetic (so divisions by phi values don't NaN). */
    {
        PetscReal phi_i_un = PetscRealPart(sol[0]);
        PetscReal phi_s_un = PetscRealPart(sol[3]);
        PetscReal phi_a_un = 1.0 - phi_i_un - phi_s_un;
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (phi_i_un < lo || phi_i_un > hi ||
            phi_s_un < lo || phi_s_un > hi ||
            phi_a_un < lo || phi_a_un > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    if (ice < 0.0) ice = 0.0;
    if (ice > 1.0) ice = 1.0;
    if (sed < 0.0) sed = 0.0;
    if (sed > 1.0) sed = 1.0;
    if (air < 0.0) air = 0.0;
    if (air > 1.0) air = 1.0;

    PetscScalar tem  = sol[1], tem_t = sol_t[1];
    PetscScalar grad_tem[dim];
    for (l = 0; l < dim; l++) grad_tem[l] = grad_sol[1][l];

    PetscScalar rhov = sol[2], rhov_t = sol_t[2];
    PetscScalar grad_rhov[dim];
    for (l = 0; l < dim; l++) grad_rhov[l] = grad_sol[2][l];

    PetscReal thcond, cp, rho, difvap, rhoI_vs, fi, fa, fs, mob;
    ThermalCond(user, ice, sed, &thcond,  NULL);
    HeatCap    (user, ice, sed, &cp,      NULL);
    Density    (user, ice, sed, &rho,     NULL);
    VaporDiffus(user, tem,      &difvap,  NULL);
    RhoVS_I    (user, tem,      &rhoI_vs, NULL);
    Fice(user, ice, sed, &fi, NULL);
    Fair(user, ice, sed, &fa, NULL);
    Fsed(user, ice, sed, &fs, NULL);
    Mobility(user, ice, sed, &mob);

    /* D_pen = difvap_pen (factor) * difvap: penalised diffusivity in the
     * SOLID phases (ice + sed). Vapor diffuses at the full physical rate in
     * the air phase, and is throttled by difvap_pen wherever phi_ice + phi_sed
     * is large. Without this, vapor in the solid interior would drift
     * unboundedly via vap_src accumulated by tiny numerical ice motions. */
    PetscReal D_pen = user->difvap_pen * difvap;
    PetscReal k_pen = user->k_pen;

    PetscReal g_solid, g_pen;
    PetscReal rhov_eq = ice * rhoI_vs + sed * rhoI_vs + air * rhov;
    /* g_pen: penalty weight concentrated near ice+sed = 1 (deep solid).
     * Zero at the diffuse ice-air interface so Gibbs-Thomson can emerge there;
     * unity in the solid interior so rhov stays pinned to rhov_eq. */
    PenaltyWeight     (ice + sed, &g_pen,  NULL);
    /* g_solid = g(ice + sed): smooth Heaviside that is 0 in pure air and 1 in
     * either pure ice or pure sed. Used to select between full vapor
     * diffusivity (in air) and the penalised D_pen (in solid). */
    SmoothHeavisidePoly(ice + sed, &g_solid, NULL);
    difvap = D_pen * g_solid + difvap * (1.0 - g_solid);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar air_eff = (air > air_lim) ? air : air_lim;
    PetscReal C3  = 3.0 * mob / (eps * EtaT);
    PetscReal loc = ice*ice * air*air;

    /* Gibbs-Thomson curvature correction to the equilibrium vapor density.
     * rhoI_vs_eff = rhoI_vs * (1 + d0_GT * kappa), where kappa = -div(grad_ice/|grad_ice|).
     * The localizer ice^2 * air^2 in sub_src restricts the correction to the
     * diffuse interface, so the (numerically uninteresting) bulk value of
     * kappa doesn't matter physically. Regularization scale 0.01/eps keeps
     * kappa bounded in bulk where |grad_ice| -> 0.
     * Default d0_GT = 0 disables the correction (rhoI_vs_eff = rhoI_vs). */
    PetscScalar kappa = 0.0;
    if (user->d0_GT != 0.0) {
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps, &kappa, NULL, NULL);
    }
    PetscReal rhoI_vs_eff  = rhoI_vs * (1.0 + user->d0_GT * (PetscReal)kappa);

    /* Hoist loop-invariant scalars out of the per-node loop */
    PetscReal AC_ice_bulk  = C3 * ((Etased + Etaa)*fi - Etaa*fs - Etased*fa);
    PetscReal AC_sed_bulk  = C3 * (-Etaa*fi - Etai*fa + (Etai + Etaa)*fs);
    PetscReal sub_src      = alph_sub * loc * (rhov - rhoI_vs_eff) / rho_ice;
    PetscReal rho_lat_air_t = xi_T * rho * lat_sub * air_t;
    PetscReal vap_pen      = xi_v * k_pen * g_pen * (rhov - rhov_eq);
    /* Mass-exchange source: Stefan condition in the diffuse-interface limit.
     * Couples vapor to *ice* motion only, not to sed motion. Sediment doesn't
     * sublimate, so sed_t must not generate or consume vapor; using
     *     vap_src = ρ_ice · air_t = -ρ_ice·(ice_t + sed_t)
     * would (wrongly) treat sed-AC interface motion as a sublimation event.
     * The form below is identical in the post-pin regime (sed_t = 0) and
     * eliminates the spurious vapor source at the sed boundary during the
     * pre-pin AC relaxation.
     * No xi_v scaling — this is the physical mass-balance closure; xi_v
     * scales only the regularisation terms (diffusion + equilibrium penalty). */
    PetscReal vap_src      = -rho_ice * ice_t;

    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal R_ice = 0.0, R_sed = 0.0, R_tem = 0.0, R_vap = 0.0;

        if (user->flag_relax) {

            /* ================================================================
             * RELAXATION MODE  (step < n_relax)
             * Run 3-phase Allen-Cahn for ice and sediment — no sublimation,
             * no thermal or vapor spatial coupling.
             * T and rhov carry only their time-derivative terms so ∂T/∂t = 0
             * and ∂rhov/∂t = 0, keeping the Jacobian non-singular while
             * holding both fields stationary at IC values.
             * ================================================================ */

            /* Ice: 3-phase AC, no sublimation */
            R_ice = N0[a] * ice_t;
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
            R_ice += AC_ice_bulk * N0[a];

            /* Sediment: 3-phase AC */
            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += AC_sed_bulk * N0[a];

            /* Temperature and vapor: time-derivative only (forces ∂T/∂t = ∂rhov/∂t = 0)
             * This keeps the Jacobian non-singular while the fields stay fixed. */
            R_tem = rho * cp * N0[a] * tem_t;
            R_vap = N0[a] * rhov_t;

        } else if (!user->flag_sed_frozen) {

            /* ================================================================
             * THREE-PHASE FORMULATION
             * Active from t = 0 until t = t_sed_freeze (> 0).
             * Both ice and sediment evolve under full 3-phase Allen-Cahn.
             * ================================================================ */

            /* Ice: 3-phase Allen-Cahn */
            R_ice = N0[a] * ice_t;
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
            R_ice += AC_ice_bulk * N0[a];
            R_ice -= sub_src * N0[a];

            /* Sediment: 3-phase Allen-Cahn */
            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += AC_sed_bulk * N0[a];

        } else {

            /* ================================================================
             * FROZEN-SEDIMENT FORMULATION  (t >= t_sed_freeze)
             *
             * Sediment is held rigid: phi_sed_t = 0. The ice equation, however,
             * keeps the full 3-phase Allen-Cahn form. The three-phase potential
             * is symmetric in the (ice, sed, air) phases, so ice = 0 remains a
             * stable equilibrium at a sed-air interface even when sed is frozen.
             *
             * (The prior 2-phase Kim-Steinbach form used here was derived under
             *  the assumption that the frozen sed boundary always neighbours
             *  ice; at sed-air boundaries it produced large spurious ice. The
             *  3-phase ice eq above does not have that failure mode.)
             * ================================================================ */

            /* Ice: 3-phase Allen-Cahn (identical to the 3-phase branch) */
            R_ice = N0[a] * ice_t;
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
            R_ice += AC_ice_bulk * N0[a];
            R_ice -= sub_src * N0[a];

            /* Sediment: RHS = 0 (time derivative only — forces phi_sed_t = 0) */
            R_sed = N0[a] * sed_t;
        }

        if (!user->flag_relax) {
            /* --- Thermal energy balance (same in both formulations) --- */
            R_tem  = rho * cp * N0[a] * tem_t;
            for (l = 0; l < dim; l++)
                R_tem += xi_T * thcond * (N1[a][l] * grad_tem[l]);
            R_tem += rho_lat_air_t * N0[a];

            /* --- Vapor: penalised diffusivity + interface equilibrium penalty
             *     (same in both formulations) --- */
            R_vap  = N0[a] * rhov_t;
            for (l = 0; l < dim; l++)
                R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
            R_vap += vap_pen * N0[a];
            R_vap -= vap_src * N0[a];
        }

        R[a][0] = R_ice;
        R[a][1] = R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;
    }
    return 0;
}


/* =========================================================================
 * Avenue 2: disabled — retained for future use. Use -flag_avenue 1.
 * ========================================================================= */
PetscErrorCode Residual_A2(IGAPoint pnt,
                            PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U,
                            PetscScalar *Re, void *ctx)
{
    (void)pnt; (void)shift; (void)V; (void)t; (void)U; (void)Re; (void)ctx;
    return 0;
}


/* =========================================================================
 * Avenue 3: disabled — retained for future use. Use -flag_avenue 1.
 * ========================================================================= */
PetscErrorCode Residual_A3(IGAPoint pnt,
                            PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U,
                            PetscScalar *Re, void *ctx)
{
    (void)pnt; (void)shift; (void)V; (void)t; (void)U; (void)Re; (void)ctx;
    return 0;
}


/* =========================================================================
 * Dispatcher: routes to the avenue selected by user->flag_avenue.
 *   1 = Avenue 1 (AC, penalty vapor, freeze → zero sed RHS)  [active]
 *   2 = Avenue 2 (disabled — retained for future use)
 *   3 = Avenue 3 (disabled — retained for future use)
 * ========================================================================= */
PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;
    switch (user->flag_avenue) {
        case 1:  return Residual_A1(pnt, shift, V, t, U, Re, ctx);
        case 2:  /* fallthrough — disabled */
        case 3:  /* fallthrough — disabled */
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Avenues 2 and 3 are disabled. Use -flag_avenue 1.");
        default:
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Unknown -flag_avenue %d (valid: 1)", user->flag_avenue);
    }
}


/* =========================================================================
 * Jacobian_A1 — analytical element Jacobian for Avenue 1
 *
 * J[a][i][b][j] = ∂R[a][i]/∂u[b][j]  +  shift × ∂R[a][i]/∂u_t[b][j]
 * DOF layout: [0]=ice, [1]=tem, [2]=rhov, [3]=sed
 *
 * Notes:
 *  - Mobility (mob) is treated as frozen at the current iterate; its
 *    derivatives through C2/C3 are omitted (modified-Newton for those terms).
 *  - All material-property derivatives are included.
 *  - N1[a]·grad_phi terms that appear in the [vap,*] and [tem,ice] off-diagonal
 *    blocks are computed via a per-a dot-product precomputation.
 * ========================================================================= */
static PetscErrorCode Jacobian_A1(IGAPoint pnt,
                                   PetscReal shift, const PetscScalar *V,
                                   PetscReal t, const PetscScalar *U,
                                   PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    PetscReal eps      = user->eps;
    PetscReal Etai     = user->Etai;
    PetscReal Etaa     = user->Etaa;
    PetscReal Etased   = user->Etam;
    PetscReal EtaT     = Etai*Etased + Etai*Etaa + Etased*Etaa;
    PetscReal rho_ice  = user->rho_ice;
    PetscReal lat_sub  = user->lat_sub;
    PetscReal air_lim  = user->air_lim;
    PetscReal xi_v     = user->xi_v;
    PetscReal xi_T     = user->xi_T;
    PetscReal alph_sub = user->alph_sub;

    if (pnt->atboundary) return 0;

    (void)V;  /* time derivatives not needed in the Jacobian */

    PetscScalar sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    /* Hessian of ice for the GT-curvature derivatives. Needed only when
     * user->d0_GT != 0; cheap to read unconditionally and skip the
     * derivative work below if disabled. */
    PetscScalar hess_sol[4][dim][dim];
    IGAPointFormHess(pnt, U, &hess_sol[0][0][0]);

    PetscScalar ice  = sol[0];
    PetscScalar sed  = sol[3];
    PetscScalar tem  = sol[1];
    PetscScalar rhov = sol[2];

    PetscScalar grad_tem[dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) grad_tem[l]  = grad_sol[1][l];
    for (l = 0; l < dim; l++) grad_rhov[l] = grad_sol[2][l];

    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];
    PetscScalar hess_ice[dim * dim];
    {
        PetscInt _k, _l;
        for (_k = 0; _k < dim; _k++)
            for (_l = 0; _l < dim; _l++)
                hess_ice[_k * dim + _l] = hess_sol[0][_k][_l];
    }

    PetscScalar air = 1.0 - ice - sed;

    if (ice < 0.0) ice = 0.0;
    if (ice > 1.0) ice = 1.0;
    if (sed < 0.0) sed = 0.0;
    if (sed > 1.0) sed = 1.0;
    if (air < 0.0) air = 0.0;
    if (air > 1.0) air = 1.0;

    /* Material properties (values needed for Jacobian coefficients) */
    PetscReal thcond, cp, rho, difvap_raw, rhoI_vs, fi, fa, fs, mob;
    ThermalCond(user, ice, sed, &thcond,    NULL);
    HeatCap    (user, ice, sed, &cp,        NULL);
    Density    (user, ice, sed, &rho,       NULL);
    VaporDiffus(user, tem,      &difvap_raw, NULL);
    RhoVS_I    (user, tem,      &rhoI_vs,   NULL);
    Fice(user, ice, sed, &fi, NULL);
    Fair(user, ice, sed, &fa, NULL);
    Fsed(user, ice, sed, &fs, NULL);
    Mobility(user, ice, sed, &mob);

    /* Effective (penalised) diffusivity — same formula as in Residual_A1.
     * Penalty applies in solid (ice + sed); air retains full physical
     * vapor diffusivity. */
    PetscReal D_pen = user->difvap_pen * difvap_raw;
    PetscReal k_pen = user->k_pen;
    PetscReal g_solid, g_pen;
    SmoothHeavisidePoly(ice + sed, &g_solid, NULL);
    PenaltyWeight     (ice + sed, &g_pen,   NULL);
    PetscReal difvap = D_pen * g_solid + difvap_raw * (1.0 - g_solid);

    PetscReal rhov_eq  = ice * rhoI_vs + sed * rhoI_vs + air * rhov;
    PetscReal air_eff  = (air > air_lim) ? air : air_lim;
    PetscInt  air_eff_active = (air > air_lim);

    PetscReal C3  = 3.0 * mob / (eps * EtaT);
    PetscReal loc = ice*ice * air*air;

    /* ---- Derivatives of material functions ---- */
    PetscReal dfi_dice, dfs_dice, dfa_dice;
    Fice(user, ice, sed, NULL, &dfi_dice);
    Fsed(user, ice, sed, NULL, &dfs_dice);
    Fair(user, ice, sed, NULL, &dfa_dice);

    PetscReal dfi_dsed, dfs_dsed, dfa_dsed;
    ChemPot_dsed(user, ice, sed, &dfi_dsed, &dfs_dsed, &dfa_dsed);

    PetscReal dcond_ice, d_difvap_raw, d_rhovs;
    ThermalCond(user, ice, sed, NULL, &dcond_ice);
    VaporDiffus(user, tem,      NULL, &d_difvap_raw);
    RhoVS_I    (user, tem,      NULL, &d_rhovs);

    PetscReal dg_solid, dg_pen;
    SmoothHeavisidePoly(ice + sed, NULL, &dg_solid);
    PenaltyWeight     (ice + sed, NULL, &dg_pen);

    /* Gibbs-Thomson curvature correction — see Residual_A1 for the physics.
     * We also need d kappa / d (grad_ice) and d kappa / d (hess_ice) here for
     * the chain rule into the J[ice, ice] block. */
    PetscScalar kappa = 0.0;
    PetscScalar dkappa_dg[3] = {0.0, 0.0, 0.0};
    PetscScalar dkappa_dH[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    if (user->d0_GT != 0.0) {
        Curvature(dim, grad_ice, hess_ice, 0.01 / eps,
                  &kappa, dkappa_dg, dkappa_dH);
    }
    PetscReal rhoI_vs_eff      = rhoI_vs * (1.0 + user->d0_GT * (PetscReal)kappa);
    PetscReal d_rhovs_eff_dtem = d_rhovs * (1.0 + user->d0_GT * (PetscReal)kappa);

    /* d(difvap_eff)/d(ice) and d(difvap_eff)/d(sed) via the smoothed switch
     * on (ice + sed). Both are equal: g_solid = g(ice + sed), so
     * d g_solid / d ice = d g_solid / d sed = dg_solid. Formula:
     *   difvap_eff = D_pen * g_solid + difvap_raw * (1 - g_solid)
     *   d/dphi    = (D_pen - difvap_raw) * dg_solid                  */
    PetscReal d_difvap_eff_dice = (D_pen - difvap_raw) * dg_solid;
    PetscReal d_difvap_eff_dsed = d_difvap_eff_dice;
    /* d(difvap_eff)/d(tem). D_pen = difvap_pen * difvap_raw, both branches
     * scale linearly with d_difvap_raw. */
    PetscReal d_difvap_eff_dtem = d_difvap_raw * (user->difvap_pen * g_solid
                                                  + (1.0 - g_solid));

    /* d(loc)/d(ice) = 2*ice*air*(air - ice),  d(loc)/d(sed) = -2*ice^2*air */
    PetscReal dloc_dice = 2.0*ice*air*(air - ice);
    PetscReal dloc_dsed = -2.0*ice*ice*air;

    /* d(rhov_eq)/d(x): rhov_eq = (ice+sed)*rhoI_vs + air*rhov */
    PetscReal d_rhovel_dice  = rhoI_vs - rhov;
    PetscReal d_rhovel_dsed  = rhoI_vs - rhov;
    PetscReal d_rhovel_drhov = air;
    PetscReal d_rhovel_dtem  = (ice + sed) * d_rhovs;

    const PetscReal *N0, (*N1)[dim], (*N2)[dim][dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    /* N2 (second-derivative shape funs) is only needed for the Gibbs-Thomson
     * curvature Jacobian; reading unconditionally is cheap. */
    IGAPointGetShapeFuns(pnt, 2, (const PetscReal**)&N2);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;

    for (a = 0; a < nen; a++) {
        /* Per-a dot products with field gradients (for asymmetric off-diagonal terms) */
        PetscReal N1a_grad_tem  = 0.0;
        PetscReal N1a_grad_rhov = 0.0;
        for (l = 0; l < dim; l++) {
            N1a_grad_tem  += N1[a][l] * grad_tem[l];
            N1a_grad_rhov += N1[a][l] * grad_rhov[l];
        }

        for (b = 0; b < nen; b++) {
            PetscReal Na_Nb   = N0[a] * N0[b];
            PetscReal N1a_N1b = 0.0;
            for (l = 0; l < dim; l++) N1a_N1b += N1[a][l] * N1[b][l];

            if (user->flag_relax) {
                /* ==============================================================
                 * RELAXATION: 3-phase AC for ice/sed; mass matrices for T/rhov
                 * ============================================================== */
                /* [ice, ice] */
                J[a][0][b][0] += shift * Na_Nb
                               + 3.0*mob*eps * N1a_N1b
                               + C3 * ((Etased+Etaa)*dfi_dice - Etaa*dfs_dice
                                       - Etased*dfa_dice) * Na_Nb;
                /* [ice, sed] */
                J[a][0][b][3] += C3 * ((Etased+Etaa)*dfi_dsed - Etaa*dfs_dsed
                                        - Etased*dfa_dsed) * Na_Nb;
                /* [tem, tem] */
                J[a][1][b][1] += shift * rho * cp * Na_Nb;
                /* [vap, vap] */
                J[a][2][b][2] += shift * Na_Nb;
                /* [sed, ice] */
                J[a][3][b][0] += C3 * (-Etaa*dfi_dice - Etai*dfa_dice
                                        + (Etai+Etaa)*dfs_dice) * Na_Nb;
                /* [sed, sed] */
                J[a][3][b][3] += shift * Na_Nb
                               + 3.0*mob*eps * N1a_N1b
                               + C3 * (-Etaa*dfi_dsed - Etai*dfa_dsed
                                        + (Etai+Etaa)*dfs_dsed) * Na_Nb;

            } else if (!user->flag_sed_frozen) {
                /* ==============================================================
                 * THREE-PHASE (before t_sed_freeze)
                 * ============================================================== */
                /* [ice, ice] */
                J[a][0][b][0] += shift * Na_Nb
                               + 3.0*mob*eps * N1a_N1b
                               + C3 * ((Etased+Etaa)*dfi_dice - Etaa*dfs_dice
                                       - Etased*dfa_dice) * Na_Nb
                               - alph_sub * dloc_dice * (rhov - rhoI_vs_eff) / rho_ice * Na_Nb;
                /* [ice, tem]: d(rhoI_vs_eff)/d(tem) in sublimation term */
                J[a][0][b][1] += alph_sub * loc * d_rhovs_eff_dtem / rho_ice * Na_Nb;
                /* [ice, rhov] */
                J[a][0][b][2] += -alph_sub * loc / rho_ice * Na_Nb;
                /* [ice, sed] */
                J[a][0][b][3] += C3 * ((Etased+Etaa)*dfi_dsed - Etaa*dfs_dsed
                                        - Etased*dfa_dsed) * Na_Nb
                               - alph_sub * dloc_dsed * (rhov - rhoI_vs_eff) / rho_ice * Na_Nb;
                /* [sed, ice] */
                J[a][3][b][0] += C3 * (-Etaa*dfi_dice - Etai*dfa_dice
                                        + (Etai+Etaa)*dfs_dice) * Na_Nb;
                /* [sed, sed] */
                J[a][3][b][3] += shift * Na_Nb
                               + 3.0*mob*eps * N1a_N1b
                               + C3 * (-Etaa*dfi_dsed - Etai*dfa_dsed
                                        + (Etai+Etaa)*dfs_dsed) * Na_Nb;

            } else {
                /* ==============================================================
                 * FROZEN SEDIMENT (t >= t_sed_freeze)
                 * Ice rows match the 3-phase branch exactly; only the sediment
                 * row degenerates to the mass matrix that pins sed_t = 0.
                 * ============================================================== */
                /* [ice, ice] */
                J[a][0][b][0] += shift * Na_Nb
                               + 3.0*mob*eps * N1a_N1b
                               + C3 * ((Etased+Etaa)*dfi_dice - Etaa*dfs_dice
                                       - Etased*dfa_dice) * Na_Nb
                               - alph_sub * dloc_dice * (rhov - rhoI_vs_eff) / rho_ice * Na_Nb;
                /* [ice, tem] */
                J[a][0][b][1] += alph_sub * loc * d_rhovs_eff_dtem / rho_ice * Na_Nb;
                /* [ice, rhov] */
                J[a][0][b][2] += -alph_sub * loc / rho_ice * Na_Nb;
                /* [ice, sed] */
                J[a][0][b][3] += C3 * ((Etased+Etaa)*dfi_dsed - Etaa*dfs_dsed
                                        - Etased*dfa_dsed) * Na_Nb
                               - alph_sub * dloc_dsed * (rhov - rhoI_vs_eff) / rho_ice * Na_Nb;
                /* [sed, sed]: mass matrix only — frozen sediment */
                J[a][3][b][3] += shift * Na_Nb;
            }

            /* ==============================================================
             * Gibbs-Thomson correction to sub_src introduces a new dependence
             * of sub_src on grad(phi_ice) and hess(phi_ice). Chain rule:
             *   d sub_src / d phi_ice^b = ... (terms above)
             *                          + alph_sub * loc * (-rhoI_vs * d0_GT)
             *                            * (sum_l dkappa/dg_l * N1[b][l]
             *                              + sum_kl dkappa/dH_{kl} * N2[b][k][l])
             *                            / rho_ice
             * Since R_ice -= sub_src * N0[a], the contribution to J[ice, ice]
             * is the negative of d sub_src / d phi_ice^b times N0[a]:
             *   J[a][0][b][0] += alph_sub * loc * rhoI_vs * d0_GT / rho_ice
             *                  * N0[a] * (grad and hess parts).
             * Active only when d0_GT != 0 and outside the relax window
             * (sub_src isn't in the residual during relax). */
            if (user->d0_GT != 0.0 && !user->flag_relax) {
                PetscReal grad_part = 0.0;
                for (l = 0; l < dim; l++)
                    grad_part += (PetscReal)dkappa_dg[l] * N1[b][l];
                PetscReal hess_part = 0.0;
                PetscInt _k, _l;
                for (_k = 0; _k < dim; _k++)
                    for (_l = 0; _l < dim; _l++)
                        hess_part += (PetscReal)dkappa_dH[_k * dim + _l] * N2[b][_k][_l];
                J[a][0][b][0] += alph_sub * loc * rhoI_vs * user->d0_GT
                              / rho_ice * N0[a] * (grad_part + hess_part);
            }

            if (!user->flag_relax) {
                /* ==============================================================
                 * THERMAL (both 2-phase and 3-phase)
                 * ============================================================== */
                /* [tem, ice]: shift from air_t = -ice_t; spatial from thcond(ice) */
                J[a][1][b][0] += shift * (-xi_T * rho * lat_sub) * Na_Nb
                               + xi_T * dcond_ice * N1a_grad_tem * N0[b];
                /* [tem, tem]: shift×rho*cp + thcond stiffness */
                J[a][1][b][1] += shift * rho * cp * Na_Nb
                               + xi_T * thcond * N1a_N1b;
                /* [tem, sed]: shift from air_t = -sed_t */
                J[a][1][b][3] += shift * (-xi_T * rho * lat_sub) * Na_Nb;

                /* ==============================================================
                 * VAPOR (both 2-phase and 3-phase)
                 * ============================================================== */
                /* [vap, ice]: shift from rho_ice*air_t (mass-balance source, no xi_v);
                 *             penalty + diffusivity variation (regularisation, with xi_v) */
                J[a][2][b][0] += shift * rho_ice * Na_Nb
                               + xi_v * k_pen * (dg_pen * (rhov - rhov_eq)
                                                - g_pen * d_rhovel_dice) * Na_Nb
                               + xi_v * (d_difvap_eff_dice * air_eff
                                        - (air_eff_active ? difvap : 0.0))
                                        * N1a_grad_rhov * N0[b];
                /* [vap, tem]: from difvap(T) and rhoI_vs(T) in rhov_eq */
                J[a][2][b][1] += xi_v * k_pen * g_pen * (-d_rhovel_dtem) * Na_Nb
                               + xi_v * d_difvap_eff_dtem * air_eff * N1a_grad_rhov * N0[b];
                /* [vap, rhov]: shift + diffusion stiffness + penalty */
                J[a][2][b][2] += shift * Na_Nb
                               + xi_v * difvap * air_eff * N1a_N1b
                               + xi_v * k_pen * g_pen * (1.0 - d_rhovel_drhov) * Na_Nb;
                /* [vap, sed]: vap_src = -ρ_ice * ice_t no longer depends on sed_t,
                 *             so the shift contribution drops. Three remaining
                 *             contributions: (1) penalty derivative via dg_pen
                 *             and rhov_eq's sed dependence; (2) difvap_eff
                 *             depends on sed through g_solid = g(ice+sed);
                 *             (3) air = 1 - ice - sed in the air_eff factor. */
                J[a][2][b][3] += xi_v * k_pen * (dg_pen * (rhov - rhov_eq)
                                                - g_pen * d_rhovel_dsed) * Na_Nb
                               + xi_v * (d_difvap_eff_dsed * air_eff
                                        - (air_eff_active ? difvap : 0.0))
                                        * N1a_grad_rhov * N0[b];
            }
        }
    }
    return 0;
}


/* =========================================================================
 * Jacobian dispatcher — routes to Jacobian_A1 for flag_avenue = 1.
 * ========================================================================= */
PetscErrorCode Jacobian(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;
    switch (user->flag_avenue) {
        case 1:  return Jacobian_A1(pnt, shift, V, t, U, Je, ctx);
        default:
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Analytical Jacobian only implemented for -flag_avenue 1");
    }
}


/* =========================================================================
 * Integration — computes per-element scalar integrals for monitoring
 * ========================================================================= */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                            PetscScalar *S, void *ctx)
{
    PetscFunctionBegin;
    PetscScalar sol[4];
    IGAPointFormValue(pnt, U, &sol[0]);

    PetscReal ice  = sol[0];
    PetscReal sed  = sol[3];
    PetscReal air  = 1.0 - sed - ice;
    PetscReal temp = sol[1];
    PetscReal rhov = sol[2];

    S[0] = ice;
    S[1] = SQ(air) * SQ(sed) * SQ(ice);
    S[2] = air;
    S[3] = temp;
    S[4] = rhov * air;
    S[5] = air*air * ice*ice;
    S[6] = sed;
    S[7] = sed*sed * air*air;   /* sed-air interface */
    S[8] = sed*sed * ice*ice;   /* ice-sed interface */

    PetscFunctionReturn(0);
}
