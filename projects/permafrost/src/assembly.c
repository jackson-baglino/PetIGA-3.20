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

    PetscScalar ice = sol[0], ice_t = sol_t[0];
    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];

    PetscScalar sed = sol[3], sed_t = sol_t[3];
    PetscScalar grad_sed[dim];
    for (l = 0; l < dim; l++) grad_sed[l] = grad_sol[3][l];

    PetscScalar air   = 1.0 - sed - ice;
    PetscScalar air_t = -ice_t - sed_t;

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

    /* D_pen = difvap_pen (factor) * difvap: penalised diffusivity in the air phase */
    PetscReal D_pen = user->difvap_pen * difvap;
    PetscReal k_pen = user->k_pen;

    PetscReal g_phia, g_phiiphis;
    PetscReal rhov_eq = ice * rhoI_vs + sed * rhoI_vs + air * rhov;
    SmoothHeavisidePoly(ice + sed, &g_phiiphis, NULL);
    SmoothHeavisidePoly(ice,       &g_phia,     NULL);
    difvap = difvap * g_phia + D_pen * (1.0 - g_phia);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar air_eff = (air > air_lim) ? air : air_lim;
    PetscReal C3  = 3.0 * mob / (eps * EtaT);
    PetscReal loc = ice*ice * air*air;

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
            R_ice += C3 * ((Etased + Etaa)*fi - Etaa*fs - Etased*fa) * N0[a];

            /* Sediment: 3-phase AC */
            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += C3 * (-Etaa*fi - Etai*fa + (Etai + Etaa)*fs) * N0[a];

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
            R_ice += C3 * ((Etased + Etaa)*fi - Etaa*fs - Etased*fa) * N0[a];
            R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;

            /* Sediment: 3-phase Allen-Cahn */
            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += C3 * (-Etaa*fi - Etai*fa + (Etai + Etaa)*fs) * N0[a];

        } else {

            /* ================================================================
             * TWO-PHASE FORMULATION  (sediment frozen, t >= t_sed_freeze)
             * Sediment RHS is zeroed: only the time-derivative term remains,
             * which forces phi_sed_t = 0 (sediment stationary).
             * Ice evolves under 2-phase AC pinned to the frozen sediment
             * boundary.
             * ================================================================ */

            /* Ice: 2-phase Allen-Cahn pinned to the frozen sediment boundary */
            PetscReal C2 = 3.0 * mob / (Etai + Etaa);
            R_ice = N0[a] * ice_t;
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
            for (l = 0; l < dim; l++)
                R_ice += C2 * Etaa * eps * (N1[a][l] * grad_sed[l]);
            R_ice += N0[a] * C2 / eps * (fi - fa);
            R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;

            /* Sediment: RHS = 0 (time derivative only — forces phi_sed_t = 0) */
            R_sed = N0[a] * sed_t;
        }

        if (!user->flag_relax) {
            /* --- Thermal energy balance (same in both formulations) --- */
            R_tem  = rho * cp * N0[a] * tem_t;
            for (l = 0; l < dim; l++)
                R_tem += xi_T * thcond * (N1[a][l] * grad_tem[l]);
            R_tem += xi_T * rho * lat_sub * N0[a] * air_t;

            /* --- Vapor: penalised diffusivity + interface equilibrium penalty
             *     (same in both formulations) --- */
            R_vap  = N0[a] * rhov_t;
            for (l = 0; l < dim; l++)
                R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
            R_vap += k_pen * g_phiiphis * (rhov - rhov_eq) * N0[a];
            R_vap -= xi_v * N0[a] * rho_ice * air_t;
        }

        R[a][0] = R_ice;
        R[a][1] = R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;
    }
    return 0;
}


/* =========================================================================
 * Avenue 2: Allen-Cahn + penalty vapor; after t_sed_freeze: sediment penalty
 *
 * Before t_sed_freeze: identical to Avenue 1.
 * After t_sed_freeze (flag_sed_frozen = 1): sediment continues to evolve via
 *   full Allen-Cahn but a restoring penalty k_sed*(phi_s - phi_s0) is added,
 *   discouraging the sediment from moving away from its relaxed reference state.
 *   Set flag_2ph_ice = 1 to optionally switch ice to the 2-phase form.
 * ========================================================================= */
PetscErrorCode Residual_A2(IGAPoint pnt,
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

    PetscScalar ice = sol[0], ice_t = sol_t[0];
    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];

    PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
    PetscScalar sed = sol[3], sed_t = sol_t[3];
    PetscScalar grad_sed[dim];
    for (l = 0; l < dim; l++) grad_sed[l] = grad_sol[3][l];
    PetscScalar sed0 = user->Phi_sed0[indGP];

    PetscScalar air   = 1.0 - sed - ice;
    PetscScalar air_t = -ice_t - sed_t;

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

    /* D_pen = difvap_pen (factor) * difvap: penalised diffusivity in the air phase */
    PetscReal D_pen  = user->difvap_pen * difvap;
    PetscReal k_pen  = user->k_pen;
    PetscReal k_sed  = user->k_sed_pen;

    PetscReal g_phia, g_phiiphis;
    PetscReal rhov_eq = ice * rhoI_vs + sed * rhoI_vs + air * rhov;
    SmoothHeavisidePoly(ice + sed, &g_phiiphis, NULL);
    SmoothHeavisidePoly(ice,       &g_phia,     NULL);
    difvap = difvap * g_phia + D_pen * (1.0 - g_phia);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar air_eff = (air > air_lim) ? air : air_lim;
    PetscReal C3  = 3.0 * mob / (eps * EtaT);
    PetscReal loc = ice*ice * air*air;

    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal R_ice = 0.0, R_sed = 0.0, R_tem = 0.0, R_vap = 0.0;

        if (user->flag_relax) {

            /* ================================================================
             * RELAXATION MODE  (step < n_relax)
             * 3-phase AC for ice and sediment; no sublimation.
             * T and rhov carry time-derivative only to keep Jacobian regular.
             * ================================================================ */
            R_ice = N0[a] * ice_t;
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
            R_ice += C3 * ((Etased + Etaa)*fi - Etaa*fs - Etased*fa) * N0[a];

            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += C3 * (-Etaa*fi - Etai*fa + (Etai + Etaa)*fs) * N0[a];

            R_tem = rho * cp * N0[a] * tem_t;
            R_vap = N0[a] * rhov_t;

        } else {

            /* --- Ice: 3-phase AC before freeze; 2-phase AC after freeze --- */
            if (!user->flag_sed_frozen) {
                R_ice = N0[a] * ice_t;
                for (l = 0; l < dim; l++)
                    R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
                R_ice += C3 * ((Etased + Etaa)*fi - Etaa*fs - Etased*fa) * N0[a];
                R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;
            } else {
                PetscReal C2 = 3.0 * mob / (Etai + Etaa);
                R_ice = N0[a] * ice_t;
                for (l = 0; l < dim; l++)
                    R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);
                for (l = 0; l < dim; l++)
                    R_ice += C2 * Etaa * eps * (N1[a][l] * grad_sed[l]);
                R_ice += N0[a] * C2 / eps * (fi - fa);
                R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;
            }

            /* --- Sediment: 3-phase AC before freeze; AC + penalty after --- */
            R_sed = N0[a] * sed_t;
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob * eps * (N1[a][l] * grad_sed[l]);
            R_sed += C3 * (-Etaa*fi - Etai*fa + (Etai + Etaa)*fs) * N0[a];
            if (user->flag_sed_frozen)
                R_sed += k_sed * (sed - sed0) * N0[a];

            /* --- Thermal energy balance --- */
            R_tem  = rho * cp * N0[a] * tem_t;
            for (l = 0; l < dim; l++)
                R_tem += xi_T * thcond * (N1[a][l] * grad_tem[l]);
            R_tem += xi_T * rho * lat_sub * N0[a] * air_t;

            /* --- Vapor: penalised diffusivity + interface equilibrium penalty --- */
            R_vap  = N0[a] * (air_eff * rhov_t + air_t * rhov);
            for (l = 0; l < dim; l++)
                R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
            R_vap += k_pen * g_phiiphis * (rhov - rhov_eq) * N0[a];
            R_vap -= xi_v * N0[a] * rho_ice * air_t;
        }

        R[a][0] = R_ice;
        R[a][1] = R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;
    }
    return 0;
}


/* =========================================================================
 * Avenue 3: Cahn-Hilliard phase-field evolution (no penalty parameters)
 *
 * PDE (conservative, mass-conserving):
 *   ∂φ/∂t = ∇·(M ∇μ),   μ = (C3/ε)·f_drv(φ) − ε·∇²φ
 *
 * Weak form (biharmonic splitting, integrated by parts twice):
 *   R = N·φ_t + C3·∇N·∇f_drv + 3M·ε·ΔN·Δφ = 0
 *
 * where ΔN = Σₗ ∂²N/∂xₗ² (requires p ≥ 2, C ≥ 1) and
 *       Δφ = Σₗ ∂²φ/∂xₗ² (from IGAPointFormHess).
 *
 * Vapor: standard diffusion, no penalties.
 * Sediment: no freeze mechanism (always evolving).
 * ========================================================================= */
PetscErrorCode Residual_A3(IGAPoint pnt,
                            PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U,
                            PetscScalar *Re, void *ctx)
{
    PetscErrorCode ierr;
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    PetscReal eps     = user->eps;
    PetscReal Etai    = user->Etai;
    PetscReal Etaa    = user->Etaa;
    PetscReal Etased  = user->Etam;
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

    PetscScalar hess_flat[4 * 3 * 3];
    ierr = IGAPointFormHess(pnt, U, hess_flat); CHKERRQ(ierr);

    PetscScalar ice = sol[0], ice_t = sol_t[0];
    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];

    PetscScalar sed = sol[3], sed_t = sol_t[3];
    PetscScalar grad_sed[dim];
    for (l = 0; l < dim; l++) grad_sed[l] = grad_sol[3][l];

    PetscScalar air   = 1.0 - sed - ice;
    PetscScalar air_t = -ice_t - sed_t;

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

    PetscReal thcond, cp, rho, difvap, rhoI_vs;
    PetscReal fi, fa, fs, mob;
    ThermalCond(user, ice, sed, &thcond,  NULL);
    HeatCap    (user, ice, sed, &cp,      NULL);
    Density    (user, ice, sed, &rho,     NULL);
    VaporDiffus(user, tem,      &difvap,  NULL);
    RhoVS_I    (user, tem,      &rhoI_vs, NULL);
    Fice(user, ice, sed, &fi, NULL);
    Fair(user, ice, sed, &fa, NULL);
    Fsed(user, ice, sed, &fs, NULL);
    Mobility(user, ice, sed, &mob);

    const PetscReal *N0, (*N1)[dim], (*N2)[dim][dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    IGAPointGetShapeFuns(pnt, 2, (const PetscReal**)&N2);

    PetscScalar air_eff = (air > air_lim) ? air : air_lim;
    PetscReal loc = ice*ice * air*air;

    // ------------------------------------------------------------------
    // Laplacians of ice and sediment fields (trace of Hessian)
    // ------------------------------------------------------------------
    PetscReal lap_ice = 0.0, lap_sed = 0.0;
    for (l = 0; l < dim; l++) {
        lap_ice += hess_flat[0*dim*dim + l*dim + l];
        lap_sed += hess_flat[3*dim*dim + l*dim + l];
    }

    // ------------------------------------------------------------------
    // Partial derivatives of free energy functions
    // Fice, Fsed, Fair return value in first pointer and d/d(ice) in second
    // ChemPot_dsed returns d/d(sed) for all three
    // ------------------------------------------------------------------
    PetscReal dfi_dice, dfs_dice, dfa_dice;
    Fice(user, ice, sed, NULL, &dfi_dice);
    Fsed(user, ice, sed, NULL, &dfs_dice);
    Fair(user, ice, sed, NULL, &dfa_dice);

    PetscReal dfi_dsed, dfs_dsed, dfa_dsed;
    ChemPot_dsed(user, ice, sed, &dfi_dsed, &dfs_dsed, &dfa_dsed);

    // ------------------------------------------------------------------
    // Chemical potential differences (bulk driving forces)
    //
    // g_is = (3/eps) * (df/dphi_i - df/dphi_s)
    // g_ia = (3/eps) * (df/dphi_i - df/dphi_a)
    //
    // These replace the old EtaT-weighted combinations from Allen-Cahn.
    // The 3/eps factor is absorbed into the gradient terms below.
    //
    // Gradients of g_is and g_ia via chain rule:
    //   nabla(g_is) = dg_is/d(ice)*nabla(ice) + dg_is/d(sed)*nabla(sed)
    // ------------------------------------------------------------------
    PetscReal dg_is_dice = (3.0/eps) * (dfi_dice - dfs_dice);
    PetscReal dg_is_dsed = (3.0/eps) * (dfi_dsed - dfs_dsed);

    PetscReal dg_ia_dice = (3.0/eps) * (dfi_dice - dfa_dice);
    PetscReal dg_ia_dsed = (3.0/eps) * (dfi_dsed - dfa_dsed);

    // ------------------------------------------------------------------
    // Fourth-order surface energy coefficients (constant mobility M0)
    // From the derived weak form with Mi = Ms = Ma = M0, MT = 3*M0:
    //
    // Ice equation:
    //   phi_i coefficient: eps * M0 * (2*Etai + Etaa)
    //   phi_s coefficient: eps * M0 * (Etased - Etaa)    [subtracted]
    //
    // Sed equation:
    //   phi_i coefficient: eps * M0 * (Etai - Etaa)
    //   phi_s coefficient: eps * M0 * (2*Etased + Etaa)  [subtracted]
    // ------------------------------------------------------------------
    PetscReal coeff_ii = eps * mob * (2.0*Etai   + Etaa);
    PetscReal coeff_is = eps * mob * (Etased - Etaa);
    PetscReal coeff_si = eps * mob * (Etai   - Etaa);
    PetscReal coeff_ss = eps * mob * (2.0*Etased + Etaa);

    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal R_ice = 0.0, R_sed = 0.0, R_tem = 0.0, R_vap = 0.0;

        if (user->flag_relax) {
            /* RELAXATION MODE — Residual_A3 (Cahn-Hilliard) is not used
             * during relaxation; zero all DOFs and skip CH assembly. */
        } else {
            // Laplacian of test function
            PetscReal lap_Na = 0.0;
            for (l = 0; l < dim; l++) lap_Na += N2[a][l][l];

            // ----------------------------------------------------------
            // Ice equation — Cahn-Hilliard weak form (W1)
            //
            // R_ice = N*phi_i_t
            //       + (M0/3)*nabla(N)*nabla(g_is + g_ia)
            //       + eps*M0*(2*Etai+Etaa)*lap(N)*lap(phi_i)
            //       - eps*M0*(Etased-Etaa)*lap(N)*lap(phi_s)
            //       - N*alph_sub*loc*(rhov-rhoI_vs)/rho_ice
            // ----------------------------------------------------------
            R_ice = N0[a] * ice_t;

            // Bulk gradient term: (M0/3) * nabla(N) . nabla(g_is + g_ia)
            for (l = 0; l < dim; l++)
                R_ice += (mob/3.0) * N1[a][l]
                       * ((dg_is_dice + dg_ia_dice) * grad_ice[l]
                        + (dg_is_dsed + dg_ia_dsed) * grad_sed[l]);

            // Fourth-order terms: eps*M0*(2Etai+Etaa)*lap(N)*lap(phi_i)
            //                   - eps*M0*(Etased-Etaa)*lap(N)*lap(phi_s)
            R_ice += coeff_ii * lap_Na * lap_ice;
            R_ice -= coeff_is * lap_Na * lap_sed;

            // Sublimation source
            R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;

            // ----------------------------------------------------------
            // Sediment equation — Cahn-Hilliard weak form (W2)
            //
            // R_sed = N*phi_s_t
            //       - (2*M0/3)*nabla(N)*nabla(g_is)
            //       + (M0/3)*nabla(N)*nabla(g_ia)
            //       + eps*M0*(Etai-Etaa)*lap(N)*lap(phi_i)
            //       - eps*M0*(2*Etased+Etaa)*lap(N)*lap(phi_s)
            // ----------------------------------------------------------
            R_sed = N0[a] * sed_t;

            // Bulk gradient terms: -(2M0/3)*nabla(N).nabla(g_is)
            //                    + (M0/3)*nabla(N).nabla(g_ia)
            for (l = 0; l < dim; l++)
                R_sed += N1[a][l]
                       * ((-2.0*mob/3.0) * (dg_is_dice * grad_ice[l]
                                           + dg_is_dsed * grad_sed[l])
                        + ( mob/3.0)     * (dg_ia_dice * grad_ice[l]
                                           + dg_ia_dsed * grad_sed[l]));

            // Fourth-order terms: eps*M0*(Etai-Etaa)*lap(N)*lap(phi_i)
            //                   - eps*M0*(2Etased+Etaa)*lap(N)*lap(phi_s)
            R_sed += coeff_si * lap_Na * lap_ice;
            R_sed -= coeff_ss * lap_Na * lap_sed;

            // ----------------------------------------------------------
            // Thermal energy balance
            // ----------------------------------------------------------
            R_tem  = rho * cp * N0[a] * tem_t;
            for (l = 0; l < dim; l++)
                R_tem += xi_T * thcond * (N1[a][l] * grad_tem[l]);
            R_tem += xi_T * rho * lat_sub * N0[a] * air_t;

            // ----------------------------------------------------------
            // Vapor transport
            // Standard diffusion — xi_v scales diffusion only,
            // not the phase change source term
            // ----------------------------------------------------------
            R_vap  = N0[a] * air_eff * rhov_t;
            for (l = 0; l < dim; l++)
                R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
            R_vap -= N0[a] * rho_ice * air_t;
        }

        R[a][0] = R_ice;
        R[a][1] = R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;
    }
    return 0;
}


/* =========================================================================
 * Dispatcher: routes to the avenue selected by user->flag_avenue
 *   1 = Avenue 1 (AC, penalty vapor, freeze → zero RHS)
 *   2 = Avenue 2 (AC, penalty vapor, freeze → penalty) [default]
 *   3 = Avenue 3 (Cahn-Hilliard, no penalties)
 * ========================================================================= */
PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;
    switch (user->flag_avenue) {
        case 1:  return Residual_A1(pnt, shift, V, t, U, Re, ctx);
        case 2:  return Residual_A2(pnt, shift, V, t, U, Re, ctx);
        case 3:  return Residual_A3(pnt, shift, V, t, U, Re, ctx);
        default:
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Unknown -flag_avenue %d (valid: 1, 2, 3)", user->flag_avenue);
    }
}


/* =========================================================================
 * Jacobian — analytical Jacobian not implemented; FD approximation is used.
 * ========================================================================= */
PetscErrorCode Jacobian(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Je, void *ctx)
{
    return 0;
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
