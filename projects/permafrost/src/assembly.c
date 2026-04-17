#include "assembly.h"
#include "material_properties.h"

PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*) ctx;

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

    PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;

    PetscReal mob, mob_sed, alph_sub;
    if (user->flag_Tdep == 1) {
        mob      = user->mob[indGP];
        mob_sed  = user->mob_sed_arr[indGP];
        alph_sub = user->alph[indGP];
    } else {
        mob      = user->mob_sub;
        mob_sed  = user->mob_sed;
        alph_sub = user->alph_sub;
    }

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar ice   = sol[0], ice_t = sol_t[0];
    PetscScalar grad_ice[dim];
    for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];

    PetscScalar sed   = sol[3], sed_t = sol_t[3];
    PetscScalar grad_sed[dim];
    for (l = 0; l < dim; l++) grad_sed[l] = grad_sol[3][l];

    // Air: algebraic from constraint, no independent evolution equation
    PetscScalar air   = 1.0 - sed - ice;
    PetscScalar air_t = -ice_t - sed_t;

    if (ice < 0.0) ice = 0.0;  // Prevent negative ice from causing NaN diffusion coefficients
    if (ice > 1.0) ice = 1.0;  // Prevent unphysical >100% ice from causing NaN diffusion coefficients
    if (air < 0.0) air = 0.0;  // Prevent negative air from causing NaN diffusion coefficients
    if (air > 1.0) air = 1.0;  // Prevent unphysical >100% air from causing NaN diffusion coefficients

    PetscScalar tem   = sol[1], tem_t = sol_t[1];
    PetscScalar grad_tem[dim];
    for (l = 0; l < dim; l++) grad_tem[l] = grad_sol[1][l];

    PetscScalar rhov  = sol[2], rhov_t = sol_t[2];
    PetscScalar grad_rhov[dim];
    for (l = 0; l < dim; l++) grad_rhov[l] = grad_sol[2][l];

    PetscReal thcond, cp, rho, difvap, rhoI_vs, fi, fa, fs;
    ThermalCond(user, ice, sed, &thcond,  NULL);
    HeatCap    (user, ice, sed, &cp,      NULL);
    Density    (user, ice, sed, &rho,     NULL);
    VaporDiffus(user, tem,      &difvap,  NULL);
    RhoVS_I    (user, tem,      &rhoI_vs, NULL);
    Fice(user, ice, sed, &fi, NULL);
    Fair(user, ice, sed, &fa, NULL);
    Fsed(user, ice, sed, &fs, NULL);

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    // air_eff: regularized air for vapor mass matrix and diffusion
    // Prevents negative diffusion coefficient when air goes slightly negative
    PetscScalar air_eff = (air > air_lim) ? air : air_lim;

    // Sublimation localization: nonzero only at ice-air interface,
    // suppressed at ice-sediment interface by (1-sed)^2
    // PetscReal loc = ice*ice * (1 - ice - sed)*(1 - ice - sed);

    // Ice: Allen-Cahn with Lagrange multiplier to enfroce \phi_i_t + \phi_a_t = 0
    // PetscReal C = 3.0 * mob / (Etai + Etaa);
    // PetscReal C = 3.0 * mob / Etai;

    // Phase-dependent effective mobilities: M(φ) = m₀ · φ · (1 − φ)
    // Localizes diffusion to the diffuse interface; zero in bulk single-phase regions
    PetscReal mob_eff     = mob     * ice * (1.0 - ice);
    PetscReal mob_sed_eff = mob_sed * sed * (1.0 - sed);

    // Ice: Allen-Cahn with Lagrange multiplier to enforce \phi_i_t + \phi_a_t + \phi_s_t = 0
    // Separate C3 per phase so each uses its own effective mobility
    PetscReal C3_ice = 3.0 * mob_eff     / (eps * EtaT);
    PetscReal C3_sed = 3.0 * mob_sed_eff / (eps * EtaT);

    PetscReal loc = ice*ice * air*air;
    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;
    for (a = 0; a < nen; a++) {
        PetscReal R_ice = 0.0, R_sed = 0.0, R_tem = 0.0, R_vap = 0.0;

        if (user->flag_tIC == 1) {
            R_ice = 0.0; R_sed = 0.0; R_tem = 0.0; R_vap = 0.0;
        } else {
            /* Ice Evolution Equation (2-phase; with Lagrange multiplier to enfroce \phi_i_t + \phi_a_t = 0) */
            // // Time derivative
            // R_ice = N0[a] * ice_t;

            // // phi_i Laplacian mobility term
            // for (l = 0; l < dim; l++)
            //     R_ice += 3.0 * mob * eps * (N1[a][l] * grad_ice[l]);

            // // phi_s Laplacian mobility term (known field, contributes through constraint)
            // for (l = 0; l < dim; l++)
            //     R_ice += C * Etaa * eps * (N1[a][l] * grad_sed[l]);

            // // Bulk driving force
            // R_ice += N0[a] * C / eps * (fi - fa);

            // // Sublimation source — disabled pending vapor equation validation
            // R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;

            /* Ice Evolution Equation (3-phase; with Lagrange multiplier to enforce \phi_i_t + \phi_a_t + \phi_s_t = 0) */
            // Time derivative
            R_ice = N0[a] * ice_t;

            // phi_i Laplacian term (phase-dependent mobility)
            for (l = 0; l < dim; l++)
                R_ice += 3.0 * mob_eff * eps * (N1[a][l] * grad_ice[l]);

            // Bulk driving force
            R_ice += C3_ice*((Etased + Etaa)*fi - Etaa*fs - Etased*fa) * N0[a];

            // Sublimation source term
            R_ice -= N0[a] * alph_sub * loc * (rhov - rhoI_vs) / rho_ice;

            /* Sediment Evolution Equation */
            // Time derivative
            R_sed = N0[a] * sed_t;

            // Laplacian term (phase-dependent mobility)
            for (l = 0; l < dim; l++)
                R_sed += 3.0 * mob_sed_eff * eps * (N1[a][l] * grad_sed[l]);

            // Bulk driving force
            R_sed += C3_sed*(-Etaa*fi - Etai*fa + (Etai + Etaa)*fs) * N0[a];

            /* Thermal energy balance */
            // Time derivative
            // R_tem  = N0[a] * tem_t;             // Does not solve temperature equation
            R_tem  = rho * cp * N0[a] * tem_t;

            // Diffusion term with effective thermal conductivity
            for (l = 0; l < dim; l++)
                R_tem += xi_T * thcond * (N1[a][l] * grad_tem[l]);

            // Source term from phase change (sublimation latent heat)
            R_tem += xi_T * rho * lat_sub * N0[a] * air_t;

            // Vapor: time derivative and diffusion use air_eff to prevent
            // negative diffusion coefficient; source from phase change via rho_ice
            // Time derivative
            // R_vap  = N0[a] * rhov_t;                // Does not solve vapor equation
            R_vap  = N0[a] * air_eff * rhov_t;

            // Diffusion term with effective diffusivity
            for (l = 0; l < dim; l++)
                R_vap += xi_v * difvap * air_eff * (N1[a][l] * grad_rhov[l]);
            
            // Source term from phase change (sublimation/condensation)
            R_vap -= xi_v * N0[a] * rho_ice * air_t;
        }

        R[a][0] = R_ice;
        R[a][1] = R_tem;
        R[a][2] = R_vap;
        R[a][3] = R_sed;   /* ∂sed/∂t = 0  (mob_sed available; full RHS to be added) */
    }

    return 0;
}


PetscErrorCode Jacobian(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Je, void *ctx)
{
//     AppCtx *user = (AppCtx*) ctx;

//     PetscInt l, dim = user->dim;
//     PetscReal eps     = user->eps;
//     PetscReal Etai    = user->Etai;
//     PetscReal rho_ice = user->rho_ice;
//     PetscReal lat_sub = user->lat_sub;
//     PetscReal air_lim = user->air_lim;
//     PetscReal xi_v    = user->xi_v;
//     PetscReal xi_T    = user->xi_T;

//     PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;

//     PetscReal mob, mob_sed, alph_sub;
//     if (user->flag_Tdep == 1) {
//         mob      = user->mob[indGP];
//         mob_sed  = user->mob_sed_arr[indGP];
//         alph_sub = user->alph[indGP];
//     } else {
//         mob      = user->mob_sub;
//         mob_sed  = user->mob_sed;
//         alph_sub = user->alph_sub;
//     }

//     if (pnt->atboundary) return 0;

//     PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
//     IGAPointFormValue(pnt, V, &sol_t[0]);
//     IGAPointFormValue(pnt, U, &sol[0]);
//     IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

//     PetscScalar ice   = sol[0], ice_t = sol_t[0];
//     PetscScalar grad_ice[dim];
//     for (l = 0; l < dim; l++) grad_ice[l] = grad_sol[0][l];

//     PetscScalar sed   = sol[3], sed_t = sol_t[3];

//     PetscScalar air   = 1.0 - sed - ice;
//     PetscScalar air_t = -ice_t - sed_t;

//     PetscScalar tem   = sol[1], tem_t = sol_t[1];
//     PetscScalar grad_tem[dim];
//     for (l = 0; l < dim; l++) grad_tem[l] = grad_sol[1][l];

//     PetscScalar rhov  = sol[2], rhov_t = sol_t[2];
//     PetscScalar grad_rhov[dim];
//     for (l = 0; l < dim; l++) grad_rhov[l] = grad_sol[2][l];

//     PetscReal thcond, dthcond_ice, cp, dcp_ice, rho, drho_ice;
//     PetscReal difvap, d_difvap, rhoI_vs, drhoI_vs;
//     PetscReal fice_ice, fair_ice;
//     ThermalCond(user, ice, sed, &thcond,  &dthcond_ice);
//     HeatCap    (user, ice, sed, &cp,      &dcp_ice);
//     Density    (user, ice, sed, &rho,     &drho_ice);
//     VaporDiffus(user, tem,      &difvap,  &d_difvap);
//     RhoVS_I    (user, tem,      &rhoI_vs, &drhoI_vs);
//     Fice(user, ice, sed, NULL, &fice_ice);

//     const PetscReal *N0, (*N1)[dim];
//     IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
//     IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

//     PetscReal loc      = ice*ice * air*air * (1.0-sed)*(1.0-sed);
//     PetscReal dloc_ice = (2.0*ice*air*air - 2.0*ice*ice*air) * (1.0-sed)*(1.0-sed);

//     (void)mob_sed;  /* reserved for sediment Jacobian entries */
//     PetscInt a, b, nen = pnt->nen;
//     PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;

//     for (a = 0; a < nen; a++) {
//         for (b = 0; b < nen; b++) {
//             if (user->flag_tIC == 1) continue;

//             // Ice
//             J[a][0][b][0] += shift * N0[a] * N0[b];
//             for (l = 0; l < dim; l++)
//                 J[a][0][b][0] += 3.0 * mob * eps * (N1[a][l] * N1[b][l]);
//             J[a][0][b][0] += N0[a] * mob * 3.0/eps/Etai * (fice_ice - fair_ice) * N0[b];
//             J[a][0][b][0] -= N0[a] * alph_sub * dloc_ice * N0[b] * (rhov - rhoI_vs) / rho_ice;
//             J[a][0][b][1] += N0[a] * alph_sub * loc * drhoI_vs * N0[b] / rho_ice;
//             J[a][0][b][2] -= N0[a] * alph_sub * loc * N0[b] / rho_ice;

//             // Temperature
//             J[a][1][b][1] += shift * N0[a] * N0[b];
//             // J[a][1][b][1] += shift * rho * cp * N0[a] * N0[b];
//             // J[a][1][b][0] += drho_ice * N0[b] * cp * N0[a] * tem_t;
//             // J[a][1][b][0] += rho * dcp_ice * N0[b] * N0[a] * tem_t;
//             // for (l = 0; l < dim; l++)
//             //     J[a][1][b][0] += xi_T * dthcond_ice * N0[b] * (N1[a][l] * grad_tem[l]);
//             // for (l = 0; l < dim; l++)
//             //     J[a][1][b][1] += xi_T * thcond * (N1[a][l] * N1[b][l]);
//             // J[a][1][b][0] += xi_T * drho_ice * N0[b] * lat_sub * N0[a] * air_t;
//             // J[a][1][b][0] -= xi_T * rho * lat_sub * N0[a] * shift * N0[b];

//             // Vapor — time derivative
//             J[a][2][b][2] += shift * N0[a] * N0[b];

//             // Sediment — time derivative only (RHS to be added)
//             J[a][3][b][3] += shift * N0[a] * N0[b];
//             // if (air > air_lim) {
//             //     J[a][2][b][2] += N0[a] * air * shift * N0[b];
//             //     J[a][2][b][0] -= N0[a] * rhov_t * N0[b];
//             // } else {
//             //     J[a][2][b][2] += N0[a] * air_lim * shift * N0[b];
//             // }

//             // // Vapor — diffusion (air_eff in both branches, consistent with residual)
//             // if (air > air_lim) {
//             //     for (l = 0; l < dim; l++)
//             //         J[a][2][b][0] -= xi_v * difvap * (N1[a][l] * grad_rhov[l]) * N0[b];
//             //     for (l = 0; l < dim; l++)
//             //         J[a][2][b][2] += xi_v * difvap * air * (N1[a][l] * N1[b][l]);
//             //     for (l = 0; l < dim; l++)
//             //         J[a][2][b][1] += xi_v * d_difvap * N0[b] * air * (N1[a][l] * grad_rhov[l]);
//             // } else {
//             //     for (l = 0; l < dim; l++)
//             //         J[a][2][b][2] += xi_v * difvap * air_lim * (N1[a][l] * N1[b][l]);
//             //     for (l = 0; l < dim; l++)
//             //         J[a][2][b][1] += xi_v * d_difvap * N0[b] * air_lim * (N1[a][l] * grad_rhov[l]);
//             // }

//             // // Vapor — source term: -xi_v * rho_ice * air_t = +xi_v * rho_ice * ice_t
//             // J[a][2][b][0] += N0[a] * xi_v * rho_ice * shift * N0[b];
//         }
//     }

    return 0;
}


PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                            PetscScalar *S, void *ctx)
{
    PetscFunctionBegin;
    // AppCtx *user = (AppCtx*)ctx;
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

    PetscFunctionReturn(0);
}