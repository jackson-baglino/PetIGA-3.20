/* =========================================================================
 * assembly_sed.c — 3-phase (ice / sediment / air) weak form, "Avenue 2" (A2).
 *
 * DOF layout: 0 = phi_i (ice), 1 = T, 2 = rho_v (vapor), 3 = phi_s (sediment).
 * phi_a = 1 - phi_i - phi_s is algebraic. Sediment is a STATIC field
 * (d phi_s/dt = 0), but carried as a DOF so it is stored/output.
 *
 * Governing equations (user-supplied, beta-eliminated Lagrange multiplier):
 *
 *  Ice:   dphi_i/dt = -3 M0/(Si+Sa) [ (1/eps)(dFtri/dphi_i - dFtri/dphi_a)
 *                                     - eps(Si+Sa) lap(phi_i) - eps Sa lap(phi_s) ]
 *                     + alpha_sub phi_i^2 phi_a^2 (rho_v - rho_vs^I(T))/rho_i
 *  Sed:   dphi_s/dt = 0
 *  Temp:  rho(phi) cp(phi) dT/dt = xi_T div(K(phi) grad T)
 *                                  + xi_T rho_ice Lsub dphi_i/dt
 *         (latent-heat density = rho_ice, not the mixture rho of the written
 *          equation — chosen after an A1 A/B test; see the R_tem comment below)
 *  Vapor: d(phi_a rho_v)/dt = xi_v div(phi_a Dv(T) grad rho_v) - rho_i dphi_i/dt
 *
 * Weak residual (R = 0), integrated by parts on the Laplacian/divergence terms
 * (no-flux boundaries), with g = dFtri/dphi_i - dFtri/dphi_a (TripleWell):
 *
 *  R_ice = N phi_i_t + (3 M0/((Si+Sa) eps)) g N
 *        + 3 M0 eps (grad_phi_i . grad_N)
 *        + 3 M0 (Sa/(Si+Sa)) eps (grad_phi_s . grad_N)
 *        - (alpha_sub/rho_i) loc (rho_v - rho_vs) N
 *  R_sed = N phi_s_t
 *  R_tem = rho cp N T_t + xi_T K (grad_T . grad_N) - xi_T rho Lsub phi_i_t N
 *  R_vap = phi_aef N rho_v_t + xi_v Dv phi_aef (grad_rhov . grad_N)
 *        + (xi_v rho_i - rho_v) phi_i_t N
 *
 * This reduces EXACTLY to Residual_A1 when phi_s = 0 (TripleWell -> the 2-phase
 * double well; the cross grad_phi_s term vanishes; mixture props -> 2-phase).
 *
 * Jacobian_A2 (below) is the analytic Jacobian; verify with -snes_test_jacobian
 * at t=0 (the IC is unclamped there, so it matches the FD Jacobian to round-off).
 * ========================================================================= */

#include "assembly.h"
#include "material_properties.h"

PetscErrorCode Residual_A2(IGAPoint pnt,
                           PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U,
                           PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    PetscReal eps     = user->eps;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;
    PetscReal Sigma_i = user->Sigma_i;
    PetscReal Sigma_a = user->Sigma_a;
    PetscReal Ssum    = Sigma_i + Sigma_a;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar phi   = sol[0],   phi_t   = sol_t[0];
    PetscScalar tem   = sol[1],   tem_t   = sol_t[1];
    PetscScalar rhov  = sol[2],   rhov_t  = sol_t[2];
    PetscScalar phi_s = sol[3],   phi_s_t = sol_t[3];
    PetscScalar grad_phi  [dim], grad_tem [dim];
    PetscScalar grad_rhov [dim], grad_phi_s[dim];
    for (l = 0; l < dim; l++) {
        grad_phi  [l] = grad_sol[0][l];
        grad_tem  [l] = grad_sol[1][l];
        grad_rhov [l] = grad_sol[2][l];
        grad_phi_s[l] = grad_sol[3][l];
    }
    PetscScalar phi_a = 1.0 - phi - phi_s;

    /* SNES domain-error catch: back off the line search if an iterate leaves
     * the configured phase bounds. Sediment is static, so it is not checked. */
    {
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (PetscRealPart(phi)   < lo || PetscRealPart(phi)   > hi ||
            PetscRealPart(phi_a) < lo || PetscRealPart(phi_a) > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    /* Clamped fractions for material-property evaluation. */
    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_sc = PetscRealPart(phi_s);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;  if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_sc < 0.0) phi_sc = 0.0;  if (phi_sc > 1.0) phi_sc = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;  if (phi_ac > 1.0) phi_ac = 1.0;

    /* Material properties (3-phase mixtures over ice/sed; air = 1-ice-sed). */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond3(user, phi_c, phi_sc, &thcond,  NULL, NULL);
    HeatCap3    (user, phi_c, phi_sc, &cp,      NULL, NULL);
    Density3    (user, phi_c, phi_sc, &rho,     NULL, NULL);
    VaporDiffus (user, tem,           &dif_vap, NULL);
    Mobility    (user, phi_c,         &mob_sub);

    /* Flat-interface saturation vapor density. */
    PetscReal rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, NULL);

    /* Triple-well driving force g = dFtri/dphi_i - dFtri/dphi_a. */
    PetscScalar gtw;
    TripleWell(user, phi_c, phi_sc, &gtw, NULL, NULL);

    /* Sublimation localizer (ice-air interface) and vapor-storage floor. */
    PetscReal loc     = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal phi_aef = (phi_ac > air_lim) ? phi_ac : air_lim;

    /* Ice-equation prefactors. */
    PetscReal c_gi = 3.0 * mob_sub * eps;                 /* grad_phi_i coeff */
    PetscReal c_gs = 3.0 * mob_sub * eps * Sigma_a / Ssum;/* grad_phi_s coeff */
    PetscReal c_g  = 3.0 * mob_sub / (Ssum * eps);        /* triple-well coeff */

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a, nen = pnt->nen;

    /* Axisymmetric r-weight (dV = 2*pi*r dr dz; the 2*pi cancels in R=0). */
    PetscReal rw = 1.0;
    if (user->axisym) {
        PetscReal xphys[3] = {0.0, 0.0, 0.0};
        IGAPointFormPoint(pnt, xphys);
        rw = xphys[1];
    }

    for (a = 0; a < nen; a++) {
        PetscReal gN_gphi  = 0.0, gN_gtem  = 0.0;
        PetscReal gN_grhov = 0.0, gN_gphis = 0.0;
        for (l = 0; l < dim; l++) {
            gN_gphi  += N1[a][l] * PetscRealPart(grad_phi  [l]);
            gN_gtem  += N1[a][l] * PetscRealPart(grad_tem  [l]);
            gN_grhov += N1[a][l] * PetscRealPart(grad_rhov [l]);
            gN_gphis += N1[a][l] * PetscRealPart(grad_phi_s[l]);
        }

        /* -decouple_phase_change 1: zero every phase-change coupling. */
        const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

        /* Ice */
        R[a][0] = rw * ( N0[a] * phi_t
                + c_g * PetscRealPart(gtw) * N0[a]
                + c_gi * gN_gphi
                + c_gs * gN_gphis
                - pc * (user->alph_sub / rho_ice) * loc
                  * (PetscRealPart(rhov) - rho_vs) * N0[a] );

        /* Temperature (xi_T scaling). Latent-heat density is rho_ice, NOT the
         * mixture rho of the written equation: an A1 A/B test showed mixture rho
         * gave ~3x worse vapor-mass conservation, and rho_ice is the physically
         * correct energy-per-unit-ice-mass form. Same choice as A1. */
        R[a][1] = rw * ( rho * cp * N0[a] * tem_t
                + user->xi_T * thcond * gN_gtem
                - pc * user->xi_T * rho_ice * lat_sub * phi_t * N0[a] );

        /* Vapor (phi_a = 1-phi_i-phi_s storage/diffusion, xi_v scaling). */
        R[a][2] = rw * ( phi_aef * N0[a] * rhov_t
                + user->xi_v * dif_vap * phi_aef * gN_grhov
                + pc * (user->xi_v * rho_ice - PetscRealPart(rhov))
                  * phi_t * N0[a] );

        /* Sediment: d phi_s/dt = 0. */
        R[a][3] = rw * ( N0[a] * phi_s_t );
    }
    return 0;
}


/* =========================================================================
 * Jacobian_A2 — analytic Jacobian of Residual_A2.
 *   J[a][k][b][l] = dR[a][k]/du[b][l] + shift * dR[a][k]/du_t[b][l]
 * DOF: 0=phi_i, 1=T, 2=rho_v, 3=phi_s.
 *
 * Derived block-by-block from Residual_A2 (see that function). Material and
 * triple-well derivatives come from ThermalCond3/HeatCap3/Density3/TripleWell.
 * Exact in the interior (all phi in [0,1]); in the thin clamp bands the
 * property derivatives are treated as unclamped (a Newton-convergence
 * approximation, not a solution error). Verify with -snes_test_jacobian at t=0,
 * where the IC is unclamped so this matches the FD Jacobian to round-off.
 * ========================================================================= */
PetscErrorCode Jacobian_A2(IGAPoint pnt,
                           PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U,
                           PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;

    PetscInt l, dim = user->dim;
    PetscReal eps     = user->eps;
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;
    PetscReal Sigma_i = user->Sigma_i, Sigma_a = user->Sigma_a;
    PetscReal Ssum    = Sigma_i + Sigma_a;

    if (pnt->atboundary) return 0;

    PetscScalar sol_t[4], sol[4], grad_sol[4][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar phi   = sol[0],  phi_t  = sol_t[0];
    PetscScalar tem   = sol[1],  tem_t  = sol_t[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    PetscScalar phi_s = sol[3];
    PetscScalar grad_tem[dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) { grad_tem[l] = grad_sol[1][l]; grad_rhov[l] = grad_sol[2][l]; }
    PetscScalar phi_a = 1.0 - phi - phi_s;

    /* Clamped fractions (match the residual). */
    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_sc = PetscRealPart(phi_s);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;  if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_sc < 0.0) phi_sc = 0.0;  if (phi_sc > 1.0) phi_sc = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;  if (phi_ac > 1.0) phi_ac = 1.0;

    /* Material properties and their derivatives w.r.t. (phi_i, phi_s). */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    PetscReal dK_di, dK_ds, dcp_di, dcp_ds, drho_di, drho_ds, d_dif_vap;
    ThermalCond3(user, phi_c, phi_sc, &thcond,  &dK_di,   &dK_ds);
    HeatCap3    (user, phi_c, phi_sc, &cp,       &dcp_di,  &dcp_ds);
    Density3    (user, phi_c, phi_sc, &rho,      &drho_di, &drho_ds);
    VaporDiffus (user, tem,           &dif_vap,  &d_dif_vap);
    Mobility    (user, phi_c,         &mob_sub);

    PetscReal rho_vs, d_rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, &d_rho_vs);

    PetscScalar gtw, dg_di, dg_ds;
    TripleWell(user, phi_c, phi_sc, &gtw, &dg_di, &dg_ds);

    /* Localizer loc = phi_i^2 phi_a^2 and its derivatives (phi_a=1-phi_i-phi_s). */
    PetscReal loc     = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal dloc_di = 2.0 * phi_c * phi_ac * (phi_ac - phi_c);
    PetscReal dloc_ds = -2.0 * phi_c * phi_c * phi_ac;

    /* Vapor storage/diffusion floor and its derivative (=-1 for both DOFs
     * above the floor, 0 when floored). */
    PetscBool above   = (phi_ac > air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal phi_aef = above ? phi_ac : air_lim;
    PetscReal dphiaef = above ? -1.0 : 0.0;

    PetscReal c_gi = 3.0 * mob_sub * eps;
    PetscReal c_gs = 3.0 * mob_sub * eps * Sigma_a / Ssum;
    PetscReal c_g  = 3.0 * mob_sub / (Ssum * eps);
    PetscReal ar   = user->alph_sub / rho_ice;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;

    PetscReal rw = 1.0;
    if (user->axisym) { PetscReal xp[3] = {0.0,0.0,0.0}; IGAPointFormPoint(pnt, xp); rw = xp[1]; }

    const PetscReal pc      = user->decouple_phase_change ? 0.0 : 1.0;
    const PetscReal Rrhov   = PetscRealPart(rhov);
    const PetscReal Rrhov_t = PetscRealPart(rhov_t);
    const PetscReal Rtem_t  = PetscRealPart(tem_t);
    const PetscReal Rphi_t  = PetscRealPart(phi_t);
    const PetscReal dsrc    = Rrhov - rho_vs;   /* rho_v - rho_vs */

    for (a = 0; a < nen; a++) {
        PetscReal gNa_gtem = 0.0, gNa_grhov = 0.0;
        for (l = 0; l < dim; l++) {
            gNa_gtem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            gNa_grhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }
        for (b = 0; b < nen; b++) {
            PetscReal NaNb = N0[a] * N0[b];
            PetscReal gNagNb = 0.0;
            for (l = 0; l < dim; l++) gNagNb += N1[a][l] * N1[b][l];

            /* ---- Row 0: ice ---- */
            J[a][0][b][0] += rw * ( shift * NaNb
                           + c_g * PetscRealPart(dg_di) * NaNb
                           + c_gi * gNagNb
                           - pc * ar * dloc_di * dsrc * NaNb );
            J[a][0][b][1] += rw * ( pc * ar * loc * d_rho_vs * NaNb );
            J[a][0][b][2] += rw * ( -pc * ar * loc * NaNb );
            J[a][0][b][3] += rw * ( c_g * PetscRealPart(dg_ds) * NaNb
                           + c_gs * gNagNb
                           - pc * ar * dloc_ds * dsrc * NaNb );

            /* ---- Row 1: temperature ---- */
            J[a][1][b][0] += rw * ( (drho_di * cp + rho * dcp_di) * Rtem_t * NaNb
                           + user->xi_T * dK_di * gNa_gtem * N0[b]
                           - pc * user->xi_T * rho_ice * lat_sub * shift * NaNb );
            J[a][1][b][1] += rw * ( shift * rho * cp * NaNb
                           + user->xi_T * thcond * gNagNb );
            J[a][1][b][3] += rw * ( (drho_ds * cp + rho * dcp_ds) * Rtem_t * NaNb
                           + user->xi_T * dK_ds * gNa_gtem * N0[b] );

            /* ---- Row 2: vapor ---- */
            J[a][2][b][0] += rw * ( dphiaef * Rrhov_t * NaNb
                           + user->xi_v * dif_vap * dphiaef * gNa_grhov * N0[b]
                           + pc * (user->xi_v * rho_ice - Rrhov) * shift * NaNb );
            J[a][2][b][1] += rw * ( user->xi_v * d_dif_vap * phi_aef * gNa_grhov * N0[b] );
            J[a][2][b][2] += rw * ( shift * phi_aef * NaNb
                           + user->xi_v * dif_vap * phi_aef * gNagNb
                           - pc * Rphi_t * NaNb );
            J[a][2][b][3] += rw * ( dphiaef * Rrhov_t * NaNb
                           + user->xi_v * dif_vap * dphiaef * gNa_grhov * N0[b] );

            /* ---- Row 3: sediment (d phi_s/dt = 0) ---- */
            J[a][3][b][3] += rw * ( shift * NaNb );
        }
    }
    return 0;
}
