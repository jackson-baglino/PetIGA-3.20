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
 *                                  + xi_T rho(phi) Lsub dphi_i/dt
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
 * The Jacobian for A2 is (for now) formed by finite differences
 * (IGAFormIJacobianFD, selected in permafrost2.c for dof=4). The analytic
 * Jacobian_A2 is the next step, to be verified with -snes_test_jacobian.
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

        /* Temperature (mixture rho for latent heat, xi_T scaling). */
        R[a][1] = rw * ( rho * cp * N0[a] * tem_t
                + user->xi_T * thcond * gN_gtem
                - pc * user->xi_T * rho * lat_sub * phi_t * N0[a] );

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
