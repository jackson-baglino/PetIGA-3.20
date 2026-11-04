#include "assembly.h"
#include "material_properties.h"

/* Weak form of the 2-phase (ice / T / vapor) system. DOFs: 0 = phi (ice),
 * 1 = T, 2 = rhov. phi_a = 1 - phi. Derivations: docs/. */

/* d/dphi of the half-normalised well (1/2)*phi^2*(1-phi)^2, which makes the
 * interface width parameter equal eps. Do not double it. */
static void DoubleWellDeriv(PetscReal phi, PetscReal *f1, PetscReal *df1)
{
    if (f1)  *f1  = phi * (1.0 - phi) * (1.0 - 2.0 * phi);
    if (df1) *df1 = 1.0 - 6.0 * phi + 6.0 * phi * phi;
}

/* Wall energy interpolant h = phi^2*(3-2*phi) and its derivatives. */
static void WallH(PetscReal phi, PetscReal *h, PetscReal *dh, PetscReal *d2h)
{
    if (h)   *h   = phi * phi * (3.0 - 2.0 * phi);
    if (dh)  *dh  = 6.0 * phi * (1.0 - phi);
    if (d2h) *d2h = 6.0 * (1.0 - 2.0 * phi);
}

typedef enum { FACE_NATURAL, FACE_WALL, FACE_OPEN } FaceKind;

/* What the phi equation does on this boundary point's face:
 *   FACE_WALL    regolith, dphi/dn = cos(theta)*phi*(1-phi)/eps
 *   FACE_OPEN    ice passes through; nothing is imposed on phi
 *   FACE_NATURAL dphi/dn = 0 */
static FaceKind FaceKindAt(IGAPoint pnt, const AppCtx *user)
{
    PetscInt axis = -1, side = -1;
    IGAPointAtBoundary(pnt, &axis, &side);
    if (axis < 0 || axis > 2 || side < 0 || side > 1) return FACE_NATURAL;
    if (user->wall_face[axis][side]) return (user->costhet != 0.0) ? FACE_WALL : FACE_NATURAL;
    if (user->open_face[axis][side]) return FACE_OPEN;
    return FACE_NATURAL;
}

/* Axisymmetric weight r (the 2*pi cancels in R = 0). */
static PetscReal AxisymWeight(IGAPoint pnt, const AppCtx *user)
{
    PetscReal xphys[3] = {0.0, 0.0, 0.0};
    if (!user->axisym) return 1.0;
    IGAPointFormPoint(pnt, xphys);
    return xphys[1];
}


PetscErrorCode Residual_A1(IGAPoint pnt,
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

    /* Boundary faces: these terms add to the volume form. */
    if (pnt->atboundary) {
        FaceKind kind = FaceKindAt(pnt, user);
        if (kind == FACE_NATURAL) return 0;

        PetscScalar sol_b[3];
        IGAPointFormValue(pnt, U, &sol_b[0]);
        PetscReal phi_b = PetscRealPart(sol_b[0]);

        PetscReal mob_b;
        Mobility(user, phi_b, &mob_b);

        const PetscReal *N0b;
        IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0b);

        PetscScalar (*Rb)[3] = (PetscScalar (*)[3])Re;
        PetscReal rwb = AxisymWeight(pnt, user);
        PetscReal flux;   /* 3*M*eps * dphi/dn */

        if (kind == FACE_WALL) {
            /* phi is clamped so the wall term is inert outside [0,1]; unclamped
             * it changes sign there and drives phi further out of range. */
            PetscReal phi_w = PetscMin(PetscMax(phi_b, 0.0), 1.0);
            PetscReal dh;
            WallH(phi_w, NULL, &dh, NULL);
            flux = 3.0 * mob_b * user->costhet * (dh / 6.0);
        } else {
            /* Open face: keep the boundary term of the integration by parts
             * with dphi/dn taken from the solution itself. */
            PetscScalar grad_b[3][dim];
            IGAPointFormGrad(pnt, U, &grad_b[0][0]);
            PetscReal dphi_dn = 0.0;
            for (l = 0; l < dim; l++) dphi_dn += PetscRealPart(grad_b[0][l]) * pnt->normal[l];
            flux = 3.0 * mob_b * eps * dphi_dn;
        }

        for (PetscInt a = 0; a < pnt->nen; a++) {
            Rb[a][0] = -rwb * flux * N0b[a];
            Rb[a][1] = 0.0;
            Rb[a][2] = 0.0;
        }
        return 0;
    }

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar phi   = sol[0],  phi_t  = sol_t[0];
    PetscScalar tem   = sol[1],  tem_t  = sol_t[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    PetscScalar grad_phi [dim];
    PetscScalar grad_tem [dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_phi [l] = grad_sol[0][l];
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar phi_a = 1.0 - phi;

    /* Out-of-bounds trial iterate: tell SNES so the line search backs off. */
    {
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (PetscRealPart(phi)   < lo || PetscRealPart(phi)   > hi ||
            PetscRealPart(phi_a) < lo || PetscRealPart(phi_a) > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;
    if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;
    if (phi_ac > 1.0) phi_ac = 1.0;

    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, phi_c,  &thcond,  NULL);
    HeatCap    (user, phi_c,  &cp,      NULL);
    Density    (user, phi_c,  &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, phi_c,  &mob_sub);

    PetscReal rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, NULL);

    PetscReal f1;
    DoubleWellDeriv(phi_c, &f1, NULL);
    PetscReal loc     = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal phi_aef = (phi_ac > air_lim) ? phi_ac : air_lim;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar (*R)[3] = (PetscScalar (*)[3])Re;
    PetscInt a, nen = pnt->nen;

    PetscReal rw = AxisymWeight(pnt, user);

    for (a = 0; a < nen; a++) {
        PetscReal gN_gphi  = 0.0;
        PetscReal gN_gtem  = 0.0;
        PetscReal gN_grhov = 0.0;
        for (l = 0; l < dim; l++) {
            gN_gphi  += N1[a][l] * PetscRealPart(grad_phi [l]);
            gN_gtem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            gN_grhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }

        /* -decouple_phase_change 1 zeroes every phase-change coupling. */
        const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

        R[a][0] = rw * ( N0[a] * phi_t
                + 3.0 * mob_sub * eps * gN_gphi
                + (3.0 * mob_sub / eps) * f1 * N0[a]
                - pc * (user->alph_sub / rho_ice) * loc
                  * (PetscRealPart(rhov) - rho_vs) * N0[a] );

        R[a][1] = rw * ( rho * cp * N0[a] * tem_t
                + user->xi_T * thcond * gN_gtem
                - pc * user->xi_T * rho_ice * lat_sub * phi_t * N0[a] );

        /* xi_v scales vapor diffusion and the mass-exchange source together. */
        R[a][2] = rw * ( phi_aef * N0[a] * rhov_t
                + user->xi_v * dif_vap * phi_aef * gN_grhov
                + pc * (user->xi_v * rho_ice - PetscRealPart(rhov))
                  * phi_t * N0[a] );
    }
    return 0;
}


PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    return Residual_A1(pnt, shift, V, t, U, Re, ctx);
}


/* J[a][i][b][j] = dR[a][i]/du[b][j] + shift * dR[a][i]/du_t[b][j] */
static PetscErrorCode Jacobian_A1(IGAPoint pnt,
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

    if (pnt->atboundary) {
        FaceKind kind = FaceKindAt(pnt, user);
        if (kind == FACE_NATURAL) return 0;

        PetscScalar sol_b[3];
        IGAPointFormValue(pnt, U, &sol_b[0]);
        PetscReal phi_b = PetscRealPart(sol_b[0]);

        PetscReal mob_b;
        Mobility(user, phi_b, &mob_b);

        const PetscReal *N0b, (*N1b)[dim];
        IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0b);
        IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1b);

        PetscInt nen_b = pnt->nen;
        PetscScalar (*Jb)[3][nen_b][3] = (PetscScalar (*)[3][nen_b][3])Je;
        PetscReal rwb = AxisymWeight(pnt, user);

        if (kind == FACE_WALL) {
            /* Same clamp as the residual; zero derivative outside [0,1]. */
            PetscReal phi_w = PetscMin(PetscMax(phi_b, 0.0), 1.0);
            const PetscReal dclamp = (phi_b > 0.0 && phi_b < 1.0) ? 1.0 : 0.0;
            PetscReal d2h;
            WallH(phi_w, NULL, NULL, &d2h);
            for (PetscInt a = 0; a < nen_b; a++)
                for (PetscInt b = 0; b < nen_b; b++)
                    Jb[a][0][b][0] += -rwb * 3.0 * mob_b * user->costhet
                                    * (d2h / 6.0) * dclamp * N0b[a] * N0b[b];
        } else {
            for (PetscInt a = 0; a < nen_b; a++)
                for (PetscInt b = 0; b < nen_b; b++) {
                    PetscReal dNb_dn = 0.0;
                    for (l = 0; l < dim; l++) dNb_dn += N1b[b][l] * pnt->normal[l];
                    Jb[a][0][b][0] += -rwb * 3.0 * mob_b * eps * N0b[a] * dNb_dn;
                }
        }
        return 0;
    }

    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscScalar phi   = sol[0],  phi_t  = sol_t[0];
    PetscScalar tem   = sol[1];
    PetscScalar rhov  = sol[2],  rhov_t = sol_t[2];
    PetscScalar grad_tem [dim], grad_rhov[dim];
    for (l = 0; l < dim; l++) {
        grad_tem [l] = grad_sol[1][l];
        grad_rhov[l] = grad_sol[2][l];
    }
    PetscScalar phi_a = 1.0 - phi;

    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;
    if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;
    if (phi_ac > 1.0) phi_ac = 1.0;

    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, phi_c,  &thcond,  NULL);
    HeatCap    (user, phi_c,  &cp,      NULL);
    Density    (user, phi_c,  &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, phi_c,  &mob_sub);

    PetscReal dthcond_dphi, d_dif_vap;
    ThermalCond(user, phi_c, NULL, &dthcond_dphi);
    VaporDiffus(user, tem,   NULL, &d_dif_vap);

    PetscReal rho_vs, d_rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, &d_rho_vs);

    PetscReal df1;
    DoubleWellDeriv(phi_c, NULL, &df1);

    PetscReal loc      = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal dloc_dph = 2.0 * phi_c * phi_ac * (phi_ac - phi_c);

    PetscBool phi_a_above_lim = (phi_ac > air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal phi_aef = phi_a_above_lim ? phi_ac : air_lim;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[3][nen][3] = (PetscScalar (*)[3][nen][3])Je;

    PetscReal rw = AxisymWeight(pnt, user);

    const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

    for (a = 0; a < nen; a++) {
        PetscReal gNa_gtem  = 0.0;
        PetscReal gNa_grhov = 0.0;
        for (l = 0; l < dim; l++) {
            gNa_gtem  += N1[a][l] * PetscRealPart(grad_tem [l]);
            gNa_grhov += N1[a][l] * PetscRealPart(grad_rhov[l]);
        }

        for (b = 0; b < nen; b++) {
            PetscReal NaNb   = N0[a] * N0[b];
            PetscReal gNagNb = 0.0;
            for (l = 0; l < dim; l++) gNagNb += N1[a][l] * N1[b][l];

            /* R_ice / phi */
            J[a][0][b][0] += rw * ( shift * NaNb
                           + 3.0 * mob_sub * eps * gNagNb
                           + (3.0 * mob_sub / eps) * df1 * NaNb
                           - pc * (user->alph_sub / rho_ice) * dloc_dph
                             * (PetscRealPart(rhov) - rho_vs) * NaNb );

            /* R_ice / T */
            J[a][0][b][1] += rw * pc * ( (user->alph_sub / rho_ice) * loc * d_rho_vs * NaNb );

            /* R_ice / rhov */
            J[a][0][b][2] -= rw * pc * ( (user->alph_sub / rho_ice) * loc * NaNb );

            /* R_tem / phi */
            J[a][1][b][0] += rw * ( user->xi_T * dthcond_dphi * gNa_gtem * N0[b]
                           - pc * user->xi_T * rho_ice * lat_sub * shift * NaNb );

            /* R_tem / T */
            J[a][1][b][1] += rw * ( shift * rho * cp * NaNb
                           + user->xi_T * thcond * gNagNb );

            /* R_vap / phi */
            if (phi_a_above_lim) {
                J[a][2][b][0] += rw * ( -NaNb * PetscRealPart(rhov_t)
                               - user->xi_v * dif_vap * N0[b] * gNa_grhov );
            }
            J[a][2][b][0] += rw * pc * ( (user->xi_v * rho_ice - PetscRealPart(rhov)) * shift * NaNb );

            /* R_vap / T */
            J[a][2][b][1] += rw * ( user->xi_v * d_dif_vap * phi_aef * gNa_grhov * N0[b] );

            /* R_vap / rhov */
            J[a][2][b][2] += rw * ( phi_aef * shift * NaNb
                           + user->xi_v * dif_vap * phi_aef * gNagNb
                           - PetscRealPart(phi_t) * NaNb );
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


/* Monitor integrals: phi, phi^2*phi_a^2, phi_a, T, rhov*phi_a. */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                           PetscScalar *S, void *ctx)
{
    PetscFunctionBegin;
    (void)n;
    AppCtx *user = (AppCtx *)ctx;
    PetscScalar sol[3];
    IGAPointFormValue(pnt, U, &sol[0]);

    PetscReal phi  = PetscRealPart(sol[0]);
    PetscReal tem  = PetscRealPart(sol[1]);
    PetscReal rhov = PetscRealPart(sol[2]);
    PetscReal phi_a = 1.0 - phi;

    /* Axisymmetric: full 2*pi*r measure, so the integrals are 3D volumes. */
    PetscReal rw = 1.0;
    if (user && user->axisym) {
        PetscReal xphys[3] = {0.0, 0.0, 0.0};
        IGAPointFormPoint(pnt, xphys);
        rw = 2.0 * PETSC_PI * xphys[1];
    }

    S[0] = rw * phi;
    S[1] = rw * phi*phi * phi_a*phi_a;
    S[2] = rw * phi_a;
    S[3] = rw * tem;
    S[4] = rw * rhov * phi_a;

    PetscFunctionReturn(0);
}
