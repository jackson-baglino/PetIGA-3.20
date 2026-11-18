#include "assembly.h"
#include "material_properties.h"

/* The derived model (D1): residual and Jacobian written from the blind rederivation by beth-mathematician,
 * 2026-11-18, not from Residual_A1. Selected with -residual_model derived; default off.
 *
 * Strong form (runs/model-development/beth-mathematician/<run>/answer.json, result_expr):
 *   phi (T07, 20261118-124517):
 *     phi_t + M [ 6 gamma/eps phi(1-phi)(1-2phi) - 6 gamma eps lap(phi) ]
 *           - 30 M rho_i R_v T phi^2 (1-phi)^2 (rho_v - rho_vs(T)) / rho_vs(T) = 0
 *   T (T04, 20261118-120651):
 *     C(phi) T_t - div(k(phi) grad T) - rho_i L phi_t = 0,   C = phi rho_i c_i + (1-phi) rho_a c_a
 *   rho_v (T05, 20261118-123403):
 *     (1-phi) rho_v,t - div(D(T) (1-phi) grad rho_v) + (rho_i - rho_v) phi_t = 0
 *
 * Choices that are NOT from the derivation, kept identical to A1 so the two can be compared:
 *   - the xi_T and xi_v time-scaling factors (Moure and Fu's acceleration), placed as in A1;
 *   - the air-fraction floor max(1 - phi, air_lim) in the vapor row, and the clamp of phi to [0, 1] in the
 *     coefficients;
 *   - the constitutive laws k(phi), D(T), rho_vs(T) and the constant L of material_properties.c;
 *   - rho_v c_v in C is left out (about 1e-6 of rho_a c_a), and T04's L(T) is taken constant.
 * Options: -derived_cp_mix 0 uses A1's rho(phi)*cp(phi); -derived_frozen 1 takes T/rho_vs(T) in the coupling
 * at temp0, as A1's constant alph_sub does. With both, and M = mob_sub/(2 gamma), D1 equals A1 term by term
 * (studies/derived_model/README.md).
 *
 * Boundaries: natural only. lunar_main.c refuses wall faces, open faces and -axisym with this model.
 * T is in Celsius in the solver, as in A1; the coupling uses T + 273.15. */

#define R_V (K_BOLTZ / M_H2O_KG) /* specific gas constant of water vapor [J kg^-1 K^-1] */

typedef struct {
    PetscReal phi, phi_t, tem, tem_t, rhov, rhov_t;
    PetscReal phi_c, phi_ac, phi_aef;
    PetscBool above_lim, phi_inside;
} PointState;

/* Heat capacity per unit volume and its derivative in phi. */
static void HeatCapVol(AppCtx *user, PetscReal phi_c, PetscReal *C, PetscReal *dC)
{
    if (user->derived_cp_mix) {
        if (C)  *C  = phi_c * user->rho_ice * user->cp_ice + (1.0 - phi_c) * user->rho_air * user->cp_air;
        if (dC) *dC = user->rho_ice * user->cp_ice - user->rho_air * user->cp_air;
    } else {
        PetscScalar cp, dcp, rho, drho;
        HeatCap(user, phi_c, &cp, &dcp);
        Density(user, phi_c, &rho, &drho);
        if (C)  *C  = PetscRealPart(rho * cp);
        if (dC) *dC = PetscRealPart(drho * cp + rho * dcp);
    }
}

/* Coupling prefactor K = 30 M rho_i R_v T / rho_vs(T) [m^3 kg^-1 s^-1] and dK/dT. */
static void Coupling(AppCtx *user, PetscReal tem, PetscReal *K, PetscReal *dK)
{
    PetscReal TK, rvs, drvs = 0.0;
    if (user->derived_frozen) {
        TK  = user->temp0 + 273.15;
        rvs = user->rhovs0;
    } else {
        PetscScalar r, dr;
        TK = tem + 273.15;
        RhoVS_I(user, tem, &r, &dr);
        rvs = PetscRealPart(r); drvs = PetscRealPart(dr);
    }
    PetscReal k = 30.0 * user->derived_M * user->rho_ice * R_V * TK / rvs;
    if (K)  *K  = k;
    if (dK) *dK = user->derived_frozen ? 0.0 : k * (1.0 / TK - drvs / rvs);
}

PetscErrorCode Residual_D1(IGAPoint pnt,
                           PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U,
                           PetscScalar *Re, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;
    (void)shift; (void)t;
    if (pnt->atboundary) return 0;

    PetscInt dim = user->dim;
    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscReal phi = PetscRealPart(sol[0]), phi_t = PetscRealPart(sol_t[0]);
    PetscReal tem = PetscRealPart(sol[1]), tem_t = PetscRealPart(sol_t[1]);
    PetscReal rhov = PetscRealPart(sol[2]), rhov_t = PetscRealPart(sol_t[2]);

    if (phi < user->phase_lo || phi > user->phase_hi || 1.0 - phi < user->phase_lo || 1.0 - phi > user->phase_hi) {
        if (user->snes) SNESSetFunctionDomainError(user->snes);
        return 0;
    }
    PetscReal phi_c   = PetscMin(PetscMax(phi, 0.0), 1.0);
    PetscReal phi_ac  = 1.0 - phi_c;
    PetscReal phi_aef = PetscMax(phi_ac, user->air_lim);

    const PetscReal eps = user->eps, gam = user->gamma_ia_bulk, M = user->derived_M;
    const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

    PetscScalar k, D, rvs;
    ThermalCond(user, phi_c, &k, NULL);
    VaporDiffus(user, tem, &D, NULL);
    RhoVS_I(user, tem, &rvs, NULL);
    PetscReal C, K;
    HeatCapVol(user, phi_c, &C, NULL);
    Coupling(user, tem, &K, NULL);

    PetscReal dwell = phi_c * (1.0 - phi_c) * (1.0 - 2.0 * phi_c);   /* d/dphi of phi^2(1-phi)^2 / 2 */
    PetscReal hp30  = phi_c * phi_c * phi_ac * phi_ac;               /* h'(phi)/30 */
    PetscReal drive = K * hp30 * (rhov - PetscRealPart(rvs));

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    PetscScalar (*R)[3] = (PetscScalar (*)[3])Re;

    for (PetscInt a = 0; a < pnt->nen; a++) {
        PetscReal gphi = 0.0, gtem = 0.0, grhov = 0.0;
        for (PetscInt l = 0; l < dim; l++) {
            gphi  += N1[a][l] * PetscRealPart(grad_sol[0][l]);
            gtem  += N1[a][l] * PetscRealPart(grad_sol[1][l]);
            grhov += N1[a][l] * PetscRealPart(grad_sol[2][l]);
        }
        R[a][0] = N0[a] * phi_t + 6.0 * gam * M * eps * gphi + (6.0 * gam * M / eps) * dwell * N0[a]
                - pc * drive * N0[a];
        R[a][1] = C * N0[a] * tem_t + user->xi_T * PetscRealPart(k) * gtem
                - pc * user->xi_T * user->rho_ice * user->lat_sub * phi_t * N0[a];
        R[a][2] = phi_aef * N0[a] * rhov_t + user->xi_v * PetscRealPart(D) * phi_aef * grhov
                + pc * (user->xi_v * user->rho_ice - rhov) * phi_t * N0[a];
    }
    return 0;
}

/* J[a][i][b][j] = dR[a][i]/du[b][j] + shift dR[a][i]/du_t[b][j], every term of Residual_D1 (inside the
 * clamps; the clamp and floor switches have zero derivative, as in A1). */
PetscErrorCode Jacobian_D1(IGAPoint pnt,
                           PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U,
                           PetscScalar *Je, void *ctx)
{
    AppCtx *user = (AppCtx*)ctx;
    (void)t;
    if (pnt->atboundary) return 0;

    PetscInt dim = user->dim;
    PetscScalar sol_t[3], sol[3], grad_sol[3][dim];
    IGAPointFormValue(pnt, V, &sol_t[0]);
    IGAPointFormValue(pnt, U, &sol[0]);
    IGAPointFormGrad (pnt, U, &grad_sol[0][0]);

    PetscReal phi = PetscRealPart(sol[0]), phi_t = PetscRealPart(sol_t[0]);
    PetscReal tem = PetscRealPart(sol[1]), tem_t = PetscRealPart(sol_t[1]);
    PetscReal rhov = PetscRealPart(sol[2]), rhov_t = PetscRealPart(sol_t[2]);

    PetscReal inside  = (phi > 0.0 && phi < 1.0) ? 1.0 : 0.0;
    PetscReal phi_c   = PetscMin(PetscMax(phi, 0.0), 1.0);
    PetscReal phi_ac  = 1.0 - phi_c;
    PetscBool above   = (phi_ac > user->air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal phi_aef = above ? phi_ac : user->air_lim;

    const PetscReal eps = user->eps, gam = user->gamma_ia_bulk, M = user->derived_M;
    const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

    PetscScalar k, dk, D, dD, rvs, drvs;
    ThermalCond(user, phi_c, &k, &dk);
    VaporDiffus(user, tem, &D, &dD);
    RhoVS_I(user, tem, &rvs, &drvs);
    PetscReal C, dC, K, dK;
    HeatCapVol(user, phi_c, &C, &dC);
    Coupling(user, tem, &K, &dK);

    PetscReal ddwell = (1.0 - 6.0 * phi_c + 6.0 * phi_c * phi_c) * inside;
    PetscReal hp30   = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal dhp30  = 2.0 * phi_c * phi_ac * (phi_ac - phi_c) * inside;
    PetscReal sat    = rhov - PetscRealPart(rvs);
    dk *= inside; dC *= inside;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    PetscInt nen = pnt->nen;
    PetscScalar (*J)[3][nen][3] = (PetscScalar (*)[3][nen][3])Je;

    for (PetscInt a = 0; a < nen; a++) {
        PetscReal gtem = 0.0, grhov = 0.0;
        for (PetscInt l = 0; l < dim; l++) {
            gtem  += N1[a][l] * PetscRealPart(grad_sol[1][l]);
            grhov += N1[a][l] * PetscRealPart(grad_sol[2][l]);
        }
        for (PetscInt b = 0; b < nen; b++) {
            PetscReal NN = N0[a] * N0[b], GG = 0.0;
            for (PetscInt l = 0; l < dim; l++) GG += N1[a][l] * N1[b][l];

            J[a][0][b][0] += shift * NN + 6.0 * gam * M * eps * GG + (6.0 * gam * M / eps) * ddwell * NN
                           - pc * K * dhp30 * sat * NN;
            J[a][0][b][1] += -pc * (dK * sat - K * PetscRealPart(drvs)) * hp30 * NN;
            J[a][0][b][2] += -pc * K * hp30 * NN;

            J[a][1][b][0] += dC * tem_t * NN + user->xi_T * PetscRealPart(dk) * gtem * N0[b]
                           - pc * user->xi_T * user->rho_ice * user->lat_sub * shift * NN;
            J[a][1][b][1] += shift * C * NN + user->xi_T * PetscRealPart(k) * GG;

            if (above)
                J[a][2][b][0] += -NN * rhov_t - user->xi_v * PetscRealPart(D) * N0[b] * grhov;
            J[a][2][b][0] += pc * (user->xi_v * user->rho_ice - rhov) * shift * NN;
            J[a][2][b][1] += user->xi_v * PetscRealPart(dD) * phi_aef * grhov * N0[b];
            J[a][2][b][2] += phi_aef * shift * NN + user->xi_v * PetscRealPart(D) * phi_aef * GG
                           - pc * phi_t * NN;
        }
    }
    return 0;
}
