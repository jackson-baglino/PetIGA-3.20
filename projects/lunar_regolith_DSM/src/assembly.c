#include "assembly.h"
#include "material_properties.h"

/* =========================================================================
 * Weak-form residuals for the 2-phase (ice / T / vapor) system.
 * phi_a = 1 - phi_i is algebraic; no sediment DOF.
 *
 * Notation: phi = phi_i,  M = mob_sub.
 *
 * f1(phi) = phi*(1-phi)*(1-2*phi)   [d/dphi of the double well (1/2)*phi^2*(1-phi)^2]
 * loc(phi) = phi^2*(1-phi)^2         [phase-change localization at interface]
 *
 * Allen-Cahn (ice):
 *   R_ice = N*phi_t  +  3*M*eps * grad_N.grad_phi  +  (3*M/eps)*f1 * N
 *         - (alph_sub/rho_ice) * loc * (rhov - rho_vs) * N
 *
 * Temperature:
 *   R_tem = rho*cp * N*T_t  +  k * grad_N.grad_T  -  rho_ice*L * phi_t * N
 *
 * Vapor (M&F 2024 eq. 26, temporal scaling):
 *   (1/xi_v) * d(phi_a*rhov)/dt = div(D_v*phi_a*grad_rhov) + rho_ice*phi_a_t
 * Multiplying through by xi_v and expanding the storage product gives
 *   R_vap = phi_a_eff * N*rhov_t  -  rhov * phi_t * N
 *         + xi_v*D_v*phi_a_eff * grad_N.grad_rhov
 *         + xi_v*rho_ice * phi_t * N
 *   where phi_a_eff = max(1-phi, air_lim)
 *   xi_v scales diffusion AND the rho_ice mass-exchange source TOGETHER:
 *   in the quasi-steady limit xi_v cancels between them, so the vapor field
 *   and interface velocity stay physical while the fast diffusion timescale
 *   is slowed by 1/xi_v (allows large dt). The -rhov*phi_t part belongs to
 *   the storage term and stays unscaled. The apparent ice->vapor mass
 *   "loss" of (1-xi_v)*Delta_m_vapor is ~rho_v/rho_i ~ 1e-6 relative —
 *   the approximation M&F accept for the temporal scaling.
 * ========================================================================= */

/* f1(phi) = phi*(1-phi)*(1-2*phi)  and  df1/dphi = 1 - 6*phi + 6*phi^2
 *
 * NOTE THE NORMALISATION: this is d/dphi of (1/2)*phi^2*(1-phi)^2, NOT of
 * phi^2*(1-phi)^2 (which would be twice this). The half-normalised well is
 * deliberate and must not be "corrected": it is what makes the equilibrium
 * profile of eps^2 phi'' = f1(phi) come out as
 *
 *     phi(x) = 0.5*(1 + tanh(x/(2*eps)))
 *
 * i.e. interface-thickness parameter W = eps exactly. That is the Karma-Plapp
 * convention the matched asymptotics assume, and it is the W that
 * lambda_sub = a1*eps/d0_sub in <project>_main.c is calibrated against.
 * Doubling f1 would shrink the real profile to W = eps/sqrt(2) while the
 * Gibbs-Thomson and kinetic calibration kept using eps, throwing d0 and beta
 * off by sqrt(2). Verified against run data: measured W/eps = 1.015. */
static void DoubleWellDeriv(PetscReal phi, PetscReal *f1, PetscReal *df1)
{
    if (f1)  *f1  = phi * (1.0 - phi) * (1.0 - 2.0 * phi);
    if (df1) *df1 = 1.0 - 6.0 * phi + 6.0 * phi * phi;
}

/* =========================================================================
 * Wall free energy: prescribed contact angle at a regolith boundary.
 *
 * The substrate is not a phase field here -- it IS the domain wall. Its
 * physics enters as a surface integral added to the free-energy functional,
 *
 *     F_wall = int_Gamma f_w(phi) dGamma
 *     f_w(phi) = gamma_as + (gamma_is - gamma_as)*h(phi)
 *              = gamma_as - gamma_ia*cos(theta)*h(phi)
 *
 * with the standard interpolation h(phi) = phi^2*(3-2*phi), so that
 * f_w(0) = gamma_as (wall|vapor) and f_w(1) = gamma_is (wall|ice), and
 * Young's equation gamma_ia*cos(theta) = gamma_as - gamma_is fixes theta.
 * h'(0) = h'(1) = 0, so the wall term is inert in both bulk phases and acts
 * only where the ice-vapor interface actually meets the wall -- it cannot
 * shift the bulk equilibria.
 *
 * Varying F_wall supplies the natural BC  dphi/dn = cos(theta)*phi(1-phi)/eps,
 * i.e. m.n = cos(theta) with m the interface normal into the ice and n the
 * outward wall normal: the contact angle measured through the ice.
 *
 * SCALING -- why neither eps nor gamma appears below.
 * The interior form is the gradient flow  dphi/dt = -(3M/eps) dFt/dphi  of the
 * DIMENSIONLESS functional Ft = int [ (eps^2/2)|grad phi|^2 + W(phi) ] with
 * W = (1/2)phi^2(1-phi)^2 (see DoubleWellDeriv above). Ft carries an
 * interfacial excess of exactly eps/6 per unit area, so physical energy is
 * (6*gamma_ia/eps)*Ft and the wall term enters Ft as ft_w = (eps/6/gamma_ia)*f_w.
 * The constant gamma_as drops out under variation and gamma_ia cancels against
 * the normalisation:
 *
 *     d ft_w / dphi = -(eps/6)*cos(theta)*h'(phi) = -eps*cos(theta)*phi*(1-phi)
 *
 * Multiplying by the residual's (3M/eps) prefactor cancels eps as well:
 *
 *     R_bnd[a][0] = -3*M*cos(theta)*phi*(1-phi)*N0[a]
 *
 * So the wall term depends only on M, cos(theta) and phi. cos(theta) is the
 * entire physical content -- do NOT "restore" a gamma or an eps here.
 *
 * Cross-check: demo/Metamorph.c writes this same boundary term as
 * -N0[a]*3*eps*mob*|grad phi|*costhet. On the equilibrium profile
 * |grad phi| = phi(1-phi)/eps, so the two agree identically -- but the h(phi)
 * form needs no 1/|grad phi| regularisation (Metamorph clamps |grad phi| to
 * 1e-5 to avoid a divide-by-zero) and its Jacobian is an exact mass-matrix
 * block rather than a normalised-gradient derivative.
 *
 * Sign: cos(theta) = +1 sheets ice flat onto the wall (theta = 0 deg);
 * cos(theta) = -1 balls it up to a point contact (theta = 180 deg).
 * ========================================================================= */

/* h(phi) = phi^2*(3-2*phi), h'(phi) = 6*phi*(1-phi), h''(phi) = 6*(1-2*phi). */
static void WallH(PetscReal phi, PetscReal *h, PetscReal *dh, PetscReal *d2h)
{
    if (h)   *h   = phi * phi * (3.0 - 2.0 * phi);
    if (dh)  *dh  = 6.0 * phi * (1.0 - phi);
    if (d2h) *d2h = 6.0 * (1.0 - 2.0 * phi);
}

/* Is this boundary quadrature point on a face flagged as regolith? */
static PetscBool WallPointActive(IGAPoint pnt, const AppCtx *user)
{
    PetscInt axis = -1, side = -1;
    if (!user->wall_any) return PETSC_FALSE;
    if (user->costhet == 0.0) return PETSC_FALSE;   /* theta = 90 deg: exact no-op */
    IGAPointAtBoundary(pnt, &axis, &side);
    if (axis < 0 || axis > 2 || side < 0 || side > 1) return PETSC_FALSE;
    return user->wall_face[axis][side];
}

/* Axisymmetric r-weight, shared by the interior and boundary forms so the two
 * cannot drift apart. dV = 2*pi*r dr dz; the constant 2*pi cancels in R = 0.
 * The azimuthal curvature mode is generated by this weight via the r-weighted
 * integration by parts (docs/axisymmetric_plan.md section 1b).
 * NOTE: -axisym combined with a wetting wall is untested. */
static PetscReal AxisymWeight(IGAPoint pnt, const AppCtx *user)
{
    PetscReal xphys[3] = {0.0, 0.0, 0.0};
    if (!user->axisym) return 1.0;
    IGAPointFormPoint(pnt, xphys);
    return xphys[1];
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
    PetscReal rho_ice = user->rho_ice;
    PetscReal lat_sub = user->lat_sub;
    PetscReal air_lim = user->air_lim;

    /* Wall (regolith) faces carry the prescribed-contact-angle surface term;
     * every other face keeps the natural Neumann dphi/dn = 0 it always had.
     * PetIGA visits an element once per enabled boundary face and once for its
     * interior, summing both into the same element vector, so this branch adds
     * to -- rather than replaces -- the volume form. Re is the per-point work
     * vector and IGAPointGetWorkVec has already zeroed it. */
    if (pnt->atboundary) {
        if (!WallPointActive(pnt, user)) return 0;

        PetscScalar sol_b[3];
        IGAPointFormValue(pnt, U, &sol_b[0]);
        PetscReal phi_b = PetscRealPart(sol_b[0]);

        /* Raw phi, deliberately NOT the clamped phi_c used for material
         * properties below: phi*(1-phi) is a polynomial, and leaving it
         * unclamped keeps the analytic Jacobian exactly consistent with this
         * residual, which is what -snes_test_jacobian checks. A trial iterate
         * far outside [0,1] is caught by the interior branch's domain guard on
         * the same element. */
        PetscReal dh;
        WallH(phi_b, NULL, &dh, NULL);

        PetscReal mob_b;
        Mobility(user, phi_b, &mob_b);

        const PetscReal *N0b;
        IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0b);

        PetscScalar (*Rb)[3] = (PetscScalar (*)[3])Re;
        PetscReal rwb = AxisymWeight(pnt, user);

        /* -(3M/eps)*eps*cos(theta)*h'(phi)/6*N = -3*M*cos(theta)*phi(1-phi)*N */
        for (PetscInt a = 0; a < pnt->nen; a++) {
            Rb[a][0] = -rwb * 3.0 * mob_b * user->costhet * (dh / 6.0) * N0b[a];
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

    /* SNES domain-error catch: if a trial Newton iterate has phi out of the
     * configured bounds, signal an invalid state so the line search backs off. */
    {
        PetscReal lo = user->phase_lo, hi = user->phase_hi;
        if (PetscRealPart(phi)   < lo || PetscRealPart(phi)   > hi ||
            PetscRealPart(phi_a) < lo || PetscRealPart(phi_a) > hi) {
            if (user->snes) SNESSetFunctionDomainError(user->snes);
            return 0;
        }
    }

    /* Clamped copies for material-property evaluation (avoids negative inputs). */
    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;
    if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;
    if (phi_ac > 1.0) phi_ac = 1.0;

    /* Material properties. */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, phi_c,  &thcond,  NULL);
    HeatCap    (user, phi_c,  &cp,      NULL);
    Density    (user, phi_c,  &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, phi_c,  &mob_sub);

    /* Flat-interface saturation vapor density. */
    PetscReal rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, NULL);

    /* Double-well derivative and localization. */
    PetscReal f1;
    DoubleWellDeriv(phi_c, &f1, NULL);
    PetscReal loc     = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal phi_aef = (phi_ac > air_lim) ? phi_ac : air_lim;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscScalar (*R)[3] = (PetscScalar (*)[3])Re;
    PetscInt a, nen = pnt->nen;

    /* Axisymmetric r-z mode: weight the integrand by the quadrature point's
     * radial coordinate (dV = 2*pi*r dr dz; the constant 2*pi cancels in
     * R = 0). The azimuthal curvature mode is generated by this weight via
     * the r-weighted integration by parts — the |grad phi|^2 term itself is
     * untouched (docs/axisymmetric_plan.md section 1b). */
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

        /* -decouple_phase_change 1: zero every phase-change coupling — the
         * sublimation source in R_ice AND the ice_t-driven sources in
         * R_tem (latent heat) / R_vap (mass exchange) — leaving pure
         * Allen-Cahn curvature relaxation with passive heat/vapor
         * diffusion. Demo/diagnostic mode. */
        const PetscReal pc = user->decouple_phase_change ? 0.0 : 1.0;

        R[a][0] = rw * ( N0[a] * phi_t
                + 3.0 * mob_sub * eps * gN_gphi
                + (3.0 * mob_sub / eps) * f1 * N0[a]
                - pc * (user->alph_sub / rho_ice) * loc
                  * (PetscRealPart(rhov) - rho_vs) * N0[a] );

        R[a][1] = rw * ( rho * cp * N0[a] * tem_t                       /* storage */
                + user->xi_T * thcond * gN_gtem                        /* conduction */
                - pc * user->xi_T * rho_ice * lat_sub * phi_t * N0[a] ); /* latent heat */

        R[a][2] = rw * ( phi_aef * N0[a] * rhov_t                      /* vapor storage */
                + user->xi_v * dif_vap * phi_aef * gN_grhov            /* vapor diffusion (xi_v-scaled) */
                + pc * (user->xi_v * rho_ice - PetscRealPart(rhov))
                  * phi_t * N0[a] );                                    /* xi_v*source - storage cross-term */
    }
    return 0;
}


/* =========================================================================
 * Residual dispatcher
 * ========================================================================= */
PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U,
                        PetscScalar *Re, void *ctx)
{
    return Residual_A1(pnt, shift, V, t, U, Re, ctx);
}


/* =========================================================================
 * Jacobian.  J[a][i][b][j] = dR[a][i]/du[b][j] + shift * dR[a][i]/du_t[b][j]
 * DOF layout: 0=phi(ice), 1=T, 2=rhov.
 * ========================================================================= */
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

    /* Exact Jacobian of the wall term in Residual_A1's boundary branch:
     *   d/dphi [ -3*M*cos(theta)*phi*(1-phi) ] = -3*M*cos(theta)*(1-2*phi)
     * i.e. -3*M*cos(theta)*h''(phi)/6. No shift term: the wall energy has no
     * phi_t dependence. */
    if (pnt->atboundary) {
        if (!WallPointActive(pnt, user)) return 0;

        PetscScalar sol_b[3];
        IGAPointFormValue(pnt, U, &sol_b[0]);
        PetscReal phi_b = PetscRealPart(sol_b[0]);

        PetscReal d2h;
        WallH(phi_b, NULL, NULL, &d2h);

        PetscReal mob_b;
        Mobility(user, phi_b, &mob_b);

        const PetscReal *N0b;
        IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0b);

        PetscInt nen_b = pnt->nen;
        PetscScalar (*Jb)[3][nen_b][3] = (PetscScalar (*)[3][nen_b][3])Je;
        PetscReal rwb = AxisymWeight(pnt, user);

        for (PetscInt a = 0; a < nen_b; a++)
            for (PetscInt b = 0; b < nen_b; b++)
                Jb[a][0][b][0] += -rwb * 3.0 * mob_b * user->costhet
                                * (d2h / 6.0) * N0b[a] * N0b[b];
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

    /* Clamped copies for material-property evaluation. */
    PetscReal phi_c  = PetscRealPart(phi);
    PetscReal phi_ac = PetscRealPart(phi_a);
    if (phi_c  < 0.0) phi_c  = 0.0;
    if (phi_c  > 1.0) phi_c  = 1.0;
    if (phi_ac < 0.0) phi_ac = 0.0;
    if (phi_ac > 1.0) phi_ac = 1.0;

    /* Material properties and their derivatives. */
    PetscReal thcond, cp, rho, dif_vap, mob_sub;
    ThermalCond(user, phi_c,  &thcond,  NULL);
    HeatCap    (user, phi_c,  &cp,      NULL);
    Density    (user, phi_c,  &rho,     NULL);
    VaporDiffus(user, tem,    &dif_vap, NULL);
    Mobility   (user, phi_c,  &mob_sub);

    PetscReal dthcond_dphi, d_dif_vap;
    ThermalCond(user, phi_c, NULL, &dthcond_dphi);
    VaporDiffus(user, tem,   NULL, &d_dif_vap);

    /* Flat-interface saturation vapor density and its T-derivative. */
    PetscReal rho_vs, d_rho_vs;
    RhoVS_I(user, PetscRealPart(tem), &rho_vs, &d_rho_vs);

    /* Double-well derivative and its phi-derivative. */
    PetscReal df1;
    DoubleWellDeriv(phi_c, NULL, &df1);

    /* Localization loc = phi^2*(1-phi)^2 and dloc/dphi. */
    PetscReal loc      = phi_c * phi_c * phi_ac * phi_ac;
    PetscReal dloc_dph = 2.0 * phi_c * phi_ac * (phi_ac - phi_c);

    /* phi_a floor for vapor storage/diffusion. */
    PetscBool phi_a_above_lim = (phi_ac > air_lim) ? PETSC_TRUE : PETSC_FALSE;
    PetscReal phi_aef = phi_a_above_lim ? phi_ac : air_lim;

    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    PetscInt a, b, nen = pnt->nen;
    PetscScalar (*J)[3][nen][3] = (PetscScalar (*)[3][nen][3])Je;

    /* Axisymmetric r-weight — must match Residual exactly (see comment
     * there); applied to every block below. */
    PetscReal rw = AxisymWeight(pnt, user);

    /* Must mirror the pc factor in Residual (see comment there). */
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

            /* ============ R_ice / phi ============ */
            J[a][0][b][0] += rw * ( shift * NaNb
                           + 3.0 * mob_sub * eps * gNagNb
                           + (3.0 * mob_sub / eps) * df1 * NaNb
                           - pc * (user->alph_sub / rho_ice) * dloc_dph
                             * (PetscRealPart(rhov) - rho_vs) * NaNb );

            /* ============ R_ice / T ============ */
            J[a][0][b][1] += rw * pc * ( (user->alph_sub / rho_ice) * loc * d_rho_vs * NaNb );

            /* ============ R_ice / rhov ============ */
            J[a][0][b][2] -= rw * pc * ( (user->alph_sub / rho_ice) * loc * NaNb );

            /* ============ R_tem / phi ============ */
            J[a][1][b][0] += rw * ( user->xi_T * dthcond_dphi * gNa_gtem * N0[b]
                           - pc * user->xi_T * rho_ice * lat_sub * shift * NaNb );

            /* ============ R_tem / T ============ */
            J[a][1][b][1] += rw * ( shift * rho * cp * NaNb
                           + user->xi_T * thcond * gNagNb );

            /* ============ R_vap / phi ============ */
            if (phi_a_above_lim) {
                J[a][2][b][0] += rw * ( -NaNb * PetscRealPart(rhov_t)
                               - user->xi_v * dif_vap * N0[b] * gNa_grhov );
            }
            /* Exact Jacobian of the xi_v-scaled residual: the source coefficient
             * in R[a][2] is (xi_v*rho_ice - rhov), so the vap/ice coupling is
             * O(xi_v*rho_ice) ~ 1 instead of O(rho_ice) ~ 1e3 — well-conditioned. */
            J[a][2][b][0] += rw * pc * ( (user->xi_v * rho_ice - PetscRealPart(rhov)) * shift * NaNb );

            /* ============ R_vap / T ============ */
            J[a][2][b][1] += rw * ( user->xi_v * d_dif_vap * phi_aef * gNa_grhov * N0[b] );

            /* ============ R_vap / rhov ============ */
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


/* =========================================================================
 * Per-element scalar integrals for the monitor table.
 *   S[0] = phi_i (ice volume fraction)
 *   S[1] = phi_i^2 * phi_a^2  (ice-air interface measure)
 *   S[2] = phi_a (air volume fraction)
 *   S[3] = T
 *   S[4] = rhov * phi_a
 * ========================================================================= */
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

    /* Axisymmetric: include the FULL 2*pi*r measure so the reported
     * integrals are true 3D volumes (TOT_ICE in m^3, etc.) and mass
     * conservation checks remain meaningful. */
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
