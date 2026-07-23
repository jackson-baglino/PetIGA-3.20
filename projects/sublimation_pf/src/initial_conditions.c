#include "initial_conditions.h"
#include "material_properties.h"


/* =========================================================================
 * GrevilleAbscissae
 *
 * Fills `g[0..n-1]` with the Greville abscissae (parametric DOF locations,
 * in [0,1]) of the given IGA axis: g_i = mean(U[i+1..i+p]) for the axis's
 * knot vector U and degree p. For an open-uniform knot vector this reduces
 * to g_i = i/N (degree 1) but is non-uniform for degree >= 2 -- this is
 * the degree-agnostic replacement for the `i/(mx+per)` index mapping used
 * by IC functions on -geom_file domains.
 *
 * Caller must PetscFree(*g).
 * =========================================================================*/
static PetscErrorCode GrevilleAbscissae(IGA iga, PetscInt dir, PetscReal **g, PetscInt *n)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    IGAAxis axis;
    PetscInt p, m;
    PetscReal *U;
    ierr = IGAGetAxis(iga, dir, &axis);     CHKERRQ(ierr);
    ierr = IGAAxisGetDegree(axis, &p);      CHKERRQ(ierr);
    ierr = IGAAxisGetKnots(axis, &m, &U);   CHKERRQ(ierr);

    /* Normalize to [0,1]: a -geom_file axis's knots already span [0,1] (the
     * igakit builder's convention), but the default (-Nx/-Ny, no -geom_file)
     * axis built by IGAAxisInitUniform() spans [0,Lx]/[0,Ly] directly -- callers
     * (FormInitialMultiGrains2D) multiply the returned abscissae by Lx/Ly
     * themselves, so this must always return [0,1] regardless of which path
     * built the axis, or physical coordinates get scaled by Lx/Ly twice. */
    PetscReal Ui = U[0], Uf = U[m];
    PetscInt nb = m - p; /* number of basis functions / DOFs along this axis */
    PetscReal *greville;
    ierr = PetscMalloc1(nb, &greville); CHKERRQ(ierr);
    for (PetscInt i = 0; i < nb; i++) {
        PetscReal sum = 0.0;
        for (PetscInt k = 1; k <= p; k++) sum += U[i + k];
        greville[i] = (sum / (PetscReal)p - Ui) / (Uf - Ui);
    }

    *g = greville;
    *n = nb;
    PetscFunctionReturn(0);
}


/**
 * @brief 1D initial condition: a centered ice slab surrounded by air.
 *
 * Sets up a diffuse-interface ice slab occupying [x_lo, x_hi] = [0.35*Lx, 0.65*Lx],
 * with air filling the rest of the domain. Temperature and vapor density are
 * initialized from user->temp0, grad_temp0, hum0.
 *
 * Two geometric variants are selected via user->flag_tIC:
 *   flag_tIC == 0  ->  centered slab (ice in [0.35 Lx, 0.65 Lx])
 *   flag_tIC == 2  ->  flat interface (ice in [0, 0.5 Lx], air in [0.5 Lx, Lx])
 */
PetscErrorCode FormInitialCondition1D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD,
                "--- INITIAL CONDITIONS (1D, dim=%d) ---\n", user->dim);

    user->n_act = 0;

    /* Build node DM for the primary field vector */
    DM            da;
    Field        *u;   /* 1D: plain pointer, not pointer-to-pointer */
    DMDALocalInfo info;

    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);            CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);           CHKERRQ(ierr);

    const PetscReal Lx  = user->Lx;
    const PetscReal eps = user->eps;

    /* Slab extents — choose variant based on flag_tIC */
    PetscReal x_lo, x_hi;
    if (user->flag_tIC == 2) {
        /* Flat interface: left half is ice */
        x_lo = 0.0;
        x_hi = 0.5 * Lx;
    } else {
        /* Default: centered slab */
        x_lo = 0.35 * Lx;
        x_hi = 0.65 * Lx;
    }

    PetscInt k = -1;
    if (user->periodic == 1) k = user->p - 1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + k);

        /* Diffuse slab: tanh transition at x_lo and x_hi */
        PetscReal ice = 0.5 * (tanh(0.5 * (x - x_lo) / eps)
                              - tanh(0.5 * (x - x_hi) / eps));
        ice = PetscMin(PetscMax(ice, 0.0), 1.0);

        u[i].ice = ice;
        u[i].tem = user->temp0 + user->grad_temp0[0] * (x - 0.5 * Lx);

        PetscScalar rho_vs_loc, temp_loc = u[i].tem;
        RhoVS_I(user, temp_loc, &rho_vs_loc, NULL);
        { PetscReal _pa = PetscMax(0.0, 1.0 - ice);
          u[i].rhov = rho_vs_loc * (user->hum0 * _pa + (1.0 - _pa)); }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/* -------------------------------------------------------------------------
 * 2D Ice Slab Initial Condition
 *   - Ice: circular blob of radius Ly, air fills the rest
 *   - Temperature: temp0 + grad_temp0 * (r - 0.5*L)
 *   - Vapor: hum0 * rho_vs(T)
 * -------------------------------------------------------------------------*/
PetscErrorCode FormInitialIceSlab2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD,
                "--- INITIAL CONDITIONS (2D Ice Slab, R=Ly) ---\n");

    const PetscReal Lx  = user->Lx;
    const PetscReal Ly  = user->Ly;
    const PetscReal eps = user->eps;

    /* Circular ice blob; air fills the rest of the domain. */
    const PetscReal R       = Ly;
    const PetscReal x_ice_c = 0.40 * Lx;
    const PetscReal y_ice_c = 0.5 * Ly;

    user->n_act = 0;

    DM            da;
    Field       **u;
    DMDALocalInfo info;

    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);            CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);           CHKERRQ(ierr);

    PetscInt per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {

            PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
            PetscReal y = Ly * (PetscReal)j / (PetscReal)(info.my + per);

            /* Signed distance from the ice surface (negative inside ice) */
            PetscReal dx       = x - x_ice_c;
            PetscReal dy       = y - y_ice_c;
            PetscReal dist_ice = PetscSqrtReal(dx*dx + dy*dy) - R;

            /* Ice: interior of the circle */
            PetscReal ice = 0.5 - 0.5 * PetscTanhReal(0.5 * dist_ice / eps);
            ice = PetscMin(PetscMax(ice, 0.0), 1.0);

            PetscReal tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            PetscScalar rho_vs_loc;
            RhoVS_I(user, tem, &rho_vs_loc, NULL);

            PetscReal air = PetscMax(0.0, 1.0 - ice);

            u[j][i].ice  = ice;
            u[j][i].tem  = tem;
            u[j][i].rhov = user->hum0 * rho_vs_loc * air + rho_vs_loc * (1.0 - air);
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialSingleIceGrain2D
 *
 * Single pure ice circle (no sediment core) centered in the domain.
 * Tanh profile with half-width eps.
 *
 * Parameters: user->RCice (grain radius), Lx, Ly, eps, temp0, hum0.
 * =========================================================================*/
PetscErrorCode FormInitialSingleIceGrain2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx    = user->Lx;
    const PetscReal Ly        = user->Ly;
    const PetscReal eps       = user->eps;
    const PetscReal RCice     = user->RCice;
    /* Equilibrium logistic profile: phi = 1/(1+exp(-(R-d)/eps))
     * = 0.5 - 0.5*tanh(0.5*(d-R)/eps), so tc = 0.5/eps. Initializing at the
     * model's own equilibrium width (1%-99% band = 9.2*eps) removes the
     * early width-relaxation transient. The old tc = 1/(sqrt(2)*0.75*eps)
     * was 1.89x steeper — a leftover from the removed eps_model=0.75*eps
     * residual scaling — and made every run start with the IC ~7 cells wide
     * relaxing to the equilibrium 13 cells over the first ~60 steps. */
    const PetscReal tc        = 0.5 / eps;
    const PetscReal cx        = 0.5 * Lx;
    const PetscReal cy        = 0.5 * Ly;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D single ice grain) ---\n"
        "  centre = (%.4e, %.4e) m,  RCice = %.4e m\n",
        cx, cy, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice must be > 0 (got %.2e)", RCice);
    if (RCice >= 0.5 * PetscMin(Lx, Ly))
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Grain radius %.2e exceeds half the domain — increase domain or reduce RCice", RCice);

    DM da;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da, U, &u); CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

    PetscInt per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
            PetscReal y = Ly * (PetscReal)j / (PetscReal)(info.my + per);

            PetscReal dist = PetscSqrtReal(SQ(x - cx) + SQ(y - cy));
            PetscReal ice  = 0.5 - 0.5 * PetscTanhReal(tc * (dist - RCice));
            ice = PetscMin(PetscMax(ice, 0.0), 1.0);

            u[j][i].ice = ice;
            u[j][i].tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            PetscScalar rho_vs, temp_loc = u[j][i].tem;
            RhoVS_I(user, temp_loc, &rho_vs, NULL);
            { PetscReal _pa = PetscMax(0.0, 1.0 - ice);
              u[j][i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialSingleIceGrain1D
 *
 * 1D cross-section through the centre of a single ice grain.
 * Ice block centered at x = 0.5*Lx with radius RCice.
 * =========================================================================*/
PetscErrorCode FormInitialSingleIceGrain1D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx        = user->Lx;
    const PetscReal eps       = user->eps;
    const PetscReal RCice     = user->RCice;
    /* Equilibrium logistic profile: phi = 1/(1+exp(-(R-d)/eps))
     * = 0.5 - 0.5*tanh(0.5*(d-R)/eps), so tc = 0.5/eps. Initializing at the
     * model's own equilibrium width (1%-99% band = 9.2*eps) removes the
     * early width-relaxation transient. The old tc = 1/(sqrt(2)*0.75*eps)
     * was 1.89x steeper — a leftover from the removed eps_model=0.75*eps
     * residual scaling — and made every run start with the IC ~7 cells wide
     * relaxing to the equilibrium 13 cells over the first ~60 steps. */
    const PetscReal tc        = 0.5 / eps;
    const PetscReal cx        = 0.5 * Lx;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (1D single ice grain) ---\n"
        "  centre = %.4e m,  RCice = %.4e m\n", cx, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice must be > 0 (got %.2e)", RCice);
    if (RCice >= 0.5 * Lx)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Grain radius %.2e exceeds half the domain — increase Lx or reduce RCice", RCice);

    DM da;
    Field *u;
    DMDALocalInfo info;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);            CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);           CHKERRQ(ierr);

    PetscInt k = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        PetscReal x    = Lx * (PetscReal)i / (PetscReal)(info.mx + k);
        PetscReal dist = PetscAbsReal(x - cx);
        PetscReal ice  = 0.5 - 0.5 * PetscTanhReal(tc * (dist - RCice));
        ice = PetscMin(PetscMax(ice, 0.0), 1.0);

        u[i].ice = ice;
        u[i].tem = user->temp0 + user->grad_temp0[0] * (x - 0.5 * Lx);

        PetscScalar rho_vs, temp_loc = u[i].tem;
        RhoVS_I(user, temp_loc, &rho_vs, NULL);
        { PetscReal _pa = PetscMax(0.0, 1.0 - ice);
          u[i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialTwoIceGrainsBoundary2D
 *
 * Two ice semicircles centered ON the x=0 and x=Lx boundaries, vertically
 * centered (y = 0.5*Ly). Grain 0 (radius RCice0, smaller) sits at x=0;
 * grain 1 (radius RCice1, larger) sits at x=Lx. Each is a tanh distance
 * profile; the two are summed and clamped to [0,1].
 *
 * Centering each grain exactly on its boundary makes the profile symmetric
 * under reflection across that boundary, so it is automatically consistent
 * with the natural-Neumann (zero-flux) BC on phi_i. The shared vapor pore
 * spans the full domain width Lx, giving an Ostwald-ripening timescale
 * test between the two grains of different curvature.
 *
 * Physical coordinates: for a plain rectangular domain (geom_bump_R == 0)
 * node (i,j) maps to (Lx*i/mx, Ly*j/my) as usual. For the -geom_file
 * sediment-grain geometry (build_geometry_sediment_grain.py), the bottom
 * edge is raised by a C-infinity bump y=g(x) and the surface is ruled to
 * the flat top edge, so y_phys = g(x) + (j/my)*(Ly - g(x)); x_phys = Lx*i/mx
 * exactly by construction (both curves share the u<->x parametrization).
 * geom_bump_R must match that script's R_sed (bump half-width == height)
 * for the IC to align with the actual geometry.
 * =========================================================================*/

/* C-infinity bump g(x) = height*exp(1 - 1/(1-t^2)) for |t|<1, t=(x-center)/R;
 * 0 outside -- must match build_geometry_sediment_grain.py's _bump(). */
static PetscReal SedimentBump(PetscReal x, PetscReal center, PetscReal R, PetscReal height)
{
    if (R <= 0.0) return 0.0;
    PetscReal t = (x - center) / R;
    if (PetscAbsReal(t) >= 1.0) return 0.0;
    return height * PetscExpReal(1.0 - 1.0 / (1.0 - t * t));
}

/* Sum of SedimentBump() humps along the bottom edge. If -sed_grain_x/-R were
 * not given (n_sed_grains == 0), falls back to the single-bump -geom_bump_R
 * behavior (centered at Lx/2) for backward compatibility with
 * build_geometry_sediment_grain.py / two_ice_grains_boundary. With
 * -sed_grain_x/-sed_grain_R/-sed_grain_h set, must match
 * build_geometry_multi_grain.py's SEDIMENT_GRAINS list. */
static PetscReal SedimentBumpField(const AppCtx *user, PetscReal x)
{
    if (user->n_sed_grains <= 0)
        return SedimentBump(x, 0.5 * user->Lx, user->geom_bump_R, user->geom_bump_R);

    PetscReal y = 0.0;
    for (PetscInt k = 0; k < user->n_sed_grains; k++)
        y += SedimentBump(x, user->sed_grain_x[k], user->sed_grain_R[k], user->sed_grain_h[k]);
    return y;
}

/* d/dx of SedimentBump(): g'(x) = g(x) * (-2t)/(R*(1-t^2)^2), t=(x-center)/R.
 * Vanishes at |t|->1 along with g() itself (C-infinity, compact support). */
static PetscReal SedimentBumpDeriv(PetscReal x, PetscReal center, PetscReal R, PetscReal height)
{
    if (R <= 0.0) return 0.0;
    PetscReal t = (x - center) / R;
    if (PetscAbsReal(t) >= 1.0) return 0.0;
    PetscReal g = SedimentBump(x, center, R, height);
    return g * (-2.0 * t) / (R * SQ(1.0 - t * t));
}

/* d/dx of SedimentBumpField() -- local slope of the actual floor curve,
 * used by the ice-shell distance-to-surface calculation below. */
static PetscReal SedimentBumpFieldDeriv(const AppCtx *user, PetscReal x)
{
    if (user->n_sed_grains <= 0)
        return SedimentBumpDeriv(x, 0.5 * user->Lx, user->geom_bump_R, user->geom_bump_R);

    PetscReal dy = 0.0;
    for (PetscInt k = 0; k < user->n_sed_grains; k++)
        dy += SedimentBumpDeriv(x, user->sed_grain_x[k], user->sed_grain_R[k], user->sed_grain_h[k]);
    return dy;
}

/* Sum of ceiling bumps (-top_grain_x/-R/-h) pushing DOWN from Ly.
 * Returns total downward displacement; caller computes y_top = Ly - TopBumpField().
 * Must match build_geometry_multi_grain.py's TOP_GRAINS list. */
static PetscReal TopBumpField(const AppCtx *user, PetscReal x)
{
    PetscReal h = 0.0;
    for (PetscInt k = 0; k < user->n_top_grains; k++)
        h += SedimentBump(x, user->top_grain_x[k], user->top_grain_R[k], user->top_grain_h[k]);
    return h;
}
PetscErrorCode FormInitialTwoIceGrainsBoundary2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx  = user->Lx;
    const PetscReal Ly  = user->Ly;
    const PetscReal eps       = user->eps;
    const PetscReal R0        = user->RCice0;   /* left grain (x=0), smaller  */
    const PetscReal R1        = user->RCice1;   /* right grain (x=Lx), larger */
    /* Equilibrium logistic profile: phi = 1/(1+exp(-(R-d)/eps))
     * = 0.5 - 0.5*tanh(0.5*(d-R)/eps), so tc = 0.5/eps. Initializing at the
     * model's own equilibrium width (1%-99% band = 9.2*eps) removes the
     * early width-relaxation transient. The old tc = 1/(sqrt(2)*0.75*eps)
     * was 1.89x steeper — a leftover from the removed eps_model=0.75*eps
     * residual scaling — and made every run start with the IC ~7 cells wide
     * relaxing to the equilibrium 13 cells over the first ~60 steps. */
    const PetscReal tc        = 0.5 / eps;

    const PetscReal c0x = 0.0,  c0y = 0.5 * Ly;
    const PetscReal c1x = Lx,   c1y = 0.5 * Ly;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D two ice grains, boundary-centered) ---\n"
        "  grain 0: centre = (%.4e, %.4e) m,  RCice0 = %.4e m\n"
        "  grain 1: centre = (%.4e, %.4e) m,  RCice1 = %.4e m\n",
        c0x, c0y, R0, c1x, c1y, R1);

    if (R0 <= 0.0 || R1 <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0 and RCice1 must be > 0 (got %.2e, %.2e)", R0, R1);
    if (R0 >= Lx || R1 >= Lx)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0/RCice1 must be < Lx (got %.2e, %.2e, Lx=%.2e)", R0, R1, Lx);
    if (R0 > 0.5 * Ly || R1 > 0.5 * Ly)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0/RCice1 must be <= Ly/2 (got %.2e, %.2e, Ly/2=%.2e)", R0, R1, 0.5*Ly);

    DM da;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da, U, &u); CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

    PetscInt per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            PetscReal x     = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
            PetscReal v     = (PetscReal)j / (PetscReal)(info.my + per);
            PetscReal y_bot = SedimentBumpField(user, x);
            PetscReal y_top = Ly - TopBumpField(user, x);
            PetscReal y     = y_bot + v * (y_top - y_bot);

            PetscReal d0 = PetscSqrtReal(SQ(x - c0x) + SQ(y - c0y));
            PetscReal d1 = PetscSqrtReal(SQ(x - c1x) + SQ(y - c1y));

            PetscReal ice0 = 0.5 - 0.5 * PetscTanhReal(tc * (d0 - R0));
            PetscReal ice1 = 0.5 - 0.5 * PetscTanhReal(tc * (d1 - R1));
            PetscReal ice  = PetscMin(PetscMax(ice0 + ice1, 0.0), 1.0);

            u[j][i].ice = ice;
            u[j][i].tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            PetscScalar rho_vs, temp_loc = u[j][i].tem;
            RhoVS_I(user, temp_loc, &rho_vs, NULL);
            { PetscReal _pa = PetscMax(0.0, 1.0 - ice);
              u[j][i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialMultiGrains2D
 *
 * N ice grains, each a tanh distance profile from a (cx,cy)/R given by
 * -ice_grain_cx/-ice_grain_cy/-ice_grain_R, summed and clamped to [0,1] --
 * the same construction as FormInitialTwoIceGrainsBoundary2D generalized
 * to an arbitrary number of grains.
 *
 * Physical coordinates use SedimentBumpField(), i.e. the bottom edge is the
 * sum of -sed_grain_x/-sed_grain_R bump humps (or the single -geom_bump_R
 * bump if those aren't set), matching build_geometry_multi_grain.py /
 * build_geometry_sediment_grain.py respectively.
 * =========================================================================*/
PetscErrorCode FormInitialMultiGrains2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx  = user->Lx;
    const PetscReal Ly  = user->Ly;
    const PetscReal eps       = user->eps;
    /* Equilibrium logistic profile: phi = 1/(1+exp(-(R-d)/eps))
     * = 0.5 - 0.5*tanh(0.5*(d-R)/eps), so tc = 0.5/eps. Initializing at the
     * model's own equilibrium width (1%-99% band = 9.2*eps) removes the
     * early width-relaxation transient. The old tc = 1/(sqrt(2)*0.75*eps)
     * was 1.89x steeper — a leftover from the removed eps_model=0.75*eps
     * residual scaling — and made every run start with the IC ~7 cells wide
     * relaxing to the equilibrium 13 cells over the first ~60 steps. */
    const PetscReal tc        = 0.5 / eps;
    const PetscInt  ng        = user->n_act;

    if (ng <= 0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "-ic_type multi_grains requires -ice_grain_cx/-ice_grain_cy/-ice_grain_R");

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D multi-grain) ---\n"
        "  %d ice grain(s), %d sediment bump(s)\n",
        (int)ng, (int)user->n_sed_grains);
    for (PetscInt k = 0; k < ng; k++) {
        PetscReal ax = user->ice_grain_ax[k];
        PetscReal ay = user->ice_grain_ay[k];
        PetscPrintf(PETSC_COMM_WORLD,
            "  ice grain %d: centre = (%.4e, %.4e) m,  ax = %.4e m,  ay = %.4e m\n",
            (int)k, user->cent[0][k], user->cent[1][k], ax, ay);
        if (ax <= 0.0 || ay <= 0.0)
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                    "-ice_grain_ax/ay[%d] must be > 0 (got ax=%.2e, ay=%.2e)", (int)k, ax, ay);
    }

    DM da;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da, U, &u); CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

    PetscReal *gx, *gy;
    PetscInt nbx, nby;
    ierr = GrevilleAbscissae(iga, 0, &gx, &nbx); CHKERRQ(ierr);
    ierr = GrevilleAbscissae(iga, 1, &gy, &nby); CHKERRQ(ierr);
    if (nbx != info.mx || nby != info.my)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_PLIB,
                "Greville abscissae count (%d,%d) != DOF grid size (%d,%d) "
                "-- multi_grains IC does not support -periodic 1",
                (int)nbx, (int)nby, (int)info.mx, (int)info.my);

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            PetscReal x     = Lx * gx[i];
            PetscReal v     = gy[j];
            PetscReal y_bot = SedimentBumpField(user, x);
            PetscReal y_top = Ly - TopBumpField(user, x);
            PetscReal y     = y_bot + v * (y_top - y_bot);

            PetscReal ice = 0.0;
            if (user->ic_grain_union) {
                /* Union form: sdf = min_k sdf_k is zero exactly on the sharp
                 * union surface, so phi = 0.5 lands there for ANY eps. Note
                 * (d-1)*sqrt(ax*ay) is the same per-grain signed distance the
                 * additive branch already uses -- the only change is min-vs-sum
                 * and a single tanh applied afterwards. See -ic_grain_union in
                 * NASA_types.h for why the additive form is eps-dependent at an
                 * overlapping neck. */
                PetscReal sdf = PETSC_MAX_REAL;
                for (PetscInt k = 0; k < ng; k++) {
                    PetscReal ax    = user->ice_grain_ax[k];
                    PetscReal ay    = user->ice_grain_ay[k];
                    PetscReal dx    = x - user->cent[0][k];
                    PetscReal dy    = y - user->cent[1][k];
                    PetscReal d     = PetscSqrtReal(SQ(dx / ax) + SQ(dy / ay));
                    PetscReal sdf_k = (d - 1.0) * PetscSqrtReal(ax * ay);
                    sdf = PetscMin(sdf, sdf_k);
                }
                ice = 0.5 - 0.5 * PetscTanhReal(tc * sdf);
            } else {
                for (PetscInt k = 0; k < ng; k++) {
                    PetscReal ax   = user->ice_grain_ax[k];
                    PetscReal ay   = user->ice_grain_ay[k];
                    PetscReal dx   = x - user->cent[0][k];
                    PetscReal dy   = y - user->cent[1][k];
                    PetscReal d    = PetscSqrtReal(SQ(dx / ax) + SQ(dy / ay)); /* =1 on ellipse boundary */
                    PetscReal tc_k = tc * PetscSqrtReal(ax * ay);              /* keeps interface width ~eps */
                    ice += 0.5 - 0.5 * PetscTanhReal(tc_k * (d - 1.0));
                }
            }
            for (PetscInt k = 0; k < user->n_ice_shells; k++) {
                PetscReal xs = user->ice_shell_x[k];
                PetscReal Rs = user->ice_shell_R[k];
                PetscReal ts = user->ice_shell_thickness[k];
                PetscReal dist;
                if (x < xs - Rs) {
                    /* left of the shell's segment: distance to its fixed
                     * endpoint (the floor curve is exactly 0 there) -- gives
                     * a naturally rounded cap, not a sharp/independent window */
                    dist = PetscSqrtReal(SQ(x - (xs - Rs)) + SQ(y));
                } else if (x > xs + Rs) {
                    dist = PetscSqrtReal(SQ(x - (xs + Rs)) + SQ(y));
                } else {
                    /* inside the segment: perpendicular distance to the floor
                     * curve's local tangent line (good approximation for a
                     * gently-curving bump); matches the endpoint formula
                     * continuously at x=xs+-Rs since slope->0 there too */
                    PetscReal slope = SedimentBumpFieldDeriv(user, x);
                    dist = (y - y_bot) / PetscSqrtReal(1.0 + SQ(slope));
                }
                PetscReal dn      = dist / ts;
                PetscReal tc_shell = tc * ts;
                ice += 0.5 - 0.5 * PetscTanhReal(tc_shell * (dn - 1.0));
            }
            for (PetscInt k = 0; k < user->n_ice_flats; k++) {
                PetscReal xf   = user->ice_flat_x[k];
                PetscReal Rf   = user->ice_flat_R[k];
                PetscReal Hf   = user->ice_flat_height[k];
                PetscReal dlat = PetscAbsReal(x - xf) / Rf;             /* lateral window: =1 at edge */
                PetscReal tc_lat = tc * Rf;
                PetscReal w    = 0.5 - 0.5 * PetscTanhReal(tc_lat * (dlat - 1.0));
                PetscReal flat = 0.5 - 0.5 * PetscTanhReal(tc * (y - Hf)); /* flat threshold, same
                                                                             * interface width as
                                                                             * everywhere else (no
                                                                             * R-dependent sharpening) */
                ice += w * flat;
            }
            ice = PetscMin(PetscMax(ice, 0.0), 1.0);

            u[j][i].ice = ice;
            u[j][i].tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            PetscScalar rho_vs, temp_loc = u[j][i].tem;
            RhoVS_I(user, temp_loc, &rho_vs, NULL);
            { PetscReal _pa = PetscMax(0.0, 1.0 - ice);
              u[j][i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    ierr = PetscFree(gx); CHKERRQ(ierr);
    ierr = PetscFree(gy); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialSedSlabGrain2D  —  simplest 3-phase (dof=4) initial condition.
 *
 * A flat SEDIMENT slab filling y < sed_slab_height (phi_s), plus ONE ice grain
 * (radius RCice) centred in x and sitting just above the slab. Air fills the
 * rest. This is the "simplest geometry first" test for the explicit sediment
 * phase (Effort 2).
 *
 * VALIDATION MODE: sed_slab_height <= 0 => phi_s = 0 everywhere and the ice
 * grain is centred at (Lx/2, Ly/2) — identical to FormInitialSingleIceGrain2D.
 * With phi_s = 0, Residual_A2 reduces to Residual_A1, so a dof=4 run must
 * reproduce the trusted 2-phase single-grain result.
 * =========================================================================*/
PetscErrorCode FormInitialSedSlabGrain2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx    = user->Lx;
    const PetscReal Ly    = user->Ly;
    const PetscReal eps   = user->eps;
    const PetscReal RCice = user->RCice;
    const PetscReal h_sed = user->sed_slab_height;
    const PetscReal tc    = 0.5 / eps;               /* equilibrium logistic width */
    const PetscReal cx    = 0.5 * Lx;
    /* Grain centre: validation (no slab) -> domain centre; otherwise seat the
     * grain just above the slab with a thin air gap so phi_a >= 0 at t=0. */
    const PetscReal cy    = (h_sed > 0.0) ? (h_sed + RCice + 4.0 * eps) : 0.5 * Ly;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D sediment slab + ice grain, 3-phase) ---\n"
        "  sed_slab_height = %.4e m %s\n"
        "  ice grain centre = (%.4e, %.4e) m,  RCice = %.4e m\n",
        h_sed, (h_sed > 0.0) ? "" : "(<=0: phi_s=0, 2-phase validation mode)",
        cx, cy, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice must be > 0 (got %.2e)", RCice);
    if (cy + RCice >= Ly)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "ice grain (cy=%.2e + R=%.2e) exceeds Ly=%.2e — lower "
                "sed_slab_height or RCice", cy, RCice, Ly);

    DM da;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    FieldSed **u;
    ierr = DMDAVecGetArray(da, U, &u); CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

    PetscInt per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
            PetscReal y = Ly * (PetscReal)j / (PetscReal)(info.my + per);

            /* Sediment slab (flat top at y = h_sed). Interface sharpened by
             * sed_width_factor to match the evolved ice width (width-mismatch fix). */
            PetscReal sed = (h_sed > 0.0)
                          ? 0.5 - 0.5 * PetscTanhReal(tc * user->sed_width_factor * (y - h_sed))
                          : 0.0;
            sed = PetscMin(PetscMax(sed, 0.0), 1.0);

            /* Ice grain. */
            PetscReal dist = PetscSqrtReal(SQ(x - cx) + SQ(y - cy));
            PetscReal ice  = 0.5 - 0.5 * PetscTanhReal(tc * (dist - RCice));
            ice = PetscMin(PetscMax(ice, 0.0), 1.0);
            /* Guarantee phi_a = 1 - ice - sed >= 0 at t=0. */
            if (ice > 1.0 - sed) ice = 1.0 - sed;

            u[j][i].ice = ice;
            u[j][i].sed = sed;
            u[j][i].tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            PetscScalar rho_vs, temp_loc = u[j][i].tem;
            RhoVS_I(user, temp_loc, &rho_vs, NULL);
            { PetscReal _pa = PetscMax(0.0, 1.0 - ice - sed);
              u[j][i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialSedIce1D  —  1D 3-phase (dof=4) initial condition.
 *
 * A three-layer stack along x:  SEDIMENT [0, x_sed] | ICE [x_sed, x_ice] | AIR.
 * x_sed = sed_slab_height, ice block width = 2*RCice. Cleanly isolates the
 * ice-sediment and ice-air interfaces in 1D — the sharpest test of the triple
 * well and whether spurious air appears at the ice-sediment contact.
 *
 * VALIDATION MODE (sed_slab_height <= 0): phi_s = 0 and a centred ice grain,
 * identical to FormInitialSingleIceGrain1D, so a dof=4 run reproduces the
 * 2-phase 1D result.
 * =========================================================================*/
PetscErrorCode FormInitialSedIce1D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx    = user->Lx;
    const PetscReal eps   = user->eps;
    const PetscReal RCice = user->RCice;
    const PetscReal h_sed = user->sed_slab_height;
    const PetscReal tc    = 0.5 / eps;
    const PetscReal x_sed = h_sed;                 /* sediment fills [0, x_sed] */
    const PetscReal x_ice = h_sed + 2.0 * RCice;   /* ice block [x_sed, x_ice] */
    const PetscReal cx    = 0.5 * Lx;              /* validation-mode grain centre */

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (1D sediment|ice|air stack, 3-phase) ---\n"
        "  sed_slab_height = %.4e m %s\n"
        "  ice block = [%.4e, %.4e] m,  RCice = %.4e m\n",
        h_sed, (h_sed > 0.0) ? "" : "(<=0: phi_s=0, 2-phase validation mode)",
        x_sed, x_ice, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "RCice must be > 0");
    if (h_sed > 0.0 && x_ice >= Lx)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "ice block end %.2e exceeds Lx=%.2e — lower sed_slab_height/RCice or raise Lx",
                x_ice, Lx);

    DM da;
    FieldSed *u;
    DMDALocalInfo info;
    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);            CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);           CHKERRQ(ierr);

    PetscInt per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + per);

        PetscReal sed, ice;
        if (h_sed > 0.0) {
            sed = 0.5 - 0.5 * PetscTanhReal(tc * (x - x_sed));
            /* smooth box = 1 between x_sed and x_ice, 0 outside */
            ice = (0.5 + 0.5 * PetscTanhReal(tc * (x - x_sed)))
                * (0.5 - 0.5 * PetscTanhReal(tc * (x - x_ice)));
        } else {
            sed = 0.0;
            ice = 0.5 - 0.5 * PetscTanhReal(tc * (PetscAbsReal(x - cx) - RCice));
        }
        sed = PetscMin(PetscMax(sed, 0.0), 1.0);
        ice = PetscMin(PetscMax(ice, 0.0), 1.0);
        if (ice > 1.0 - sed) ice = 1.0 - sed;   /* guarantee phi_a >= 0 */

        u[i].ice = ice;
        u[i].sed = sed;
        u[i].tem = user->temp0 + user->grad_temp0[0] * (x - 0.5 * Lx);

        PetscScalar rho_vs, temp_loc = u[i].tem;
        RhoVS_I(user, temp_loc, &rho_vs, NULL);
        { PetscReal _pa = PetscMax(0.0, 1.0 - ice - sed);
          u[i].rhov = rho_vs * (user->hum0 * _pa + (1.0 - _pa)); }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
