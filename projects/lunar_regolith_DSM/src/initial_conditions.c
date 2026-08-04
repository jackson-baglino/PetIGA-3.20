#include "initial_conditions.h"
#include "material_properties.h"


/* =========================================================================
 * GrevilleAbscissae
 *
 * Fills `g[0..count-1]` with the Greville abscissae (parametric DOF locations,
 * normalized to the physical parameter range) of the given IGA axis:
 * g_i = mean(U[i+1..i+p]) for the axis's knot vector U and degree p. For an
 * open-uniform knot vector this reduces to g_i = i/N (degree 1) but is
 * non-uniform for degree >= 2 -- this is the degree-agnostic replacement for
 * the `i/(mx+per)` index mapping used by IC functions on -geom_file domains.
 *
 * `count` is how many DOFs the caller actually has along this axis, i.e.
 * DMDALocalInfo's mx/my. Passing it in (rather than deriving m-p here) is what
 * makes this work under -periodic 1: IGAAxisInitUniform sets
 * nnp = periodic ? n-C : n+1 (petigaaxis.c:452), so a periodic axis has p
 * FEWER DOFs than basis-function slots, and the first `count` Greville values
 * are exactly the periodic DOF positions.
 *
 * Normalization uses U[p]..U[m-p], NOT U[0]..U[m]. They are identical for an
 * open knot vector, but IGAAxisInitUniform's periodic branch (petigaaxis.c:
 * 442-446) rewrites the p leading and trailing knots to lie OUTSIDE the
 * domain -- for N=8, p=2, C=1 on [0,1] it gives U[0]=-0.25, U[m]=1.25.
 * Dividing by that inflated span compressed the whole DOF grid: spacing came
 * out 0.0833 instead of the correct 0.125 = 1/N. U[p]..U[m-p] is the true
 * domain in both cases.
 *
 * Under periodicity the first abscissa is slightly negative (-h/2 for p=2),
 * which is correct: that DOF's basis function wraps around the domain.
 * Consumers must handle it, e.g. by minimum-image wrapping.
 *
 * Caller must PetscFree(*g).
 * =========================================================================*/
static PetscErrorCode GrevilleAbscissae(IGA iga, PetscInt dir, PetscInt count, PetscReal **g)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    IGAAxis axis;
    PetscInt p, m;
    PetscReal *U;
    ierr = IGAGetAxis(iga, dir, &axis);     CHKERRQ(ierr);
    ierr = IGAAxisGetDegree(axis, &p);      CHKERRQ(ierr);
    ierr = IGAAxisGetKnots(axis, &m, &U);   CHKERRQ(ierr);

    if (count < 1 || count > m - p)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "GrevilleAbscissae(dir=%d): asked for %d abscissae but the axis "
                "has only %d basis functions", (int)dir, (int)count, (int)(m - p));

    /* Normalize to [0,1]: a -geom_file axis's knots already span [0,1] (the
     * igakit builder's convention), but the default (-Nx/-Ny, no -geom_file)
     * axis built by IGAAxisInitUniform() spans [0,Lx]/[0,Ly] directly -- callers
     * (FormInitialMultiGrains2D) multiply the returned abscissae by Lx/Ly
     * themselves, so this must always return [0,1] regardless of which path
     * built the axis, or physical coordinates get scaled by Lx/Ly twice. */
    PetscReal Ui = U[p], Uf = U[m - p];
    PetscReal *greville;
    ierr = PetscMalloc1(count, &greville); CHKERRQ(ierr);
    for (PetscInt i = 0; i < count; i++) {
        PetscReal sum = 0.0;
        for (PetscInt k = 1; k <= p; k++) sum += U[i + k];
        greville[i] = (sum / (PetscReal)p - Ui) / (Uf - Ui);
    }

    *g = greville;
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

/* =========================================================================
 * WallBottom / WallTop -- the actual wall LINES, in physical y.
 *
 * The bump fields above only ever produced deviations from a flat floor at
 * y=0 and a flat ceiling at y=Ly. A tapered (wedge) domain needs a LINEAR
 * ramp, which compact-support bumps cannot represent at all, so each wall
 * gains an affine baseline:
 *
 *     y_bot(x) = wall_bot_y0 + wall_bot_slope*x + SedimentBumpField(x)
 *     y_top(x) = wall_top_y0 + wall_top_slope*x - TopBumpField(x)
 *
 * These MUST match build_geometry_multi_grain.py's build_surface(), which
 * cuts the mesh from the same two curves -- the IC is seeded by mapping
 * (u,v) through them, so any disagreement puts the ice somewhere the mesh
 * is not.
 *
 * Defaults (0, 0, Ly, 0) reproduce the previous flat-baseline behaviour
 * exactly, so every pre-existing geometry is untouched.
 *
 * Note these return the wall lines directly, NOT TopBumpField()'s "downward
 * displacement from Ly" convention. A rising ceiling threaded through that
 * sign convention is a sign error waiting to happen; callers should use
 * these and never hand-roll `Ly - TopBumpField(...)` again.
 * =========================================================================*/
static PetscReal WallBottom(const AppCtx *user, PetscReal x)
{
    return user->wall_bot_y0 + user->wall_bot_slope * x + SedimentBumpField(user, x);
}

static PetscReal WallTop(const AppCtx *user, PetscReal x)
{
    return user->wall_top_y0 + user->wall_top_slope * x - TopBumpField(user, x);
}

/* d/dx of WallBottom() -- local slope of the floor curve including the
 * affine baseline, used by the ice-shell distance-to-surface calculation. */
static PetscReal WallBottomDeriv(const AppCtx *user, PetscReal x)
{
    return user->wall_bot_slope + SedimentBumpFieldDeriv(user, x);
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
            PetscReal y_bot = WallBottom(user, x);
            PetscReal y_top = WallTop(user, x);
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

    /* The wedge band (below) is a standalone ice body, so a run may legitimately
     * have zero ice_grain_* circles. Only demand grains when nothing else would
     * put ice in the domain. */
    const PetscBool has_band = (PetscBool)(user->wedge_band_r2 > user->wedge_band_r1);

    if (ng <= 0 && !has_band)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "-ic_type multi_grains requires -ice_grain_cx/-ice_grain_cy/-ice_grain_R "
                "(or -wedge_band_r1/-wedge_band_r2 for a wedge-bridging band)");

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D multi-grain) ---\n"
        "  %d ice grain(s), %d sediment bump(s)\n",
        (int)ng, (int)user->n_sed_grains);
    if (has_band)
        PetscPrintf(PETSC_COMM_WORLD,
            "  wedge band: apex = (%.4e, %.4e) m,  r1 = %.4e m,  r2 = %.4e m\n",
            user->wedge_apex_x, user->wedge_apex_y,
            user->wedge_band_r1, user->wedge_band_r2);
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
    ierr = GrevilleAbscissae(iga, 0, info.mx, &gx); CHKERRQ(ierr);
    ierr = GrevilleAbscissae(iga, 1, info.my, &gy); CHKERRQ(ierr);

    /* Under -periodic 1 the grain field must wrap: a grain near x=0 has to be
     * seen by DOFs near x=Lx, or the microstructure is discontinuous across
     * the face that the solver (and the k_eff cell problem downstream) treats
     * as glued. Only meaningful with the union SDF -- the additive branch
     * sums per-grain tanh profiles and is unused for packings. */
    const PetscBool wrap = (PetscBool)(user->periodic == 1);

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            PetscReal x     = Lx * gx[i];
            PetscReal v     = gy[j];
            PetscReal y_bot = WallBottom(user, x);
            PetscReal y_top = WallTop(user, x);
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
                PetscReal sdf  = PETSC_MAX_REAL;
                PetscInt  kmin = -1;              /* grain achieving the min */
                for (PetscInt k = 0; k < ng; k++) {
                    PetscReal ax    = user->ice_grain_ax[k];
                    PetscReal ay    = user->ice_grain_ay[k];
                    PetscReal dx    = x - user->cent[0][k];
                    PetscReal dy    = y - user->cent[1][k];
                    if (wrap) {   /* minimum image: nearest periodic copy.
                                   * floor(t+0.5) == round(t); PETSc 3.20 has
                                   * no PetscRoundReal. */
                        dx -= Lx * PetscFloorReal(dx / Lx + 0.5);
                        dy -= Ly * PetscFloorReal(dy / Ly + 0.5);
                    }
                    PetscReal d     = PetscSqrtReal(SQ(dx / ax) + SQ(dy / ay));
                    PetscReal sdf_k = (d - 1.0) * PetscSqrtReal(ax * ay);
                    if (sdf_k < sdf) { sdf = sdf_k; kmin = k; }
                }

                /* min_k sdf_k is the distance to the nearest grain SURFACE, which
                 * equals the distance to the union BOUNDARY only when that nearest
                 * surface point is itself on the boundary. Where grains overlap it
                 * can be buried inside a neighbour, and then the formula reports a
                 * far smaller depth than the truth, leaving phi < 1 in bulk ice.
                 *
                 * Measured on the Molaro pair (R = 72.5/101 um, overlapped to form
                 * a neck): at the neck centre min_k sdf_k = -1.88 um against a true
                 * -16.41 um, giving phi = 0.87 / 0.97 at eps = 0.858 / 0.603 um
                 * instead of 1. The nearest point of grain 0's surface there lies
                 * 3.2 um inside grain 1.
                 *
                 * Trigger on OCCLUSION, not on how many grains contain the point:
                 * a point inside a single grain is affected too, whenever its
                 * nearest surface point is buried in a neighbour. Correction takes
                 * the nearest surface point of each of the two grains, discards it
                 * if buried, and falls back to the crease points where the circles
                 * cross -- those are always on the union boundary.
                 *
                 * Circles only (ax == ay); ellipses keep the approximation.
                 * Never fires for non-overlapping packings (production uses a 4 um
                 * contact gap), so it costs nothing there. */
                if (sdf < 0.0 && kmin >= 0) {
                    PetscReal ax0 = user->ice_grain_ax[kmin];
                    PetscReal ay0 = user->ice_grain_ay[kmin];
                    if (PetscAbsReal(ax0 - ay0) <= 1e-12 * ax0) {
                        PetscReal R0  = ax0;
                        PetscReal c0x = user->cent[0][kmin], c0y = user->cent[1][kmin];
                        PetscReal px = x, py = y;
                        if (wrap) {
                            PetscReal sx = px - c0x, sy = py - c0y;
                            px = c0x + sx - Lx * PetscFloorReal(sx / Lx + 0.5);
                            py = c0y + sy - Ly * PetscFloorReal(sy / Ly + 0.5);
                        }
                        PetscReal ddx = px - c0x, ddy = py - c0y;
                        PetscReal L   = PetscSqrtReal(ddx*ddx + ddy*ddy);
                        if (L > 0.0) {
                            PetscReal qx = c0x + R0*ddx/L, qy = c0y + R0*ddy/L;
                            /* is that surface point buried inside another grain? */
                            PetscInt  occ = -1;
                            for (PetscInt j = 0; j < ng && occ < 0; j++) {
                                if (j == kmin) continue;
                                PetscReal axj = user->ice_grain_ax[j];
                                PetscReal ayj = user->ice_grain_ay[j];
                                if (PetscAbsReal(axj - ayj) > 1e-12 * axj) continue;
                                PetscReal cjx = user->cent[0][j], cjy = user->cent[1][j];
                                if (wrap) {
                                    PetscReal sx = cjx - c0x, sy = cjy - c0y;
                                    cjx = c0x + sx - Lx * PetscFloorReal(sx / Lx + 0.5);
                                    cjy = c0y + sy - Ly * PetscFloorReal(sy / Ly + 0.5);
                                }
                                if (SQ(qx-cjx) + SQ(qy-cjy) < SQ(axj)) occ = j;
                            }
                            if (occ >= 0) {
                                PetscReal R1  = user->ice_grain_ax[occ];
                                PetscReal c1x = user->cent[0][occ], c1y = user->cent[1][occ];
                                if (wrap) {
                                    PetscReal sx = c1x - c0x, sy = c1y - c0y;
                                    c1x = c0x + sx - Lx * PetscFloorReal(sx / Lx + 0.5);
                                    c1y = c0y + sy - Ly * PetscFloorReal(sy / Ly + 0.5);
                                }
                                PetscReal best = PETSC_MAX_REAL;
                                /* the occluder's own nearest surface point, if visible */
                                PetscReal e1x = px - c1x, e1y = py - c1y;
                                PetscReal L1  = PetscSqrtReal(e1x*e1x + e1y*e1y);
                                if (L1 > 0.0) {
                                    PetscReal q1x = c1x + R1*e1x/L1, q1y = c1y + R1*e1y/L1;
                                    if (SQ(q1x-c0x) + SQ(q1y-c0y) >= SQ(R0))
                                        best = PetscMin(best, PetscAbsReal(R1 - L1));
                                }
                                /* crease points -- always on the union boundary */
                                PetscReal Dx = c1x - c0x, Dy = c1y - c0y;
                                PetscReal Dd = PetscSqrtReal(Dx*Dx + Dy*Dy);
                                if (Dd > 0.0 && Dd < R0 + R1 && Dd > PetscAbsReal(R0 - R1)) {
                                    PetscReal aa = (Dd*Dd + R0*R0 - R1*R1) / (2.0*Dd);
                                    PetscReal h2 = R0*R0 - aa*aa;
                                    if (h2 > 0.0) {
                                        PetscReal h  = PetscSqrtReal(h2);
                                        PetscReal mx = c0x + aa*Dx/Dd, my = c0y + aa*Dy/Dd;
                                        PetscReal ex = -Dy/Dd, ey = Dx/Dd;
                                        for (PetscInt sg = -1; sg <= 1; sg += 2) {
                                            PetscReal cxp = mx + sg*h*ex, cyp = my + sg*h*ey;
                                            best = PetscMin(best,
                                                PetscSqrtReal(SQ(px-cxp) + SQ(py-cyp)));
                                        }
                                    }
                                }
                                if (best < PETSC_MAX_REAL) sdf = -best;
                            }
                        }
                    }
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
                    PetscReal slope = WallBottomDeriv(user, x);
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
            /* Wedge-bridging band: the annulus r1 <= |X - apex| <= r2.
             *
             * In a wedge the two walls are rays from the apex, and a circular
             * arc centred on that apex is perpendicular to every such ray --
             * so an apex-centred annulus is exactly the shape that meets BOTH
             * walls at 90 degrees, the model's natural (Neumann) contact angle.
             * Initialising anything else, a circle in particular, starts the
             * run with a contact-angle relaxation transient that swamps the
             * slow wedge-driven migration we are trying to measure.
             *
             * max(r1-rho, rho-r2) is the exact signed distance to the annulus
             * (negative inside), so the same tanh/tc used for every other
             * feature gives the equilibrium interface width here too. */
            if (has_band) {
                PetscReal dxb = x - user->wedge_apex_x;
                PetscReal dyb = y - user->wedge_apex_y;
                PetscReal rho = PetscSqrtReal(SQ(dxb) + SQ(dyb));
                PetscReal sdf = PetscMax(user->wedge_band_r1 - rho,
                                         rho - user->wedge_band_r2);
                ice += 0.5 - 0.5 * PetscTanhReal(tc * sdf);
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
 * ReadGrainsFromFile
 *
 * Populates user->cent[0..1][k] and user->radius[k] from the whitespace-
 * delimited grain list at user->grains_file, setting user->n_act. Written for
 * the packings produced by preprocess/generate_packing.py:
 *
 *     Lx Ly          <- optional 2-token header (domain size in metres)
 *     x  y  r        <- one row per grain, metres
 *
 * A 4-token row (x y z r) from the legacy 3D-flavoured generators is also
 * accepted; z is read and discarded because the model is 2D.
 *
 * Row width is fixed by the FIRST data row and every later row must match.
 * The legacy reader in dry_snow_metamorphism inferred the width from the
 * first line and then accepted any row with >= 3 fields, so a 3-field row in
 * a 4-field file left the radius holding whatever the previous iteration had
 * -- silently producing grains of the wrong size. Here a ragged file is an
 * error.
 *
 * Parsed on rank 0 and broadcast, so every rank sees identical geometry
 * regardless of filesystem visibility.
 * =========================================================================*/
PetscErrorCode ReadGrainsFromFile(AppCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       ng = 0;
    PetscFunctionBegin;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    if (!user->grains_file[0])
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "-ic_type multi_grains_file requires -grains_file <path>");

    if (rank == 0) {
        FILE *fp = fopen(user->grains_file, "r");
        if (!fp)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                    "Cannot open -grains_file '%s'", user->grains_file);

        char     line[512];
        PetscInt ncol = 0;          /* 3 (x y r) or 4 (x y z r); 0 until known */
        PetscInt lineno = 0;

        while (fgets(line, (int)sizeof(line), fp)) {
            lineno++;

            /* Skip blank lines and '#' comments */
            char *s = line;
            while (*s == ' ' || *s == '\t') s++;
            if (*s == '\0' || *s == '\n' || *s == '\r' || *s == '#') continue;

            PetscReal v[4];
            PetscInt  nf = (PetscInt)sscanf(s, "%lf %lf %lf %lf",
                                            &v[0], &v[1], &v[2], &v[3]);

            /* A leading 2-token line is the "Lx Ly" domain header. */
            if (nf == 2) {
                if (ng == 0 && ncol == 0) continue;
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                        "%s:%d: 2 fields, but a domain header is only allowed "
                        "before the first grain", user->grains_file, (int)lineno);
            }
            if (nf < 3)
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                        "%s:%d: expected 3 (x y r) or 4 (x y z r) fields, got %d",
                        user->grains_file, (int)lineno, (int)nf);

            if (ncol == 0) {
                ncol = nf;                      /* fixed by the first data row */
            } else if (nf != ncol) {
                /* A 3-field first row followed by 4-field rows is almost
                 * always a legacy 3D "Lx Ly Lz" header, which is ambiguous
                 * with a 3-field "x y r" grain row and so cannot be detected
                 * from the row alone. Name it rather than leaving the user
                 * with a bare ragged-file complaint. */
                if (ng == 1 && ncol == 3 && nf == 4)
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                            "%s:%d: looks like the legacy 3D format -- a "
                            "'Lx Ly Lz' header (indistinguishable from an "
                            "'x y r' grain row) followed by 'x y z r' rows. "
                            "The model is 2D: drop the Lz field from the "
                            "header, or drop the header entirely.",
                            user->grains_file, (int)lineno);
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                        "%s:%d: ragged file -- %d fields here but %d in the "
                        "first grain row", user->grains_file, (int)lineno,
                        (int)nf, (int)ncol);
            }

            if (ng >= user->n_grain_max) {
                fclose(fp);
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                        "%s holds more than -n_grain_max (%d) grains; raise it",
                        user->grains_file, (int)user->n_grain_max);
            }

            PetscReal gx = v[0], gy = v[1];
            PetscReal gr = (ncol == 4) ? v[3] : v[2];   /* z (v[2]) discarded */
            if (gr <= 0.0) {
                fclose(fp);
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                        "%s:%d: radius must be > 0 (got %.6e)",
                        user->grains_file, (int)lineno, (double)gr);
            }

            user->cent[0][ng] = gx;
            user->cent[1][ng] = gy;
            user->cent[2][ng] = 0.0;
            user->radius[ng]  = gr;
            ng++;
        }
        fclose(fp);

        if (ng == 0)
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                    "%s contains no grain rows", user->grains_file);
    }

    ierr = MPI_Bcast(&ng, 1, MPIU_INT, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(user->cent[0], (PetscMPIInt)ng, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(user->cent[1], (PetscMPIInt)ng, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(user->radius,  (PetscMPIInt)ng, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    user->n_act = ng;
    for (PetscInt k = 0; k < ng; k++) {
        user->ice_grain_ax[k] = user->radius[k];
        user->ice_grain_ay[k] = user->radius[k];
    }

    PetscPrintf(PETSC_COMM_WORLD,
        "  -grains_file: %d grains from %s\n", (int)ng, user->grains_file);
    PetscFunctionReturn(0);
}
