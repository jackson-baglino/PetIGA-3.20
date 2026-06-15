#include "initial_conditions.h"
#include "material_properties.h"


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
    const PetscReal Ly    = user->Ly;
    const PetscReal eps   = user->eps;
    const PetscReal RCice = user->RCice;
    const PetscReal tc    = 1.0 / (PetscSqrtReal(2.0) * eps);
    const PetscReal cx    = 0.5 * Lx;
    const PetscReal cy    = 0.5 * Ly;

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

    const PetscReal Lx    = user->Lx;
    const PetscReal eps   = user->eps;
    const PetscReal RCice = user->RCice;
    const PetscReal tc    = 1.0 / (PetscSqrtReal(2.0) * eps);
    const PetscReal cx    = 0.5 * Lx;

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
PetscErrorCode FormInitialTwoIceGrainsBoundary2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx  = user->Lx;
    const PetscReal Ly  = user->Ly;
    const PetscReal eps = user->eps;
    const PetscReal R0  = user->RCice0;   /* left grain (x=0), smaller  */
    const PetscReal R1  = user->RCice1;   /* right grain (x=Lx), larger */
    const PetscReal tc  = 1.0 / (PetscSqrtReal(2.0) * eps);

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
            PetscReal x = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
            PetscReal v = (PetscReal)j / (PetscReal)(info.my + per);
            PetscReal bump = SedimentBump(x, 0.5 * Lx, user->geom_bump_R, user->geom_bump_R);
            PetscReal y = bump + v * (Ly - bump);

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
