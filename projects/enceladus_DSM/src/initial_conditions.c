#include "initial_conditions.h"
#include "material_properties.h"


/* =========================================================================
 * HOW THIS FILE IS ORGANISED
 *
 * Every initial condition does the same six things: create the node DM, walk
 * the local DOF box, turn (i,j) into physical (x,y), evaluate a phase field,
 * clamp it to [0,1], and write ice/temperature/vapour. Only the fourth step
 * differs between ICs.
 *
 * So that scaffolding lives once, in FillIC1D() / FillIC2D(), and each IC
 * supplies just a *shape function* -- a pure PetscReal(x[,y],user) that
 * returns the unclamped ice field. The public FormInitial*() entry points are
 * thin wrappers: validate arguments, print the banner, hand a shape to a
 * filler. Adding an IC means writing one shape function.
 *
 * Layout below:
 *   1. GrevilleAbscissae      -- parametric DOF locations
 *   2. wall / bump geometry   -- the curved domain boundaries
 *   3. FillIC1D / FillIC2D    -- the shared scaffolding
 *   4. shape functions        -- one per IC, plus the multi-grain feature set
 *   5. FormInitial*()         -- public entry points
 *   6. ReadGrainsFromFile     -- packing loader
 * =========================================================================*/


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
     * multiply the returned abscissae by Lx/Ly themselves, so this must always
     * return [0,1] regardless of which path built the axis, or physical
     * coordinates get scaled by Lx/Ly twice. */
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


/* =========================================================================
 * 2. WALL / BUMP GEOMETRY
 * =========================================================================*/

/* C-infinity bump g(x) = height*exp(1 - 1/(1-t^2)) for |t|<1, t=(x-center)/R;
 * 0 outside -- must match build_geometry_sediment_grain.py's _bump(). */
static PetscReal SedimentBump(PetscReal x, PetscReal center, PetscReal R, PetscReal height)
{
    if (R <= 0.0 || height == 0.0) return 0.0;
    PetscReal t = (x - center) / R;
    if (PetscAbsReal(t) >= 1.0) return 0.0;
    return height * PetscExpReal(1.0 - 1.0 / (1.0 - t * t));
}

/* Total floor displacement at x: the sum of every sediment bump. */
static PetscReal SedimentBumpField(const AppCtx *user, PetscReal x)
{
    PetscReal h = 0.0;
    for (PetscInt k = 0; k < user->n_sed_grains; k++)
        h += SedimentBump(x, user->sed_grain_x[k], user->sed_grain_R[k], user->sed_grain_h[k]);
    return h;
}

/* d/dx of SedimentBump(). */
static PetscReal SedimentBumpDeriv(PetscReal x, PetscReal center, PetscReal R, PetscReal height)
{
    if (R <= 0.0 || height == 0.0) return 0.0;
    PetscReal t = (x - center) / R;
    if (PetscAbsReal(t) >= 1.0) return 0.0;
    PetscReal s = 1.0 - t * t;
    return height * PetscExpReal(1.0 - 1.0 / s) * (-2.0 * t / (s * s)) / R;
}

static PetscReal SedimentBumpFieldDeriv(const AppCtx *user, PetscReal x)
{
    PetscReal d = 0.0;
    for (PetscInt k = 0; k < user->n_sed_grains; k++)
        d += SedimentBumpDeriv(x, user->sed_grain_x[k], user->sed_grain_R[k], user->sed_grain_h[k]);
    return d;
}

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
 * The floor is y=0 raised by the sediment bumps; the ceiling is y=Ly lowered
 * by the top bumps:
 *
 *     y_bot(x) = SedimentBumpField(x)
 *     y_top(x) = Ly - TopBumpField(x)
 *
 * These MUST match build_geometry_multi_grain.py's build_surface(), which
 * cuts the mesh from the same two curves -- the IC is seeded by mapping
 * (u,v) through them, so any disagreement puts the ice somewhere the mesh
 * is not.
 *
 * With no bumps configured these are y=0 and y=Ly, so IC_COORD_WALLS
 * degenerates to the plain rectangular mapping.
 *
 * Note these return the wall lines directly, NOT TopBumpField()'s "downward
 * displacement from Ly" convention. A rising ceiling threaded through that
 * sign convention is a sign error waiting to happen; callers should use
 * these and never hand-roll `Ly - TopBumpField(...)` again.
 * =========================================================================*/
static PetscReal WallBottom(const AppCtx *user, PetscReal x)
{
    return SedimentBumpField(user, x);
}

static PetscReal WallTop(const AppCtx *user, PetscReal x)
{
    return user->Ly - TopBumpField(user, x);
}

/* d/dx of WallBottom() -- local slope of the floor curve including the
 * affine baseline, used by the ice-shell distance-to-surface calculation. */
static PetscReal WallBottomDeriv(const AppCtx *user, PetscReal x)
{
    return SedimentBumpFieldDeriv(user, x);
}


/* =========================================================================
 * 3. SHARED SCAFFOLDING
 * =========================================================================*/

/* Interface sharpness used by every shape function below.
 *
 * The equilibrium logistic profile is phi = 1/(1+exp(-(R-d)/eps))
 * = 0.5 - 0.5*tanh(0.5*(d-R)/eps), so the tanh coefficient is 0.5/eps.
 * Initializing at the model's own equilibrium width (1%-99% band = 9.2*eps)
 * removes the early width-relaxation transient. The old tc = 1/(sqrt(2)*
 * 0.75*eps) was 1.89x steeper -- a leftover from the removed eps_model =
 * 0.75*eps residual scaling -- and made every run start with the IC ~7 cells
 * wide relaxing to the equilibrium 13 cells over the first ~60 steps. */
static inline PetscReal PhaseSlope(const AppCtx *user) { return 0.5 / user->eps; }

/* A shape function returns the UNCLAMPED ice field at a physical point.
 * Everything it needs comes from `user`; the fillers clamp to [0,1]. */
typedef PetscReal (*IceShape1D)(PetscReal x, const AppCtx *user);
typedef PetscReal (*IceShape2D)(PetscReal x, PetscReal y, const AppCtx *user);

/* How (i,j) becomes physical (x,y). */
typedef enum {
    IC_COORD_UNIFORM,   /* x = Lx*i/(mx+per),      y = Ly*j/(my+per)          */
    IC_COORD_WALLS,     /* x as above,             y ruled between the walls  */
    IC_COORD_GREVILLE   /* x = Lx*greville_x[i],   y ruled between the walls  */
} ICCoordMode;

/* The closure every IC shares: temperature is the linear background field and
 * vapour is saturated inside ice, hum0-scaled in air, blended by ice fraction. */
static void SetNodeFields(Field *f, PetscReal ice, PetscReal tem, AppCtx *user)
{
    PetscScalar rho_vs;
    RhoVS_I(user, tem, &rho_vs, NULL);

    PetscReal air = PetscMax(0.0, 1.0 - ice);

    f->ice  = ice;
    f->tem  = tem;
    f->rhov = rho_vs * (user->hum0 * air + (1.0 - air));
}

/* Walk the local 1D DOF box, evaluate `shape`, write the node fields. */
static PetscErrorCode FillIC1D(IGA iga, Vec U, AppCtx *user, IceShape1D shape)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    DM            da;
    Field        *u;   /* 1D: plain pointer, not pointer-to-pointer */
    DMDALocalInfo info;

    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);           CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);          CHKERRQ(ierr);

    const PetscReal Lx  = user->Lx;
    const PetscInt  per = (user->periodic == 1) ? user->p - 1 : -1;

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        PetscReal x   = Lx * (PetscReal)i / (PetscReal)(info.mx + per);
        PetscReal ice = PetscMin(PetscMax(shape(x, user), 0.0), 1.0);
        PetscReal tem = user->temp0 + user->grad_temp0[0] * (x - 0.5 * Lx);

        SetNodeFields(&u[i], ice, tem, user);
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

/* Walk the local 2D DOF box, evaluate `shape`, write the node fields.
 *
 * `mode` picks the index->physical map. IC_COORD_WALLS and IC_COORD_GREVILLE
 * rule y between WallBottom(x) and WallTop(x), which is what makes the IC
 * follow a curved or tapered -geom_file domain instead of a bounding box. */
static PetscErrorCode FillIC2D(IGA iga, Vec U, AppCtx *user,
                               ICCoordMode mode, IceShape2D shape)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    DM            da;
    Field       **u;
    DMDALocalInfo info;

    ierr = IGACreateNodeDM(iga, user->dof, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);           CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info);          CHKERRQ(ierr);

    const PetscReal Lx  = user->Lx;
    const PetscReal Ly  = user->Ly;
    const PetscInt  per = (user->periodic == 1) ? user->p - 1 : -1;

    PetscReal *gx = NULL, *gy = NULL;
    if (mode == IC_COORD_GREVILLE) {
        ierr = GrevilleAbscissae(iga, 0, info.mx, &gx); CHKERRQ(ierr);
        ierr = GrevilleAbscissae(iga, 1, info.my, &gy); CHKERRQ(ierr);
    }

    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {

            PetscReal x = (mode == IC_COORD_GREVILLE)
                        ? Lx * gx[i]
                        : Lx * (PetscReal)i / (PetscReal)(info.mx + per);

            PetscReal y;
            if (mode == IC_COORD_UNIFORM) {
                y = Ly * (PetscReal)j / (PetscReal)(info.my + per);
            } else {
                PetscReal v = (mode == IC_COORD_GREVILLE)
                            ? gy[j]
                            : (PetscReal)j / (PetscReal)(info.my + per);
                PetscReal y_bot = WallBottom(user, x);
                PetscReal y_top = WallTop(user, x);
                y = y_bot + v * (y_top - y_bot);
            }

            PetscReal ice = PetscMin(PetscMax(shape(x, y, user), 0.0), 1.0);
            PetscReal tem = user->temp0
                          + user->grad_temp0[0] * (x - 0.5 * Lx)
                          + user->grad_temp0[1] * (y - 0.5 * Ly);

            SetNodeFields(&u[j][i], ice, tem, user);
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDestroy(&da);                 CHKERRQ(ierr);
    if (mode == IC_COORD_GREVILLE) {
        ierr = PetscFree(gx); CHKERRQ(ierr);
        ierr = PetscFree(gy); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}


/* =========================================================================
 * 4. SHAPE FUNCTIONS
 * =========================================================================*/

/* Diffuse slab between x_lo and x_hi.
 *   flag_tIC == 2  ->  flat interface (ice fills [0, 0.5 Lx])
 *   otherwise      ->  centred slab   (ice fills [0.35 Lx, 0.65 Lx]) */
static PetscReal ShapeSlab1D(PetscReal x, const AppCtx *user)
{
    const PetscReal Lx  = user->Lx;
    const PetscReal eps = user->eps;
    const PetscReal x_lo = (user->flag_tIC == 2) ? 0.0       : 0.35 * Lx;
    const PetscReal x_hi = (user->flag_tIC == 2) ? 0.5 * Lx  : 0.65 * Lx;

    return 0.5 * (PetscTanhReal(0.5 * (x - x_lo) / eps)
                - PetscTanhReal(0.5 * (x - x_hi) / eps));
}

/* Layered ice slab on [0, Ly/2], air above -- the k_eff analytic benchmark.
 * See the FormInitialIceSlab2D banner comment for why this is built from a
 * signed distance rather than a single tanh. */
static PetscReal ShapeIceSlab2D(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal Ly = user->Ly;

    /* Signed distance to the nearest ice/air interface, negative inside the
     * ice. Periodic: interfaces at y = Ly/2 and y = 0 == Ly.
     * Non-periodic: the single interface at y = Ly/2. */
    PetscReal dist_ice;
    if (user->periodic == 1) {
        dist_ice = (y <= 0.5 * Ly)
                 ? -PetscMin(y, 0.5 * Ly - y)          /* in the ice */
                 :  PetscMin(y - 0.5 * Ly, Ly - y);    /* in the air */
    } else {
        dist_ice = y - 0.5 * Ly;
    }

    return 0.5 - 0.5 * PetscTanhReal(0.5 * dist_ice / user->eps);
}

/* Single ice circle of radius RCice at the domain centre. */
static PetscReal ShapeSingleIceGrain2D(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal dist = PetscSqrtReal(SQ(x - 0.5 * user->Lx)
                                       + SQ(y - 0.5 * user->Ly));

    return 0.5 - 0.5 * PetscTanhReal(PhaseSlope(user) * (dist - user->RCice));
}

/* 1D cross-section through the centre of a single ice grain. */
static PetscReal ShapeSingleIceGrain1D(PetscReal x, const AppCtx *user)
{
    const PetscReal dist = PetscAbsReal(x - 0.5 * user->Lx);

    return 0.5 - 0.5 * PetscTanhReal(PhaseSlope(user) * (dist - user->RCice));
}

/* Two ice semicircles centred ON the x=0 and x=Lx boundaries, summed. */
static PetscReal ShapeTwoIceGrainsBoundary2D(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal tc = PhaseSlope(user);
    const PetscReal cy = 0.5 * user->Ly;

    PetscReal d0 = PetscSqrtReal(SQ(x - 0.0)      + SQ(y - cy));
    PetscReal d1 = PetscSqrtReal(SQ(x - user->Lx) + SQ(y - cy));

    return (0.5 - 0.5 * PetscTanhReal(tc * (d0 - user->RCice0)))
         + (0.5 - 0.5 * PetscTanhReal(tc * (d1 - user->RCice1)));
}


/* ---- multi-grain feature fields -----------------------------------------
 *
 * The multi-grain IC is a SUM of independent ice bodies. Each helper below
 * contributes one family; ShapeMultiGrains2D adds them up. Keeping them
 * separate is what makes the geometry legible -- the alternative is one
 * 300-line loop where the grain, shell, flat, wedge and bridge logic are
 * interleaved and the nesting runs ten braces deep.
 * ------------------------------------------------------------------------*/

/* Minimum-image displacement under -periodic 1: the nearest periodic copy.
 * floor(t+0.5) == round(t); PETSc 3.20 has no PetscRoundReal. */
static inline PetscReal MinImage(PetscReal d, PetscReal L) {
    return d - L * PetscFloorReal(d / L + 0.5);
}

/* Additive form: sum of per-grain tanh profiles.
 * Each grain's own field crosses 0.5 at exactly its own surface for any eps,
 * which is what makes contact eps-independent for tangent (non-overlapping)
 * packings. See -ic_grain_union in enceladus_types.h for why this form is
 * eps-DEPENDENT at an overlapping neck. */
static PetscReal GrainFieldAdditive(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal tc  = PhaseSlope(user);
    PetscReal       ice = 0.0;

    for (PetscInt k = 0; k < user->n_act; k++) {
        PetscReal ax   = user->ice_grain_ax[k];
        PetscReal ay   = user->ice_grain_ay[k];
        PetscReal dx   = x - user->cent[0][k];
        PetscReal dy   = y - user->cent[1][k];
        PetscReal d    = PetscSqrtReal(SQ(dx / ax) + SQ(dy / ay)); /* =1 on ellipse boundary */
        PetscReal tc_k = tc * PetscSqrtReal(ax * ay);              /* keeps interface width ~eps */
        ice += 0.5 - 0.5 * PetscTanhReal(tc_k * (d - 1.0));
    }
    return ice;
}

/* Distance from (px,py) to the true union boundary of the overlapping pair
 * (kmin, occ), or PETSC_MAX_REAL if no correction applies.
 *
 * WHY THIS EXISTS. min_k sdf_k is the distance to the nearest grain SURFACE,
 * which equals the distance to the union BOUNDARY only when that nearest
 * surface point is itself on the boundary. Where grains overlap it can be
 * buried inside a neighbour, and then the formula reports a far smaller depth
 * than the truth, leaving phi < 1 in bulk ice.
 *
 * Measured on the Molaro pair (R = 72.5/101 um, overlapped to form a neck):
 * at the neck centre min_k sdf_k = -1.88 um against a true -16.41 um, giving
 * phi = 0.87 / 0.97 at eps = 0.858 / 0.603 um instead of 1. The nearest point
 * of grain 0's surface there lies 3.2 um inside grain 1.
 *
 * The fix takes the nearest surface point of each of the two grains, discards
 * it if buried, and falls back to the crease points where the circles cross --
 * those are always on the union boundary. */
static PetscReal UnionBoundaryDistance(PetscReal px, PetscReal py,
                                       PetscInt kmin, PetscInt occ,
                                       const AppCtx *user)
{
    const PetscBool wrap = (PetscBool)(user->periodic == 1);
    const PetscReal Lx = user->Lx, Ly = user->Ly;

    const PetscReal R0  = user->ice_grain_ax[kmin];
    const PetscReal c0x = user->cent[0][kmin], c0y = user->cent[1][kmin];

    PetscReal R1  = user->ice_grain_ax[occ];
    PetscReal c1x = user->cent[0][occ], c1y = user->cent[1][occ];
    if (wrap) {
        c1x = c0x + MinImage(c1x - c0x, Lx);
        c1y = c0y + MinImage(c1y - c0y, Ly);
    }

    PetscReal best = PETSC_MAX_REAL;

    /* the occluder's own nearest surface point, if visible */
    PetscReal e1x = px - c1x, e1y = py - c1y;
    PetscReal L1  = PetscSqrtReal(e1x * e1x + e1y * e1y);
    if (L1 > 0.0) {
        PetscReal q1x = c1x + R1 * e1x / L1, q1y = c1y + R1 * e1y / L1;
        if (SQ(q1x - c0x) + SQ(q1y - c0y) >= SQ(R0))
            best = PetscMin(best, PetscAbsReal(R1 - L1));
    }

    /* crease points -- always on the union boundary */
    PetscReal Dx = c1x - c0x, Dy = c1y - c0y;
    PetscReal Dd = PetscSqrtReal(Dx * Dx + Dy * Dy);
    if (Dd > 0.0 && Dd < R0 + R1 && Dd > PetscAbsReal(R0 - R1)) {
        PetscReal aa = (Dd * Dd + R0 * R0 - R1 * R1) / (2.0 * Dd);
        PetscReal h2 = R0 * R0 - aa * aa;
        if (h2 > 0.0) {
            PetscReal h  = PetscSqrtReal(h2);
            PetscReal mx = c0x + aa * Dx / Dd, my = c0y + aa * Dy / Dd;
            PetscReal ex = -Dy / Dd,           ey = Dx / Dd;
            for (PetscInt sg = -1; sg <= 1; sg += 2) {
                PetscReal cxp = mx + sg * h * ex, cyp = my + sg * h * ey;
                best = PetscMin(best, PetscSqrtReal(SQ(px - cxp) + SQ(py - cyp)));
            }
        }
    }
    return best;
}

/* Which OTHER circular grain buries the point q, or -1 if none does. */
static PetscInt OccludingGrain(PetscReal qx, PetscReal qy,
                               PetscReal c0x, PetscReal c0y, PetscInt kmin,
                               const AppCtx *user)
{
    const PetscBool wrap = (PetscBool)(user->periodic == 1);
    const PetscReal Lx = user->Lx, Ly = user->Ly;

    for (PetscInt jg = 0; jg < user->n_act; jg++) {
        if (jg == kmin) continue;
        PetscReal axj = user->ice_grain_ax[jg];
        PetscReal ayj = user->ice_grain_ay[jg];
        if (PetscAbsReal(axj - ayj) > 1e-12 * axj) continue;   /* circles only */

        PetscReal cjx = user->cent[0][jg], cjy = user->cent[1][jg];
        if (wrap) {
            cjx = c0x + MinImage(cjx - c0x, Lx);
            cjy = c0y + MinImage(cjy - c0y, Ly);
        }
        if (SQ(qx - cjx) + SQ(qy - cjy) < SQ(axj)) return jg;
    }
    return -1;
}

/* Union form: sdf = min_k sdf_k is zero exactly on the sharp union surface,
 * so phi = 0.5 lands there for ANY eps. (d-1)*sqrt(ax*ay) is the same
 * per-grain signed distance the additive branch uses -- the only change is
 * min-vs-sum and a single tanh applied afterwards.
 *
 * Under -periodic 1 the grain field must wrap: a grain near x=0 has to be seen
 * by DOFs near x=Lx, or the microstructure is discontinuous across the face
 * that the solver (and the k_eff cell problem downstream) treats as glued.
 *
 * The crease correction below is circles-only (ax == ay); ellipses keep the
 * approximation. It never fires for non-overlapping packings (production uses
 * a 4 um contact gap), so it costs nothing there. */
static PetscReal GrainFieldUnion(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscBool wrap = (PetscBool)(user->periodic == 1);
    const PetscReal Lx = user->Lx, Ly = user->Ly;

    PetscReal sdf  = PETSC_MAX_REAL;
    PetscInt  kmin = -1;                  /* grain achieving the min */

    for (PetscInt k = 0; k < user->n_act; k++) {
        PetscReal ax = user->ice_grain_ax[k];
        PetscReal ay = user->ice_grain_ay[k];
        PetscReal dx = x - user->cent[0][k];
        PetscReal dy = y - user->cent[1][k];
        if (wrap) {
            dx = MinImage(dx, Lx);
            dy = MinImage(dy, Ly);
        }
        PetscReal d     = PetscSqrtReal(SQ(dx / ax) + SQ(dy / ay));
        PetscReal sdf_k = (d - 1.0) * PetscSqrtReal(ax * ay);
        if (sdf_k < sdf) { sdf = sdf_k; kmin = k; }
    }

    /* Inside the union: check whether the nearest surface point is buried,
     * and if so fall back to the true union boundary. Trigger on OCCLUSION,
     * not on how many grains contain the point: a point inside a single grain
     * is affected too, whenever its nearest surface point is buried. */
    if (sdf < 0.0 && kmin >= 0) {
        PetscReal ax0 = user->ice_grain_ax[kmin];
        PetscReal ay0 = user->ice_grain_ay[kmin];

        if (PetscAbsReal(ax0 - ay0) <= 1e-12 * ax0) {       /* circles only */
            PetscReal R0  = ax0;
            PetscReal c0x = user->cent[0][kmin], c0y = user->cent[1][kmin];

            PetscReal px = x, py = y;
            if (wrap) {
                px = c0x + MinImage(px - c0x, Lx);
                py = c0y + MinImage(py - c0y, Ly);
            }
            PetscReal ddx = px - c0x, ddy = py - c0y;
            PetscReal L   = PetscSqrtReal(ddx * ddx + ddy * ddy);

            if (L > 0.0) {
                PetscReal qx = c0x + R0 * ddx / L, qy = c0y + R0 * ddy / L;
                PetscInt  occ = OccludingGrain(qx, qy, c0x, c0y, kmin, user);
                if (occ >= 0) {
                    PetscReal best = UnionBoundaryDistance(px, py, kmin, occ, user);
                    if (best < PETSC_MAX_REAL) sdf = -best;
                }
            }
        }
    }

    return 0.5 - 0.5 * PetscTanhReal(PhaseSlope(user) * sdf);
}

/* Conformal ice shells hugging the floor curve. */
static PetscReal IceShellField(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal tc    = PhaseSlope(user);
    const PetscReal y_bot = WallBottom(user, x);
    PetscReal       ice   = 0.0;

    for (PetscInt k = 0; k < user->n_ice_shells; k++) {
        PetscReal xs = user->ice_shell_x[k];
        PetscReal Rs = user->ice_shell_R[k];
        PetscReal ts = user->ice_shell_thickness[k];
        PetscReal dist;

        if (x < xs - Rs) {
            /* left of the shell's segment: distance to its fixed endpoint (the
             * floor curve is exactly 0 there) -- gives a naturally rounded cap,
             * not a sharp/independent window */
            dist = PetscSqrtReal(SQ(x - (xs - Rs)) + SQ(y));
        } else if (x > xs + Rs) {
            dist = PetscSqrtReal(SQ(x - (xs + Rs)) + SQ(y));
        } else {
            /* inside the segment: perpendicular distance to the floor curve's
             * local tangent line (good approximation for a gently-curving
             * bump); matches the endpoint formula continuously at x = xs+-Rs
             * since slope -> 0 there too */
            PetscReal slope = WallBottomDeriv(user, x);
            dist = (y - y_bot) / PetscSqrtReal(1.0 + SQ(slope));
        }

        ice += 0.5 - 0.5 * PetscTanhReal((tc * ts) * (dist / ts - 1.0));
    }
    return ice;
}

/* Laterally-windowed flat ice layers. */
static PetscReal IceFlatField(PetscReal x, PetscReal y, const AppCtx *user)
{
    const PetscReal tc  = PhaseSlope(user);
    PetscReal       ice = 0.0;

    for (PetscInt k = 0; k < user->n_ice_flats; k++) {
        PetscReal xf   = user->ice_flat_x[k];
        PetscReal Rf   = user->ice_flat_R[k];
        PetscReal Hf   = user->ice_flat_height[k];
        PetscReal dlat = PetscAbsReal(x - xf) / Rf;                /* =1 at edge */
        PetscReal w    = 0.5 - 0.5 * PetscTanhReal((tc * Rf) * (dlat - 1.0));
        /* flat threshold, same interface width as everywhere else (no
         * R-dependent sharpening) */
        PetscReal flat = 0.5 - 0.5 * PetscTanhReal(tc * (y - Hf));
        ice += w * flat;
    }
    return ice;
}

/* The multi-grain ice field: grains plus the optional conformal bodies. */
static PetscReal ShapeMultiGrains2D(PetscReal x, PetscReal y, const AppCtx *user)
{
    PetscReal ice = user->ic_grain_union ? GrainFieldUnion(x, y, user)
                                         : GrainFieldAdditive(x, y, user);
    ice += IceShellField(x, y, user);
    ice += IceFlatField(x, y, user);
    return ice;
}


/* =========================================================================
 * 5. PUBLIC ENTRY POINTS
 * =========================================================================*/

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

    ierr = FillIC1D(iga, U, user, ShapeSlab1D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}



/* -------------------------------------------------------------------------
 * 2D Ice Slab Initial Condition
 *
 *   - Ice: a slab spanning the full width of the domain, from y = 0 up to
 *          y = Ly/2. The top half is pure air. Exactly 50% ice by volume.
 *   - Temperature: temp0 + grad_temp0 . (x - L/2)
 *   - Vapor: hum0 * rho_vs(T) in air, saturated in ice
 *
 * This is the analytic benchmark for the k_eff homogenization (-keff): for a
 * two-phase layered medium the effective conductivity is known in closed form,
 *   parallel to the layers      k = arithmetic mean = (k_ice + k_air)/2
 *   perpendicular to the layers k = harmonic mean   = 2/(1/k_ice + 1/k_air)
 * which at k_ice = 2.29, k_air = 0.02 gives 1.155000 and 0.0396538 W/m/K.
 *
 * THE PERIODIC SEAM. Under -periodic 1 the domain wraps in y, so a slab
 * occupying [0, Ly/2] has TWO interfaces: the obvious one at y = Ly/2, and a
 * second at the y = 0 == Ly seam. That second interface is physically real and
 * must be resolved like the first. Writing the profile with a single tanh,
 *     phi = 0.5*(1 - tanh(0.5*(y - Ly/2)/eps)),
 * gives phi(0) = 1 and phi(Ly) = 0 -- a discontinuity across the seam. It looks
 * right in a plot and is wrong: the jump is unresolved by any mesh, and it
 * corrupts the perpendicular (harmonic) component of k_eff, which is precisely
 * the component the benchmark is checking. The same bug is on record in the
 * standalone homogenization code (effective_thermal_cond/CHANGES.md:40-56).
 *
 * So phi is built from the SIGNED DISTANCE to the nearest interface, the same
 * construction -ic_grain_union uses for grains:
 *
 *     phi = 0.5 - 0.5*tanh(0.5*sdf/eps),   sdf < 0 inside ice
 *
 * with sdf measured to whichever of the two interfaces is closer. That form is
 * exactly 50% ice for any eps: the ice and air slabs are both Ly/2 thick, so
 * sdf(y + Ly/2) = -sdf(y), hence phi(y + Ly/2) = 1 - phi(y), and the two halves
 * of the integral sum to Ly/2 identically. (A naive difference of two tanh
 * profiles does NOT have this property -- it is short by eps*ln(2), i.e. ~2.8%
 * of the ice at eps = 2e-5 on a 1 mm cell, which would show up as a bias in
 * k_parallel that looks like a solver error.)
 *
 * Under -periodic 0 there is only the one interface at y = Ly/2, and phi -> 1
 * at y = 0; that branch is also exactly 50% by antisymmetry about y = Ly/2.
 * -------------------------------------------------------------------------*/
PetscErrorCode FormInitialIceSlab2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD,
                "--- INITIAL CONDITIONS (2D Ice Slab: ice on [0, Ly/2], air above, "
                "%s) ---\n",
                (user->periodic == 1)
                    ? "periodic in y: interfaces at Ly/2 AND the y=0 seam"
                    : "non-periodic: single interface at Ly/2");

    user->n_act = 0;

    ierr = FillIC2D(iga, U, user, IC_COORD_UNIFORM, ShapeIceSlab2D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialSingleIceGrain2D
 *
 * Single pure ice circle (no sediment core) centered in the domain.
 *
 * Parameters: user->RCice (grain radius), Lx, Ly, eps, temp0, hum0.
 * =========================================================================*/
PetscErrorCode FormInitialSingleIceGrain2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx = user->Lx, Ly = user->Ly, RCice = user->RCice;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D single ice grain) ---\n"
        "  centre = (%.4e, %.4e) m,  RCice = %.4e m\n",
        0.5 * Lx, 0.5 * Ly, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice must be > 0 (got %.2e)", RCice);
    if (RCice >= 0.5 * PetscMin(Lx, Ly))
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Grain radius %.2e exceeds half the domain — increase domain or reduce RCice", RCice);

    ierr = FillIC2D(iga, U, user, IC_COORD_UNIFORM, ShapeSingleIceGrain2D); CHKERRQ(ierr);
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

    const PetscReal Lx = user->Lx, RCice = user->RCice;

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (1D single ice grain) ---\n"
        "  centre = %.4e m,  RCice = %.4e m\n", 0.5 * Lx, RCice);

    if (RCice <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice must be > 0 (got %.2e)", RCice);
    if (RCice >= 0.5 * Lx)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "Grain radius %.2e exceeds half the domain — increase Lx or reduce RCice", RCice);

    ierr = FillIC1D(iga, U, user, ShapeSingleIceGrain1D); CHKERRQ(ierr);
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
 * Physical coordinates: for a plain rectangular domain the wall baselines are
 * flat, so IC_COORD_WALLS reduces to (Lx*i/mx, Ly*j/my) as usual. For the
 * -geom_file sediment-grain geometry (build_geometry_sediment_grain.py) the
 * bottom edge is raised by a C-infinity bump y=g(x) and the surface is ruled
 * to the flat top edge, so y_phys = g(x) + (j/my)*(Ly - g(x)); x_phys =
 * Lx*i/mx exactly by construction (both curves share the u<->x
 * parametrization). geom_bump_R must match that script's R_sed (bump
 * half-width == height) for the IC to align with the actual geometry.
 * =========================================================================*/
PetscErrorCode FormInitialTwoIceGrainsBoundary2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscReal Lx = user->Lx, Ly = user->Ly;
    const PetscReal R0 = user->RCice0;   /* left grain (x=0), smaller  */
    const PetscReal R1 = user->RCice1;   /* right grain (x=Lx), larger */

    PetscPrintf(PETSC_COMM_WORLD,
        "--- INITIAL CONDITIONS (2D two ice grains, boundary-centered) ---\n"
        "  grain 0: centre = (%.4e, %.4e) m,  RCice0 = %.4e m\n"
        "  grain 1: centre = (%.4e, %.4e) m,  RCice1 = %.4e m\n",
        0.0, 0.5 * Ly, R0, Lx, 0.5 * Ly, R1);

    if (R0 <= 0.0 || R1 <= 0.0)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0 and RCice1 must be > 0 (got %.2e, %.2e)", R0, R1);
    if (R0 >= Lx || R1 >= Lx)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0/RCice1 must be < Lx (got %.2e, %.2e, Lx=%.2e)", R0, R1, Lx);
    if (R0 > 0.5 * Ly || R1 > 0.5 * Ly)
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                "RCice0/RCice1 must be <= Ly/2 (got %.2e, %.2e, Ly/2=%.2e)", R0, R1, 0.5*Ly);

    ierr = FillIC2D(iga, U, user, IC_COORD_WALLS, ShapeTwoIceGrainsBoundary2D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * FormInitialMultiGrains2D
 *
 * N ice grains, each a tanh distance profile from a (cx,cy)/R given by
 * -ice_grain_cx/-ice_grain_cy/-ice_grain_R, summed and clamped to [0,1] --
 * the same construction as FormInitialTwoIceGrainsBoundary2D generalized
 * to an arbitrary number of grains. Optional extra ice bodies (conformal
 * shells, flat layers) add on top; see the feature fields above.
 *
 * Physical coordinates use the Greville abscissae and the wall curves, i.e.
 * the bottom edge is the sum of -sed_grain_x/-sed_grain_R bump humps (or the
 * single -geom_bump_R bump if those aren't set), matching
 * build_geometry_multi_grain.py / build_geometry_sediment_grain.py.
 * =========================================================================*/
PetscErrorCode FormInitialMultiGrains2D(IGA iga, Vec U, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const PetscInt ng = user->n_act;

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

    ierr = FillIC2D(iga, U, user, IC_COORD_GREVILLE, ShapeMultiGrains2D); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


/* =========================================================================
 * 6. PACKING LOADER
 * =========================================================================*/

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
