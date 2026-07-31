#include "keff.h"

/* ---------------------------------------------------------------------------
 * keff_field.c -- project the phase field onto the corrector mesh's quadrature
 * points.
 *
 * THE GAUSS-POINT INDEXING CONVENTION (project-wide; adopted 2026-07-31)
 * ---------------------------------------------------------------------
 * Any array holding one value per local quadrature point is indexed by
 *
 *     idx = point->index + point->count * point->parent->index
 *
 * i.e. (element index within this rank's element loop) * (points per element)
 * + (point index within the element). Never by a sequential ++ counter.
 *
 * The reason is not style. The read side of this array is a PetIGA FORM
 * CALLBACK (IGASetFormMatrix / IGAComputeScalar): PetIGA owns the loop and
 * hands us one point at a time, so there is no correct place to keep a counter.
 * The formula is stateless -- every point derives its own index from its own
 * identity -- which is the only thing that works inside a callback.
 *
 * Its correctness condition is that point->count is the same for every element
 * being indexed. That is NOT automatic. From PetIGA (src/petigaelem.c:812-818):
 *
 *     nqp = IGA_Quadrature_SIZE(BD,ID,NQ);
 *     if (element->atboundary) { axis = element->boundary_id/2;
 *                                nqp /= NQ[axis]; NQ[axis] = 1; }
 *     element->iterator->count = nqp;
 *
 * so point->count is per-element and is REDUCED on boundary passes; and
 * IGAElementNextForm (petigaelem.c:427-440) can visit one element several
 * times -- once per active boundary face, then once for the interior -- with
 * parent->index unchanged. Boundary passes are therefore excluded everywhere
 * by an explicit `if (!point->atboundary)` guard, and KeffCheckGaussLayout
 * asserts the uniformity of point->count once at startup.
 *
 * (In the -keff configuration boundary passes cannot fire anyway: the mesh is
 * fully periodic, and IGAClone does not copy iga->form, so PetIGA's
 * boundary-visit table is empty. The guard and the assertion are there so this
 * stays correct if either of those ever changes.)
 * ------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------
 * KeffProjectIce
 *
 * Evaluate phi (dof 0 of the 3-dof solution) at every local quadrature point of
 * the SOLVER's IGA and store it in kc->ice[], which the cell-problem assembly
 * later reads while looping the CLONED IGA. The two loops agree because
 * IGAClone copies elem_start/elem_width and the quadrature rule, which
 * KeffCreate asserts, and because both sides use the formula above.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffProjectIce(AppCtx *app, Vec U)
{
  PetscErrorCode     ierr;
  KeffCtx           *kc = app->keff;
  Vec                localU;
  const PetscScalar *arrayU;
  IGAElement         element;
  IGAPoint           point;
  PetscScalar       *UU;
  PetscInt           nvisited = 0;

  PetscFunctionBegin;

  ierr = IGAGetLocalVecArray(app->iga, U, &localU, &arrayU); CHKERRQ(ierr);
  ierr = IGABeginElement(app->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(app->iga, element)) {
    ierr = IGAElementGetValues(element, arrayU, &UU); CHKERRQ(ierr);
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      if (point->atboundary) continue;
      {
        PetscScalar sol[3];
        PetscInt    idx = point->index + point->count * point->parent->index;

        if (idx < 0 || idx >= app->ngp)
          SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
                  "k_eff: Gauss-point index %d out of range [0,%d). Element %d, "
                  "point %d of %d. The quadrature layout is not what the flat "
                  "array assumes.", (int)idx, (int)app->ngp,
                  (int)point->parent->index, (int)point->index, (int)point->count);

        ierr = IGAPointFormValue(point, UU, &sol[0]); CHKERRQ(ierr);
        kc->ice[idx] = PetscRealPart(sol[0]);
        nvisited++;
      }
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(app->iga, &element); CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(app->iga, U, &localU, &arrayU); CHKERRQ(ierr);

  /* Every entry the assembly will read must have been written. A short count
   * means some element contributed fewer points than the flat layout reserves
   * for it, which would leave stale phi behind. */
  if (nvisited != app->ngp)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
            "k_eff: projected %d Gauss points but the flat array holds %d. The "
            "quadrature layout is not uniform across elements.",
            (int)nvisited, (int)app->ngp);

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffCheckGaussLayout
 *
 * Assert, once, that point->count is identical on every element of the CLONE
 * and equals the parent's -- the precondition for the index formula above.
 * Called from KeffCreate.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffCheckGaussLayout(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;
  IGAElement     element;
  IGAPoint       point;
  PetscInt       count0 = -1, nelem = 0, ntotal = 0;

  PetscFunctionBegin;

  ierr = IGABeginElement(kc->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(kc->iga, element)) {
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      if (point->atboundary) continue;
      if (count0 < 0) count0 = point->count;
      else if (point->count != count0) {
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        ierr = IGAEndElement(kc->iga, &element); CHKERRQ(ierr);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
                "k_eff: quadrature points per element are not uniform (%d vs %d "
                "at element %d). The flat Gauss-point index formula "
                "point->index + point->count*point->parent->index assumes they "
                "are; see keff_field.c.",
                (int)point->count, (int)count0, (int)point->parent->index);
      }
      ntotal++;
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    nelem++;
  }
  ierr = IGAEndElement(kc->iga, &element); CHKERRQ(ierr);

  if (kc->debug_phibar) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "  [keff] layout: %d elements x %d pts/elem = %d local Gauss points (ngp=%d)\n",
      (int)nelem, (int)count0, (int)ntotal, (int)app->ngp); CHKERRQ(ierr);
  }

  /* The clone must see exactly the population the solver IGA sized ngp for. */
  if (ntotal != app->ngp)
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
            "k_eff: cloned IGA has %d local Gauss points but the solver IGA "
            "reported %d (%d elements x %d points). The two meshes do not share "
            "a quadrature layout, so phi cannot be handed between them.",
            (int)ntotal, (int)app->ngp, (int)nelem, (int)count0);

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * Ice moments: [ INT phi dV, INT phi*x dV, INT phi*y dV, (INT phi*z dV) ]
 *
 * WHY MOMENTS AND NOT JUST THE MEAN. The obvious cross-check is the mean ice
 * fraction, INT phi dV, computed on both meshes. It is worthless for this
 * purpose: it is a plain weighted sum, so ANY permutation of the index mapping
 * that preserves the multiset of (phi, quadrature weight) pairings leaves it
 * exactly unchanged. On a uniform mesh that includes every whole-element shift
 * (all elements carry identical weights) and every within-element reversal
 * (Gauss weights are palindromic, so w[n-1-q] == w[q]). Both were confirmed to
 * pass a mean-only check with a deliberately scrambled index -- the check had
 * no teeth at all.
 *
 * The FIRST MOMENTS fix this because they pair each phi with its own physical
 * coordinate, obtained from the point itself via IGAPointFormGeomMap. Move phi
 * to the wrong point and the ice centroid moves. Normalised by INT phi dV it is
 * the centre of mass of the ice, in metres, which also makes a failure legible
 * rather than an abstract residual. Measured: a one-slot cyclic shift within
 * each element moves the centroid by 1.3e-3 relative while leaving the mean at
 * 2.3e-15 -- i.e. the classic off-by-one is caught only by the moments.
 *
 * SECOND MOMENTS close the remaining blind spot. PetIGA orders the (p+1)^dim
 * points of an element as a tensor grid, so index q -> count-1-q is exactly the
 * reflection through the element centre, x_{n-1-q} = 2*c_e - x_q. Under it the
 * first moment per element becomes 2*c_e*M0_e - M1_e, which cancels to within
 * the quadrature error of a smooth field -- a central reflection sails through
 * a mean-and-centroid check. The second moment does not cancel, because
 * (2c - x)^2 is not -(x^2) plus a constant multiple of x.
 *
 * So the check compares 1 + 2*dim numbers: the mean, the centroid, and the
 * per-axis mean square. All three are computed twice, once from kc->ice[] on
 * the clone and once from U on the solver mesh, sharing no code.
 * ------------------------------------------------------------------------- */

/* On the CLONE, reading phi back out of kc->ice[]. */
static PetscErrorCode KeffMomentCloneIntegrand(IGAPoint point, const PetscScalar U[],
                                               PetscInt n, PetscScalar S[], void *ctx)
{
  KeffCtx  *kc = (KeffCtx *)ctx;
  PetscReal x[3] = {0.0, 0.0, 0.0};

  PetscFunctionBegin;
  (void)U; (void)n;
  if (!point->atboundary) {
    PetscInt  idx = point->index + point->count * point->parent->index;
    PetscReal phi = kc->ice[idx];
    IGAPointFormGeomMap(point, x);
    S[0] += phi;
    for (PetscInt d = 0; d < kc->dim; d++) {
      S[1 + d]           += phi * x[d];
      S[1 + kc->dim + d] += phi * x[d] * x[d];
    }
  }
  PetscFunctionReturn(0);
}

/* On the SOLVER IGA, reading phi straight out of the 3-dof solution vector.
 * Shares no code and no array with the projection, which is the point. */
static PetscErrorCode KeffMomentSolverIntegrand(IGAPoint point, const PetscScalar U[],
                                                PetscInt n, PetscScalar S[], void *ctx)
{
  AppCtx     *app = (AppCtx *)ctx;
  PetscReal   x[3] = {0.0, 0.0, 0.0};
  PetscScalar sol[3];

  PetscFunctionBegin;
  (void)n;
  if (!point->atboundary) {
    PetscReal phi;
    IGAPointFormValue(point, U, &sol[0]);
    IGAPointFormGeomMap(point, x);
    phi = PetscRealPart(sol[0]);
    S[0] += phi;
    for (PetscInt d = 0; d < app->dim; d++) {
      S[1 + d]            += phi * x[d];
      S[1 + app->dim + d] += phi * x[d] * x[d];
    }
  }
  PetscFunctionReturn(0);
}

/* m[0]            = mean ice fraction
 * m[1 .. dim]     = ice centroid coordinates [m]
 * m[1+dim .. 2dim]= ice mean square per axis [m^2]
 * from_clone selects which mesh (and hence which phi source) is used; U is the
 * 3-dof solution, needed only on the solver path. */
PetscErrorCode KeffIceMoments(AppCtx *app, PetscBool from_clone, Vec U, PetscReal m[7])
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;
  PetscScalar    s[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const PetscInt n    = 1 + 2 * kc->dim;
  PetscReal      m0;

  PetscFunctionBegin;

  if (from_clone) {
    ierr = IGAComputeScalar(kc->iga, kc->T[0], n, &s[0],
                            KeffMomentCloneIntegrand, kc); CHKERRQ(ierr);
  } else {
    ierr = IGAComputeScalar(app->iga, U, n, &s[0],
                            KeffMomentSolverIntegrand, app); CHKERRQ(ierr);
  }

  m0   = PetscMax(PetscRealPart(s[0]), 1.0e-300);
  m[0] = PetscRealPart(s[0]) / kc->vol;
  for (PetscInt d = 0; d < kc->dim; d++) {
    m[1 + d]           = PetscRealPart(s[1 + d])           / m0;
    m[1 + kc->dim + d] = PetscRealPart(s[1 + kc->dim + d]) / m0;
  }

  PetscFunctionReturn(0);
}
