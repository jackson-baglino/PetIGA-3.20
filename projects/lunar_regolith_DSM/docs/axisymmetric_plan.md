# Axisymmetric (r–z) Mode — Implementation Plan

> ✅ **CURRENT (banner 2026-07-22).** Applies to the two-phase model as implemented. Predates the rename to `sublimation_pf` — read any `projects/permafrost` path as `projects/sublimation_pf`.


*2026-07-10. Goal: true 3D two-sphere sintering physics at 2D cost, for the
Molaro et al. (2019) Fig. 11 validation. Status: PLAN — review before
implementation.*

## 0. Why this works and what it buys

The Molaro two-grain configuration is exactly axisymmetric about the line
joining the grain centers. In cylindrical coordinates (z along that axis,
r transverse, no θ-dependence), every volume integral becomes

```
∫∫∫ f dV  =  2π ∫∫ f(r,z) · r  dr dz
```

and for scalar diffusion terms the weak form keeps its Cartesian shape —
`∫ k ∇u·∇v r dr dz` with the ordinary (r,z) gradients — because the extra
`(1/r)∂r(r·)` of the cylindrical Laplacian is generated automatically by the
r-weighted integration by parts. **The entire formulation change is:
multiply every integrand by the quadrature point's r-coordinate.** The
missing second principal curvature (spheres vs cylinders: 2/R vs 1/R;
the neck saddle's azimuthal term) then emerges from the geometry with no
further modeling.

Cost: the (z,r) domain is HALF the current planar one (r ≥ 0, axis on the
boundary), so the Molaro mesh drops from 866k to ~440k nodes. Expected
physics effect: ~2–3× faster neck growth (the remaining deficit after the
eps sweep exonerated discretization).

## 1. Formulation choices

- **Coordinates:** keep the existing mesh axes; x ≡ z (axis direction),
  y ≡ r (radial). The symmetry axis is the y = 0 boundary.
- **Grains:** centered ON the axis (`-ice_grain_cy 0,0`). A circle in the
  (z,r) plane centered on the axis IS a sphere: the in-plane distance
  `sqrt((z−cz)² + r²)` equals the 3D distance. The existing multi-grain IC
  therefore produces exact spheres with no changes to its math.
- **The 2π factor** is dropped from the residual/Jacobian (a global constant
  scaling of R = 0 changes nothing) but INCLUDED in the monitor integrals so
  reported masses/areas are true 3D quantities (units change: TOT_ICE
  becomes m³).
- **Axis boundary (r = 0):** natural Neumann is automatically correct — the
  r-weight sends the boundary measure to zero; no special treatment. Gauss
  points are element-interior, so no integrand is ever evaluated at exactly
  r = 0 (no 0/0, no zero-diagonal rows: every basis function's support
  contains points with r > 0).
- **Other boundaries:** Neumann as today; `-periodic 0`.

## 1b. Where the third curvature mode comes from (the |grad phi|^2 question)

The gradient energy density |grad phi|^2 is local and geometry-blind; the
curvature physics enters through the VARIATIONAL BALANCE — i.e. through the
volume measure the gradient term is integrated against, because
curvature-driven motion exists exactly because moving a curved interface
changes its area, and the measure is what knows about area. With
dV = 2*pi*r dr dz, integrating the weak-form term  ∫ grad(N)·grad(phi) r dr dz
by parts yields

    -(1/r) d/dr( r d(phi)/dr ) - d2(phi)/dz2
      = -( d2/dr2 + d2/dz2 )phi  -  (1/r) d(phi)/dr
        [in-plane curvature]       [azimuthal curvature mode]

The (1/r)*d(phi)/dr term IS the second principal curvature — generated
automatically by the r-weight, with the |grad phi|^2 term untouched.
Sphere check (must give kappa = 2/R everywhere): at the equator (normal
radial, r = R) the in-plane part gives phi'/R and (1/r)phi_r = phi'/R ->
2/R; at the pole (r -> 0, phi_r -> 0) the limit of (1/r)phi_r is
phi_rr = phi'/R -> again 2/R. Correct and regular at the axis.

At the neck this is the dominant 2D-vs-3D difference: kappa_neck becomes
1/x - 1/rho (the azimuthal neck-circle term is new physics 2D cannot
represent), and grain surfaces go 1/R -> 2/R, roughly doubling the
source-sink driving difference. The ONLY code path needing a hand-added
azimuthal term is Curvature() (explicit kappa from grad/Hessian for d0_GT)
— off and out of scope.

## 2. Code touch points (in implementation order)

1. **`include/NASA_types.h`** — add `PetscBool axisym;` to AppCtx.
2. **`src/permafrost2.c`** — `-axisym` option (default 0 = planar, fully
   backward compatible); print mode in the parameter header; pass through
   context as usual.
3. **`src/assembly.c`** (`Residual`, `Jacobian`) — at the top of each point
   evaluation:
   ```c
   PetscReal rw = 1.0;
   if (user->axisym) { PetscReal xphys[3];
       IGAPointFormPoint(pnt, xphys); rw = xphys[1]; }
   ```
   and multiply every `R[a][*]` and `J[a][*][b][*]` accumulation by `rw`
   (single multiplicative factor at the end of each expression — no term
   changes shape). API verified: `IGAPointFormPoint` (petiga.h:725) returns
   the mapped physical coordinates.
4. **`src/monitoring.c`** — the domain-integral point loop in `Monitor()`
   gets the same `rw` (times 2π) so TOT_ICE/TOT_AIR/TOT_RHOV/TOTAL_MASS and
   the I-A interface integral are genuine 3D volumes/areas; mass-conservation
   checks then remain meaningful. `Sigma0/flag_Tdep` branch: same weight if
   ever re-enabled (touch now, costs nothing). Interface-CFL limiter:
   untouched (pointwise |dphi|, geometry-agnostic).
5. **`src/initial_conditions.c`** — relax the multi-grain validation that
   currently rejects grains poking outside the domain in y (on-axis grains
   have cy = 0 and legitimately "extend" to r < 0; the field for r>=0 is all
   that exists). Guard: only when `-axisym 1`.
6. **`inputs/geometry/2D_molaro_axisym.opts`** — z-extent as today
   (Lx = 384.3–385 µm, per-eps calibrated overlap unchanged: the (z,r)
   cross-section union math is identical to the planar case, so the SAME
   calibrated x0 values apply); r-extent Ly_r = R1 + pad = 121 µm;
   `-ice_grain_cy 0,0`; Nx × Ny_r ≈ 1173 × 369 at s050 eps (≈ 433k nodes).
7. **`postprocess/neck_width.py`** — `--axisym` flag: neck width =
   2 × (height of the outermost phi = 0.5 crossing above y = 0) at the waist
   column; grain-center peaks logic unchanged (w(x) is then the neck
   diameter directly).
8. **Docs/log** — record formulation + validation results.

Deliberately NOT in scope:
- The GT curvature correction (`d0_GT`) stays 0. If it is ever re-enabled in
  axisym, `Curvature()` must add the azimuthal term −(∂r phi)/(r|∇phi|);
  flagged in comments, not implemented.
- Sediment bumps / shells / flats in axisym (they're planar-specific); the
  IC guard only relaxes for plain grains.

## 3. Validation sequence (each step gates the next)

- **V0 — volume check (free):** at t = 0, TOT_ICE must equal
  (4/3)π(R0³+R1³) − (lens overlap volume) to ~0.1%. Analytic value computed
  alongside; instant sanity check on the r-weight + 2π plumbing.
- **V1 — equilibrium hold (cheap, local):** single on-axis sphere,
  saturated, short run: phi bounds clean, TOT_ICE constant, no spurious
  drift near the axis. Catches axis-boundary mistakes.
- **V2 — sphere-shrinkage rate (the quantitative gate):** small on-axis
  sphere with humidity < 1: kinetic-limited rate dR/dt = −(rhov_s − rhov)/
  (rho_ice·beta_scaled)·(...) has an analytic 3D form; also comparable to the
  same test run planar (expect the 2/R vs 1/R curvature factor where
  capillarity drives). Verifies the axisym physics, not just the plumbing.
- **V3 — Molaro rerun (HPC):** 2D_molaro_axisym + molaro_T-20_h1.00,
  neck_width.py --axisym. Compare against the planar s050 curve (+17.8%):
  the 3D-curvature hypothesis predicts ~2–3× that growth. Then the residual
  gap to the experiment's ~2× is the alpha_c calibration, run as a short
  ladder on this same geometry.

## 4. Risks / open points

- `IGAPointFormPoint` cost per quadrature point: one small matvec; if
  profiling shows it matters, hoist to using `pnt->mapX[0]` directly.
- Any postprocess script that assumes the grains sit at mid-height
  (cy = Ly/2) needs the --axisym switch (only neck_width.py is in active
  use; others checked as encountered).
- Monitor units change in axisym runs (volumes not areas) — mass_plots
  scripts read the columns generically, but percent-change logic is
  unit-free, so expected to pass through; verify in V1.
- The 2-hour experiment at ~2–3× faster kinetics: interface-CFL limiter
  will throttle more near the neck; expect modestly more steps than the
  planar run (still ≪ 1000 → per-step output preserved).

## 5. Effort estimate

Items 1–3: ~1 hour including rebuild. Items 4–7: ~1–2 hours. V0–V2
validations: local, same day. V3: one HPC submission.
