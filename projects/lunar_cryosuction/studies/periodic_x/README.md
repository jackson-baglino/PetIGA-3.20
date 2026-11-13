# periodic_x: one pore cell, periodic along x (plan type 5, redesigned 13 Nov)

**Question.** A lens in a converging-diverging pore should migrate by
curvature alone (no temperature gradient) toward the throat and stop there.
Fixed vapor at the pore ends adds or removes water and biases the motion, so
the pore is one periodic cell: throat in the middle, widest at the seam,
regolith walls on y = 0 and y = Ly with the contact-angle condition. The
lens' two faces exchange vapor through the neighboring period and the water
is conserved.

## What had to change

| Where | Change | Size |
|---|---|---|
| `src/cryosuction_main.c`, `include/NASA_types.h` (44f8e40f) | `-periodic_x 1`: `IGAAxisSetPeriodic` on the x axis only, for the built-in rectangle and **before** `IGARead` for `-geom_file`; check the file's x knot vector is periodic; Nx from the element count; skip x faces in the Dirichlet loops; refuse x walls, x reservoirs, an x temperature gradient, `-periodic 1` together, 1D | ~110 lines, all behind the flag |
| `src/initial_conditions.c` (44f8e40f) | shape ICs on the built-in rectangle take x from the periodic Greville points; the multi-grain minimum image wraps x only | ~20 lines |
| `preprocess/build_ice_layout.py` (f4d5fcda) | `"periodic": true` mesh + phi file; `"walls": {"profile": ...}` cell | ~150 lines |
| `postprocess/pplib.py`, `plot_fields.py` (dc02cf4d) | read the N-column vector of a periodic run onto the N + p control net | ~40 lines |

Nothing in the residual, the Jacobian or the wall term changed: a periodic
axis is a different function space (PetIGA wraps the last p basis functions
onto the first p, `petigaaxis.c`, nnp = n - C), and the assembly loops over
elements as before. The user's guess was right: it is a boundary-condition
change, not a solver change.

Why it works with a geometry file (PetIGA source, `src/petigaio.c` IGALoad,
`src/petiga.c` IGASetUp_Stage1): IGALoad resets the IGA but not the axes'
periodic flag; it sets the knots and then sizes the node grid from the flag;
the control net is stored on `geom_sizes = n + 1` points either way. So a
periodic mesh file is an ordinary igakit file whose x knot vector is the
periodic one and whose last p control columns repeat the first p shifted by
Lx. No flag is stored in the file; the solver must set it before reading.
igakit writes such a file as is (`NURBS([U_periodic, V], control)`; its
`unclamp` is not needed because the builder places the control points
itself).

## Wall profile for type 5

gap g(x) about the throat x_t = Lx/2, walls symmetric about the axis:

- `sinhyp` (proposed): g = sqrt(wt^2 + (c Lx/pi)^2 sin^2(pi xi/Lx)),
  xi = x - x_t, c = pi sqrt(W^2 - wt^2)/Lx. It is the hyperbolic throat
  sqrt(wt^2 + c^2 xi^2) of plan_estimates section 8b with xi replaced by
  (Lx/pi) sin(pi xi/Lx): the same to O(xi^2/Lx^2) near the throat, periodic,
  C-infinity, zero slope at the seam (g = W).
- `cos`: g = wt + (W - wt)(1 - cos(2 pi xi/Lx))/2.

Tradeoff (W 150, wt 40, Lx 800 um): both reach the same steepest wall angle
(12.2 deg), but sinhyp gets there 121 um from the throat (cos: 200 um) and
its throat is 2.4x sharper (g'' 8.1e3 vs 3.4e3 1/m), so the wall angle is
4.4 deg at 20 um and 10.1 deg at 60 um from the throat (cos: 1.9 and 5.6).
A lens near the throat therefore feels a stronger curvature difference in
the sinhyp cell, and the sinhyp cell keeps the plan's hyperbolic throat.
Any smooth periodic profile has zero slope at the throat and the seam, so
the drive vanishes at both; a lens started at the seam sits on an unstable
balance point. The cost of sinhyp: the steep flank is closer to the throat,
so theta + omega reaches its largest value (72 deg at theta 60) nearer the
rest position.

## Inputs

`studies/periodic_x/scripts/make_px_inputs.py` writes five layouts
(`preprocess/layouts/px_*.json`), builds them into
`inputs/geometry/periodic_x/` and writes
`inputs/experiment/periodic_x/px_T-20_iso_theta60_2d_ac1e-2.opts`. The test
cell is small (Lx 120 um, W 40 um, wt 16 um, eps 1.5 um, 13.7 k DOF); the
production type-5 cell is built the same way with the type-5 sizes.

## Verification

`verification/verify_periodic_x.sh` (9 runs through
`scripts/Studio/run_cryosuction.sh`, one at a time); results in
`verification/verify_periodic_x.csv` and `.log` (2026-11-13: 56/56 PASS).
G1 of `studies/contact_angle/verification/verify_wall_bc.sh`: PASS (base
83590222).

| Gate | Measured | Tolerance |
|---|---|---|
| a mesh/IC: knots = PetIGA's, wrapped columns = first + Lx, wall and slope jump at the seam, phi file = N nodes, seam layout = mid layout shifted N/2 | 0 to 1.5e-14 | 1e-14 to 1e-9 |
| b flat periodic channel, 2 d: lens center drift | 4.5e-15 m | 1.5e-10 m (1e-4 eps) |
| b water Q change (sealed open channel: 4.5e-10) | 4.5e-10 | 2x reference, floor 1e-8 |
| c1 lens across the seam vs half-period shift: area, center | 3.7e-10, 6e-15 m | 1e-6, 1.5e-10 m |
| c2 cd cell with the throat on the seam vs shift: area, center | 6.8e-8, 1.5e-13 m | 1e-6, 1.5e-10 m |
| d cd cell: lens center toward the throat | +1.543 um (17.7 -> 16.2 um off) | > 0 |
| d water Q change | 8.7e-8 | 1e-6 (Newton tolerance) |
| f built-in rectangle vs geom_file mesh: area, center | 2.4e-11, 1e-15 m | 1e-6, 1.5e-10 m |
| g 4 ranks vs 1 rank: area, center | 1.3e-7, 2.6e-13 m | 1e-6, 1.5e-10 m |
| h tight Newton: Q change; c2 area, center | -5.9e-13; 5.3e-16, 2.4e-18 m | 10x below default; 1e-8, 1.5e-12 m |

Gate h is why the d and c2 tolerances are the Newton tolerance: with
`-snes_atol 1e-6` the vapor equation is accepted after one iteration, and the
water change and the shift differences then sit at 1e-8 to 1e-7; tightening
the solve takes them to 1e-13 and 1e-16, so they are solver tolerance, not
the seam.

Limits of these checks: the flat-channel lenses (b, c1) are symmetric, so
zero drift there is also required by symmetry; the drift gate tests the seam
only in that a seam artifact would break the symmetry. The cd run is 2 days,
long enough for the sign (1.5 um, ~0.8 um/d at 17 um from the throat) but not
for the stop at the throat.
