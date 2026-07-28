# Phase 4 — safety-factor and xi_v calibration

**Question.** What is the cheapest discretisation that still reproduces the
physics? The production matrix is 200 runs, so a factor of 2 in `eps` is a
factor of 4 in cost across the whole study.

## Where each arm runs

All calibration domains are **<= 0.5 mm**. Cost is (DOF) x (steps), and STEP
COUNT is what separates the two families -- not mesh size:

| family | steps | DOF | wall, 6 ranks | where |
|---|---|---|---|---|
| Molaro neck (test A) | ~125 (t_final 2 h, dtmax 100 s) | 0.4-2.3 M | 0.6 h total | **local** |
| 7-day calibration (all arms) | ~3075 (t_final 7 d, dtmax 200 s) | 0.2-6 M | ~55 h serial | **HPC** |

Two numbers behind that, both measured rather than assumed:

- **Memory** (`preprocess/estimate_memory.py`): ILU(3) fill dominates, ~3x the
  Jacobian. On 0.5 mm domains every arm fits in 64 GB except `ripen_s025`
  (10 GB) and `pack_s025` (21 GB). An earlier 1 mm `pack_s025` needed **83 GB**
  and exhausted the workstation -- always run the estimator before a local job.
- **Time**: 7.5 s/step for 1.5 M DOF on 6 ranks, ~linear in DOF. The adaptive
  stepper needs ~50 steps just to climb from dt = 1e-4 to dtmax = 200 s, so the
  7-day arms are ~3075 steps regardless of how coarse the mesh is.

So the 7-day arms fit in memory but not in patience: they all go to the
cluster.

Locally:

```bash
./scripts/Studio/run_calibration_local.sh --dry-run
./scripts/Studio/run_calibration_local.sh   # 3 Molaro arms, ~0.6 h
```

On the cluster — `heavy` (14 jobs) is everything except the three Molaro arms:

```bash
cd $PETIGA_DIR/projects && git fetch origin && \
  git checkout restructure/enceladus-lunar-split && git pull
cd enceladus_DSM && make
./scripts/HPC/submit_calibration.sh --dry-run heavy
./scripts/HPC/submit_calibration.sh heavy
```

Then, once everything has landed (local + cluster results in one tree):

```bash
./venv_enceladus/bin/python3 postprocess/calibration_report.py \
    <results-root> --tol 0.02 --mass-tol 0.01
```

---

## The two knobs push in opposite directions

**`safety`** — `eps = safety * min(K&P bounds)`. With `R_feat = 2 um` the
kinetic bound binds, so `eps = 2*safety` um exactly. Cost goes as `1/eps^2`:

| arm | safety | eps | ripen (0.5 x 0.25 mm) | pack (0.5 x 0.5 mm) | pack RAM |
|---|---|---|---|---|---|
| s025 | 0.25 | 0.5 um | 1415 x 708 | 1415 x 1415 | 21 G |
| s050 | 0.50 | 1.0 um | 708 x 354 | 708 x 708 | 5.2 G |
| s075 | 0.75 | 1.5 um | 472 x 236 | 472 x 472 | 2.3 G |
| s100 | 1.00 | 2.0 um | 354 x 177 | 354 x 354 | 1.3 G |

**`xi_v`** — this one is counter-intuitive, and the intuition matters. The
solver assembles (`src/assembly.c:21`)

    (1/xi_v) * d(phi_a*rhov)/dt = div(D_v*phi_a*grad rhov) + rho_ice*phi_a_t

so **1/xi_v amplifies STORAGE**. Smaller `xi_v` slows vapour relaxation and
permits a larger `dt` — small is the CHEAP direction, large is the physical one
(`xi_v = 1` is unscaled). The failure mode is at the small end, where the
artificially slowed vapour field stops being quasi-steady:

    tau_vap = L^2 / (xi_v * D_v)

| xi_v | tau_vap (0.5 mm test) | vs 7 d | tau_vap (2 mm production) | vs 28 d |
|---|---|---|---|---|
| 1e-4 | 132 s | 4600x | 2108 s | 1150x |
| 1e-3 | 13 s | 46000x | 211 s | 11500x |
| 1e-2 | 1.3 s | 460000x | 21 s | 115000x |

`tau_vap` grows as `L^2`, so the 2 mm production domain is **16x less
quasi-steady** than the 0.5 mm calibration domain at the same `xi_v`. Read the
`xi_v` result as a LOWER BOUND on what production needs, and re-check the
accepted value at 2 mm before committing to the 200-run matrix. This is the
one result here that does not transfer directly.

Note K&P Eq. 48 (`xi_v <= rho_vs/rho_ice ~ 1e-6`) is a *different* condition —
it is what you need to resolve the true vapour transient, which is exactly what
the temporal scaling is designed to avoid. Exceeding it by ~1000x is
intentional. `comp_eps.py` warns "quasi-steady assumption violated" here, which
is misleading: exceeding Eq. 48 pushes *towards* quasi-steady, not away.

## Three independent tests

A setting must survive all three — each probes a different failure mode.

| test | geometry | observable | probes |
|---|---|---|---|
| **A. sintering** | `2D_molaro_axisym_T-20pair_union_eps{strict,mid,loose}` | neck width vs time | whether necks form at the right rate |
| **B. Ostwald ripening** | `calib_ripen_s*` | large-grain growth at the small grain's expense; ice-mass drift | whether curvature-driven transport survives coarsening |
| **C. packing** | `calib_pack_s*` | SSA(t), k_eff(t) | the actual production observable |

Test A reuses the existing Molaro geometries whose neck width is held FIXED
across `eps` by `preprocess/calibrate_neck_geometry.py`. That matters: without
it, a coarser `eps` would start from a different neck and the comparison would
confound geometry with resolution.

Test B uses the study's grain scale (45 / 90 um), not the legacy 9 / 19 um of
`2D_two_ice_grains_boundary.opts`, so `eps/R` matches production.

Test C runs on 0.5 mm rather than the 2 mm production domain (28 grains vs
~420) — same grain size, contact gap and porosity. That is adequate for a
CONVERGENCE comparison, which contrasts one geometry against itself at
different `eps`; it is NOT an REV, so its absolute k_eff means nothing.

## Acceptance

Adopt the **coarsest** arm that holds every observable within **2%** of the
finest arm, with ice-mass drift under **1%**, in **all three** tests. If tests
disagree, the strictest one wins.

Also confirm before adopting:
- no `DIVERGED_FUNCTION_DOMAIN` (the phase bounds guard fired)
- Peclet `eps*v_n/D_v << 1`
- at the accepted `eps`, the narrowest throats still pass appreciable vapour:
  `phi_a` at the throat centre falls off sharply with `eps` (see `DESIGN.md`
  section 4 — at gap = 4 um it is 0.45 at eps = 1 um but 0.24 at eps = 2 um)

That last point is the one to watch: `s100` is 16x cheaper than `s025` but
halves the conductance of the narrowest throats, so k_eff is the observable
most likely to reject it.

## k_eff for test C

`SSA_evo.dat` gives SSA directly; k_eff needs a second pass over the snapshots:

```bash
cd ../effective_thermal_cond
OUT_ROOT=$SCRATCH/keff_calib ./scripts/run_effective_k_ice_homog.sh <run_dir>
```

then re-run the report with `--keff-root $SCRATCH/keff_calib` to fold it in.
The run directory stages its own `solver.opts` and geometry `.opts`, which the
k_eff script now reads so the mesh matches by construction.
