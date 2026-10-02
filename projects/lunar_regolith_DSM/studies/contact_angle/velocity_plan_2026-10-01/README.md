# Contact angle and interface velocity: figures and talking points

2026-10-01. Theory only; no runs submitted yet. Figures are produced by
`make_theory_figures.py` in this folder (the PNGs are gitignored, so rerun the
script to regenerate them).

The sealed runs already recover Young's angle to about 0.05°, so the wall term
is trusted. This study is about what the contact angle does to the interface
velocity, first in the rectangular channel and then in the wedge.

---

## 1. Setup: one relation for both geometries

![Channel and wedge geometries](fig1_geometry.png)

```
v_n = (σ∞ − d0·χ) / (β + Z)
```

Talking points

- Gibbs–Thomson at the interface, `σ_i = d0·χ + β·v_n`, in series with
  quasi-steady vapour diffusion to the fixed-σ∞ walls,
  `ρ_ice·v_n = D_v·ρ_vs·(σ∞ − σ_i)/length`. Eliminating σ_i gives the relation
  above.
- The contact angle enters only through the curvature χ. The geometry sets χ
  and the diffusive resistance Z.
- Z = K·ℓ in the channel and Z = K·r·ln(r_res/r) in the wedge, with
  K = ρ_ice/(ρ_vs·D_v) = 5.7e10 s/m² at −20 °C.
- At ℓ = 100 µm the diffusive resistance is Z = 5.8e6 s/m. Whether the runs
  are kinetics- or diffusion-limited is set by α_c (next section).
- Sign conventions: θ is measured through the ice; χ > 0 where the ice is
  convex into the vapour; v_n > 0 is growth.

---

## 2. Two attachment coefficients

Every batch is run twice, at a low and a high α_c. β_sub0 = 7.9408e3/α_c s/m
at −20 °C (Hertz–Knudsen).

| | α_c | β_sub0 [s/m] | kinetic share of resistance | stands for |
|---|---|---|---|---|
| low | 1e-3 | 7.94e6 | 0.60 | lunar regolith |
| high | 1e-2 | 7.94e5 | 0.22 | laboratory experiments |

Talking points

- Low α_c is the lunar case: very cold, near-saturated pores sit at low
  supersaturation, where attachment is nucleation-limited and α_c is small.
  1e-3 is the bottom of the literature band (1e-3 to 1e-1) and the value the
  enceladus campaign uses.
- High α_c is the laboratory case, and it is diffusion-limited: the wall
  distance, not the interface, sets the rate.
- 1e-2 rather than 1e-1 for the high value because `comp_eps.py` returns the
  same ε = 0.8584 µm for every α_c from 1e-4 to 1e-2, so both values reuse the
  existing meshes. At 1e-1 the heat bound binds, ε drops to 0.296 µm and the
  mesh is 3× finer each way, for a predicted velocity change of only ~15 %.
- The high-α_c runs carry `-dtmax 3.35e2` (τ_sub falls to 670 s), so they take
  about 4× the steps of the low-α_c runs.
- Thin-interface offset: the measured interface coefficient on the August
  wedge batch was β_sub0 plus an additive 8.7e5 s/m, i.e. +22 % at 2e-3. It
  is not a fixed percentage: +11 % at 1e-3, +110 % at 1e-2. The curves here
  include it; `--beta-offset 0` redraws them for a solver that removes it
  exactly. The high-α_c batch is the sensitive test of which is right.
- Predicted at θ = 60°, σ∞ = 0: channel 48 nm/day (α_c = 1e-3) against
  95 nm/day (1e-2); wedge ice change over 150 days +14 % against +28 %.
- Figures below are the α_c = 1e-3 set; the `_ac1e-2` files beside them are
  the same panels at the high value. The shapes are identical; only the
  velocity scale and the inner/outer slope ratio change.
- All runs are at −20 °C. Lunar temperatures are far lower; ε and the mesh are
  temperature-dependent, so moving T is a separate step.

---

## 3. Run matrix: four batches, four plots

| Batch | Geometry | Swept | Held fixed | Plot |
|---|---|---|---|---|
| A | Channel | θ = 30° to 150° | σ∞ = 0 | v_n against θ |
| B | Channel | σ∞ = ±3e-5 | θ = 60° | v_n against σ∞ |
| C | Wedge | θ = 30° to 150° | σ∞ = 0 | v_n against θ, both menisci |
| D | Wedge | σ∞ = ±3e-5 | θ = 60° | v_n against σ∞, both menisci |

Talking points

- All runs at −20 °C, temperature pinned on every face, vapour density fixed on
  the left and right walls, wetting walls top and bottom.
- Existing geometries: the 125 µm channel and the 300 µm wedge with wall slopes
  ±0.25 (half-angle α = 14.0°).
- The values in the table are a proposal, not yet agreed. Five angles (30, 60,
  90, 120, 150) and about seven σ∞ values would give 24 runs.
- The σ∞ range is set by the capillary term. d0 = 1.0166e-9 m and
  χ = 8e3 1/m at 60° give d0·χ = 8e-6, so the sweep has to resolve humidity in
  the fifth decimal place. A coarser sweep buries the curvature term.

---

## 4. Channel: velocity follows cos θ

![Channel predictions](fig2_channel_ac1e-3.png)

```
χ = −2·cos θ / H
```

Talking points

- Curvature is uniform along the channel, so θ fixes it outright.
- (a) Batch A should trace a cosine through zero at 90°. A nonzero σ∞ shifts
  the whole curve up or down.
- (b) Batch B should give parallel straight lines of slope 1/(β + K·ℓ),
  crossing zero at σ∞ = −2·d0·cos θ / H.
- Caveat on "velocity is a function of θ alone": Z = K·ℓ depends on the
  meniscus-to-wall distance ℓ, which shrinks or grows as the bridge moves, so
  v_n drifts slowly in time. Integrating gives
  `β·x + K·(ℓ0·x − x²/2) = (σ∞ − d0·χ)·t`.
- Compare runs at matched ℓ, or use the early-time velocity.
- Curves drawn at ℓ = 101 µm, H = 125 µm, α_c = 1e-3.

---

## 5. Wedge: curvature falls as 1/r, giving an ODE

![Wedge curvature and ODE](fig3_wedge_ode_ac1e-3.png)

```
χ_inner = −(cos θ + sin α) / (r·sin α)
χ_outer = −(cos θ − sin α) / (r·sin α)

dr_out/dt = +[σ∞ + d0·(cos θ − sin α)/(r_out·sin α)] / [β + K·r_out·ln(r_R/r_out)]
dr_in/dt  = −[σ∞ + d0·(cos θ + sin α)/(r_in·sin α)]  / [β + K·r_in·ln(r_in/r_L)]
```

Talking points

- r is the meniscus position on the centreline, measured from the virtual apex.
  r_L = 100 µm and r_R = 400 µm are the two fixed-σ∞ walls.
- Derivation: a circular arc centred on the wedge axis at distance c from the
  apex meets both walls at angle θ when c·sin α = R·cos θ. Both forms reduce to
  −1/r and +1/r for the apex-centred 90° band.
- (a) Curvature decays as 1/r for both menisci. The inner meniscus is always
  the more concave of the pair.
- (b) The ODE integrated at σ∞ = 0 for 150 days from the band IC
  (r_in = 200 µm, r_out = 300 µm). At 60° both menisci advance, the inner one
  fastest. At 120° both retreat, the outer one fastest.
- Kinetic limit at σ∞ = 0: r·dr/dt = const, so r² is linear in t.
- With σ∞ ≠ 0 each meniscus has a rest radius
  r* = −d0·(cos θ ∓ sin α)/(σ∞·sin α). It is stable for the outer meniscus and
  unstable for the inner one.

---

## 6. Wedge: the two menisci split about 90°

![Wedge predictions](fig4_wedge_ac1e-3.png)

Talking points

- (a) Batch C: the inner meniscus stalls at θ = 90° + α and the outer at
  90° − α. Between 76° and 104° one face grows while the other sublimates.
- The flat-channel curve sits between the two, crossing at 90°.
- (b) Batch D: each meniscus has its own slope and its own zero crossing at
  σ∞ = d0·χ, and both move as the bridge migrates.
- Velocities are evaluated at the initial band positions, so these are initial
  velocities; later ones follow the ODE.
- v_n here is the centreline value. The contact line moves faster along the
  wall, by the ratio of its apex distance to the centreline's.

---

## 7. Check against the runs we already have

| Reservoir run, σ∞ = 0 | θ | Predicted ice change | Measured ice change |
|---|---|---|---|
| Channel, 90 days | 60° | +6.9 % | +7.5 % |
| Channel, 90 days | 120° | −6.9 % | −7.5 % |
| Wedge, 150 days | 60° | +19.3 % | +21.3 % |
| Wedge, 150 days | 120° | −20.4 % | −21.8 % |

Talking points

- These existing runs used α_c = 2e-3, so the check is made at that value
  (`--alpha 2e-3`). No fitted parameters: β includes the offset measured on
  the August wedge batch; everything else is an input.
- The theory sits 1 to 2 points under the measurement in every case, with the
  right sign and the right channel-to-wedge ratio.
- Predicted ice change uses centreline displacement only: 2·v·t over the
  172 µm bridge for the channel (fixed ℓ), and r_out² − r_in² from the ODE for
  the wedge. Neither includes the meniscus shape correction, so a small
  consistent shortfall is not surprising.
- Measured values are from `../wedge_2026-09-19/README.md`. The channel figure
  is quoted there as about ±7.5 % under the same BCs; I assumed that was the
  90-day `grow90` run. **To confirm.**

---

## 8. Input files and submission

The four batches are built at both α_c values: 44 experiment files and eight
tests files, `batch{A,B,C,D}_*_ac1e-3_tests.txt` and `..._ac1e-2_tests.txt`.

| Batch | Runs per α_c | Experiment files |
|---|---|---|
| A | 5 | `grow90_T-20_theta{30,60,90,120,150}_ac*` |
| B | 6 | `grow90_T-20_theta60_sig{m,p}{1,2,3}e-5_ac*` |
| C | 5 | `wedgeres150_T-20_theta{30,60,90,120,150}_ac*` |
| D | 6 | `wedgeres150_T-20_theta60_sig{m,p}{1,2,3}e-5_ac*` |

- σ∞ is set through `-rhovfix_lo/hi`, which are fractions of ρ_vs(temp0), so
  σ∞ = rhovfix − 1. `-humidity` is set to the same value so the vapour IC
  starts on the reservoir value.
- The σ∞ = 0 points of B and D are the θ = 60° runs of A and C.
- Nothing here repeats an earlier run: the September batches used α_c = 2e-3,
  and the five-angle `grow90` batch of 2026-09-14 was thermally throttled
  (`-flag_BC_Tfix 0`, 0.003 % ice change) besides.
- Submit in stages, cheapest first, and check A before the rest:

```bash
./scripts/HPC/submit_batch.sh --tag velA_channel_theta_ac1e-3 \
    --tests-file studies/contact_angle/velocity_plan_2026-10-01/batchA_channel_theta_ac1e-3_tests.txt
```

- Measurement: `postprocess/meniscus_velocity.py` for the channel,
  `postprocess/wedge_gt_velocity.py` for the wedge, and
  `postprocess/gt_balance.py` for the σ = d0·κ + β·v_n balance on any run.
  Background in `docs/curvature_driven_growth.md` and
  `docs/enceladus_carryover.md`.

---

### Wall time and time limits

Estimates, not measurements. They rest on one data point (the September
channel runs: 3.2 s/step at 8 ranks) and on the enceladus scaling fit for how
the step slows with fewer ranks. Treat them as ±50 %.

The step count is set by `-dtmax`, not by the interface-CFL limiter, which
rarely binds on these slow interfaces. The experiment files set
`-dtmax = 2·τ_sub` for their own α_c (1.17e4 s at 1e-3, 1.34e3 s at 1e-2),
the value the dtmax ladder below found indistinguishable from 0.8·τ_sub.

| Runs | Steps | 200k DoF/core | 50k DoF/core |
|---|---|---|---|
| channel, α_c = 1e-3 (A, B) | 750 | ~3 h (2 ranks) | ~1 h (8 ranks) |
| wedge, α_c = 1e-3 (C, D) | 1,200 | ~4 h (3 ranks) | ~1.5 h (10 ranks) |
| channel, α_c = 1e-2 (A, B) | 5,900 | ~23 h | ~6 h |
| wedge, α_c = 1e-2 (C, D) | 9,800 | ~34 h | ~11 h |

- `scripts/lib/alloc.sh` targets 200k DoF/core. On meshes this small that
  costs about the same core-hours as 50k but roughly 4× the wall time;
  `TARGET_DOFS_PER_CORE=50000` on the submit line trades it back.
- α_c = 1e-3: `-- --time=0-06:00:00` (channel) and `0-08:00:00` (wedge) at
  200k.
- α_c = 1e-2 at 200k: the channel runs sit at the 24 h mark and the wedge
  runs exceed it, so use the 50k target or restart legs for those.

---

### Is `-dtmax` the right ceiling?

Two different limits act on the timestep, and they guard different things.

- **Interface-CFL limiter** (`-dtCFL`, `-dtCFL_dphimax 0.2`): dynamic. After
  each step it measures the largest pointwise change in φ and rejects the
  step if any point moved more than 0.2. It protects fast events (a lens
  vanishing, a neck snapping) and is silent when the interface is slow.
- **`-dtmax`**: static. The adaptive stepper grows dt whenever Newton
  converges quickly, and with a slow interface nothing else stops it. The
  ceiling is what bounds temporal accuracy in that regime.

Some ceiling is needed, but the value is not derived. The τ_sub scaling
dates from the summer, before the bounds rollback, the wall-term clamp and
the CFL limiter existed, when a large dt simply broke the solve. The only
measurement since is `dtmax_study.sh`: the contact angle moved 0.003° across
0.05 to 0.815·τ_sub. Nothing has been tested above that, and these
interfaces move about one element per 180 steps at 0.8·τ_sub.

`dtmax_ladder_tests.txt` tests it directly: the θ = 60°, α_c = 1e-3 channel
run at 2, 5 and 10·τ_sub (about 660, 270 and 130 steps against 1,700).
Compare the meniscus velocity and `solver_evo.dat` against the batch A run.

A local version runs the same test further and faster: six 30-day runs at
0.8, 2, 5, 10, 20 and 40·τ_sub (`dtlad30_T-20_theta60_ac1e-3_dt*tau`), scored
by `analyze_dtmax_ladder.py` on the per-step ice area in `SSA_evo.dat` and the
iteration counts in `solver_evo.dat`. 40·τ_sub is the largest the six-snapshot
output cadence allows.

Three local cases, each a six-rung ladder (0.8 to 40·τ_sub), all scored by
`analyze_dtmax_ladder.py`:

| Case | Geometry | Experiment files | Scored on |
|---|---|---|---|
| channel growth | `channel_2D_H31um_eps0.86um` | `dtlad20_T-20_theta60_ac1e-3_dt*tau` | ice area |
| sintering pair | `sinterpair_2D_L120um_eps0.86um` | `dtlad60_T-20_sealed_ac1e-3_dt*tau` | interface length, event time |
| ripening pair | `ripenpair_2D_L132um_eps0.86um` | `dtlad60_T-20_sealed_ac1e-3_dt*tau` | interface length, event time |

The grain pairs (15 and 30 µm radius, tangent or 12 µm apart, sealed box) end
in the small grain disappearing. That is the regime the CFL limiter was built
for, and the one where a large step is known to cost accuracy: on 2026-07-10 a
12× larger `-dtmax` stayed stable and ripple-free but put extinction 32 % late.

#### Result: small channel (batch_2026-10-01__14.25.24_dtmax_small)

All six rungs completed, no rejections, φ within bounds. Ice grew 29 % in
20 days at the reference.

| dtmax/τ_sub | steps | growth error | half-growth time | CFL caps |
|---|---|---|---|---|
| 0.8 | 428 | reference | 10.07 d | 0 |
| 2 | 210 | 0.06 % | +0.1 % | 0 |
| 5 | 124 | 0.6 % | +1.0 % | 1 |
| 10 | 97 | 2.8 % | +3.7 % | 2 |
| 20 | 85 | 10 % | +10.5 % | 5 |
| 40 | 81 | 22 % | +16.0 % | 9 |

- Stability is not the limit: every rung ran clean. Accuracy is. The error
  grows about quadratically in dt and the larger step always under-predicts
  growth.
- 2·τ_sub costs 0.06 % and halves the steps; 5·τ_sub costs 0.6 % for 3.5×
  fewer. Past 10·τ_sub the error is no longer small against the effects the
  study measures.
- The CFL limiter starts binding near dt ≈ 4e4 s (7·τ_sub) on this geometry,
  so the 20 and 40 rungs are partly limiter-controlled; about 60 of every
  run's steps are the ramp up from the 1e-4 s first step.
- **Every step takes exactly one Newton iteration, at every rung.** The φ
  residual is ~1e-17 in SI units, far below `-snes_atol 1e-6`, so the block
  is declared converged on entry; after the single iteration its residual has
  fallen by only 0.54×. The scheme is in effect linearised once per step.
  This is the known, deferred tolerance issue noted in `lunar_main.c`
  (per-field tolerances matched to each field's scale), and it is the
  concrete case for non-dimensionalising: not conditioning (5–7 Krylov
  iterations per solve is already cheap) but a convergence test that means
  something. Whether a properly converged Newton solve would hold accuracy
  at larger dt is untested.

#### Result: sintering pair (batch_2026-10-01__14.52.38_dtmax_sinter)

All six rungs completed, no rejections. The small grain was absorbed: a neck
forms within about two days, the pair becomes a pear, and the small lobe is
consumed from its outer side until nothing is left at its old centre between
days 48 and 54.

| dtmax/τ_sub | steps | max trajectory error | absorption time | CFL caps |
|---|---|---|---|---|
| 0.8 | 1,167 | reference | 43.0 d | 0 |
| 2 | 505 | 0.6 % | 0 % | 0 |
| 5 | 243 | 1.5 % | +0.1 % | 3 |
| 10 | 158 | 2.5 % | +2.1 % | 7 |
| 20 | 118 | 3.6 % | +11.5 % | 13 |
| 40 | 103 | 3.6 % | +24 % | 24 |

- "Absorption time" is when the ice height at the small grain's original
  centre falls to 10 µm, interpolated between snapshots six days apart, so it
  is good to about ±2 %. Trajectory error is on the per-step interface length.
- Same verdict as the channel: clean to 5·τ_sub, a few percent at 10, and
  wrong beyond. A large step makes the event LATE, as in the July stress test.
- The script's `t_half` is not useful here: half the interface-length change
  happens in the first two days of neck formation, while every rung is still
  ramping its timestep.
- Neck width, measured as the ice height at the original contact plane
  (13.5 µm at t = 0, 50 µm at 60 d), lags the reference by at most 0.01 µm at
  2·τ_sub, 0.11 µm at 5, 0.39 µm at 10 and about 1.0 µm at 20 and 40
  (0.0, 0.3, 1.1 and 2.9 % of the rise). Always behind, never ahead. A true
  neck (a waist in the outline) exists only before day 6, i.e. before the
  first snapshot, so the early neck-growth law is not resolved by these runs.
- One Newton iteration per step again. Krylov iterations per solve are 22–27
  on this 88k-DoF mesh against 5–7 on the 24k-DoF channel, both on one rank.

#### Result: ripening pair (batch_2026-10-01__15.49.53_dtmax_ripen)

All six rungs completed, no rejections, φ undershoot at most 6e-5. The small
grain shrinks from 15 µm equivalent radius to 6.6 µm by day 24 and is gone by
day 27 in the reference; the large one grows from 30.0 to 33.5 µm.

| dtmax/τ_sub | steps | extinction time | vs reference | CFL caps |
|---|---|---|---|---|
| 0.8 | 1,167 | 26.56 d | reference | 0 |
| 2 | 506 | 26.58 d | +0.1 % | 2 |
| 5 | 247 | 26.80 d | +0.9 % | 4 |
| 10 | 165 | 27.53 d | +3.6 % | 7 |
| 20 | 127 | 29.85 d | +12.4 % | 11 |
| 40 | 113 | 34.65 d | +30.5 % | 15 |

- Extinction time is when the per-step interface length has completed 97 % of
  its total drop; it is resolved to one step.
- The limiter throttles dt to about 1e4 s around the extinction itself at
  every rung from 10 upward, so the lateness is accumulated during the slow
  shrinkage before it, not at the event.
- The interface-length trajectory at 20 and 40·τ_sub departs from the
  reference well before the grains differ (radii at day 6 agree to 0.01 µm),
  so that integral also responds to profile shape under large steps. Radii
  and extinction time are the trustworthy measures here.

**Adopted 2026-10-02: `-dtmax = 2·τ_sub` everywhere** (the velocity-study
inputs, the geometry-file defaults, and the enceladus production inputs and
generator). The solver's startup warning now fires above 5·τ_sub.

#### All three cases together

| dtmax/τ_sub | channel growth | sintering absorption | ripening extinction | steps saved |
|---|---|---|---|---|
| 2 | 0.06 % | 0 % | +0.1 % | 2.0–2.3× |
| 5 | 0.6 % | +0.1 % | +0.9 % | 3.5–4.8× |
| 10 | 2.8 % | +2.1 % | +3.6 % | 4.4–7.4× |
| 20 | 10 % | +11.5 % | +12.4 % | 5–10× |
| 40 | 22 % | +24 % | +30.5 % | 5–11× |

Three different mechanisms give the same curve: error roughly quadratic in
dt, always in the direction of "too slow", never unstable. 2·τ_sub is free.
5·τ_sub costs under 1 %. 10 costs 2–4 %. The returns also flatten: the ~60–100
steps of ramp from the 1e-4 s first step are a fixed cost, so going from 5 to
40 saves little.

### Solver diagnostics

Every run now writes `solver_evo.dat` beside `SSA_evo.dat`:
`step t dt newton_its krylov_its krylov_per_newton rejections`, one row per
step. `krylov_per_newton` is the conditioning measure (how hard each linear
solve is); `rejections` counts step attempts thrown away. It reads counters
the time integrator already keeps, so it adds no solver work.

---

## 9. Where this is going: a mm-scale pore channel

![Pore channel with lenses and wall-adhered ice](fig5_pore_channel.png)

Talking points

- Schematic only: the wall shapes and ice placement are illustrative, not a
  simulation or a planned geometry file.
- Boundary conditions: the left end is open to a low vapour saturation
  (σ∞ < 0); the right end is a dead pore, i.e. no vapour flux.
- Two kinds of ice. Lenses span both walls and sit at the throats.
  Wall-adhered ice touches one wall only and sits in the troughs, not on the
  peaks.
- Every lens meniscus is locally the wedge problem: the walls diverge or
  converge at some local half-angle, so batches C and D are the building
  block for this picture.
- Wall-adhered ice is the sessile case. It has one meniscus and its curvature
  is set by θ and the trough shape; it is not covered by the four batches.

---

## 10. Open points before submitting

- **Velocity is not set by θ alone.** Diffusion to the wall is as large as the
  kinetic term, so v_n also depends on the meniscus-to-wall distance. Compare
  at matched distance or take early-time velocities.
- **Dynamic contact angle.** A moving interface sat about 1° off Young's angle
  at 60° in the reservoir runs. Plot against the measured angle as well as the
  prescribed one.
- **Radial-diffusion approximation.** The wedge resistance assumes vapour
  diffuses radially from an arc-shaped wall. The real walls are straight, and
  the meniscus is not an arc about the apex unless θ = 90°.
- **Sweep ranges.** Which angles, which σ∞ values, and which fixed θ for the
  saturation batches.
- **Where to measure.** Centreline (as `meniscus_velocity.py` and
  `wedge_gt_velocity.py` do) or contact line; they differ in the wedge.
