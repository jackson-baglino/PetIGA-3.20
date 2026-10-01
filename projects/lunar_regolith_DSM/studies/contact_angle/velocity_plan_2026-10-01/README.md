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

## 9. Where this is going: a mm-scale pore channel

![Pore channel with lenses and wall-adhered ice](fig5_pore_channel.png)

Talking points

- Schematic only: the wall shapes and ice placement are illustrative, not a
  simulation or a planned geometry file.
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
