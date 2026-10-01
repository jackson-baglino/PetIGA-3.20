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
- At ℓ = 100 µm, Z = 5.8e6 s/m against β_eff = 4.8e6 s/m. The runs are in a
  mixed kinetic/diffusive regime, not kinetics-limited.
- β_eff = 1.22 × the requested β_sub0. This is the thin-interface calibration
  offset measured on the 2026-08-07 wedge batch, and it is the value used in
  every curve here.
- Sign conventions: θ is measured through the ice; χ > 0 where the ice is
  convex into the vapour; v_n > 0 is growth.

---

## 2. Run matrix: four batches, four plots

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

## 3. Channel: velocity follows cos θ

![Channel predictions](fig2_channel.png)

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
- Curves drawn at ℓ = 101 µm, H = 125 µm.

---

## 4. Wedge: curvature falls as 1/r, giving an ODE

![Wedge curvature and ODE](fig3_wedge_ode.png)

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

## 5. Wedge: the two menisci split about 90°

![Wedge predictions](fig4_wedge.png)

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

## 6. Check against the runs we already have

| Reservoir run, σ∞ = 0 | θ | Predicted ice change | Measured ice change |
|---|---|---|---|
| Channel, 90 days | 60° | +6.9 % | +7.5 % |
| Channel, 90 days | 120° | −6.9 % | −7.5 % |
| Wedge, 150 days | 60° | +19.3 % | +21.3 % |
| Wedge, 150 days | 120° | −20.4 % | −21.8 % |

Talking points

- No fitted parameters. β is the value measured on the August wedge batch;
  everything else is an input.
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

## 7. Input files and submission

The four batches are built. One tests file per batch, in this folder:

| Batch | Tests file | Runs | New input files |
|---|---|---|---|
| A | `batchA_channel_theta_tests.txt` | 5 | none (`grow90_T-20_theta*`) |
| B | `batchB_channel_sigma_tests.txt` | 6 | `grow90_T-20_theta60_sig{m,p}{1,2,3}e-5` |
| C | `batchC_wedge_theta_tests.txt` | 5 | `wedgeres150_T-20_theta{30,90,150}` |
| D | `batchD_wedge_sigma_tests.txt` | 6 | `wedgeres150_T-20_theta60_sig{m,p}{1,2,3}e-5` |

- σ∞ is set through `-rhovfix_lo/hi`, which are fractions of ρ_vs(temp0), so
  σ∞ = rhovfix − 1. `-humidity` is set to the same value so the vapour IC
  starts on the reservoir value.
- The σ∞ = 0 points of B and D are the θ = 60° runs of A and C.
- Batch A is not a repeat of the 2026-09-14 five-angle `grow90` batch. That
  one ran with `-flag_BC_Tfix 0` and was thermally throttled (0.003 % ice
  change). Only θ = 60° and 120° have been run with the thermal bath since.
- Submit in stages, cheapest first, and check A before the rest:

```bash
./scripts/HPC/submit_batch.sh --tag velA_channel_theta \
    --tests-file studies/contact_angle/velocity_plan_2026-10-01/batchA_channel_theta_tests.txt
```

- Measurement: `postprocess/meniscus_velocity.py` for the channel,
  `postprocess/wedge_gt_velocity.py` for the wedge, and
  `postprocess/gt_balance.py` for the σ = d0·κ + β·v_n balance on any run.
  Background in `docs/curvature_driven_growth.md` and
  `docs/enceladus_carryover.md`.

---

## 8. Open points before submitting

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
