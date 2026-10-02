# Historical two-sphere sintering experiments: Kingery (1960) and Thomas et al. (1994)

Molaro et al. (2019) compare their Swinkels–Ashby model against three historical
neck-growth datasets in their **Fig. 8**. We have two of the source papers. This
folder covers those two: where each experiment is described, the data, and
what it would take for our model to reproduce it. The third dataset, Hobbs &
Mason (1964), is Fig. 8(a). It is left out because the paper is paywalled.

![Kingery and Thomas neck data](figures/historical_neck_data.png)

*Relative neck size x/a against time, as digitised from Molaro et al. (2019)
Fig. 8. The dashed line is the resolution floor of the mesh proposed below.
Each mesh is sized so the floor sits at 0.9× the first measured point.*

Regenerate (from `enceladus_DSM/`):

```bash
python studies/historical_sintering/analysis/plot_data.py   # figure
python studies/historical_sintering/analysis/size_runs.py   # data/run_sizing.csv
```

---

## 1. Where the experiments are in each paper

### Molaro et al. (2019), JGR Planets 124, 243–276

| what | where |
|---|---|
| Discussion of the comparison | §3.2 *Behavior of Water Ice*, pp. 253–254 |
| The figure | **Fig. 8**, p. 256: (a) Hobbs & Mason, −5 °C, a = 288 µm; **(b) Kingery, −18 °C, a = 110 µm**; **(c) Thomas et al., −20 °C, a = 120 µm** |
| Definition of "relative neck size" | §2.3, p. 249: neck **radius** / grain **radius** (x/a, their Table 1 and Fig. 2) |

Their verdict (p. 254): the SA81 model is much slower than Kingery and much
faster than Thomas, "in spite of the fact that both studies used comparable
grain sizes and temperatures". Both were run in a cold room at atmospheric
pressure. Molaro suggests uncontrolled humidity as the likely cause. One slip:
§3.2 cites "Kingery & Berg, 1955" for this dataset, but the Fig. 8 caption and
the data are Kingery (1960).

### Kingery (1960), *Regelation, surface diffusion, and ice sintering*, J. Appl. Phys. 31, 833–838

Page numbers are the journal's. In `Kingery1960.pdf`, PDF page 1 is an AIP cover sheet, so journal p. 833 is PDF p. 2.

| what | where |
|---|---|
| Method | **§IV *Experimental Observations*, p. 835, right column.** Ice spheres were made by spraying water into liquid oxygen; radii 0.1–3 mm were selected. Two or more spheres were placed in a line on a glass slide under a microscope, brought together lightly, and measured with a micrometer ocular. The whole apparatus sat in a constant-temperature room in ordinary air, with no enclosure. Most runs were at −17.8 °C; the full range was −2.2 to −25.1 °C. |
| **The data Molaro used** | **Fig. 3(c), p. 836**: neck growth at −17.8 °C, the **r = 0.011 cm (110 µm)** curve. Axes are x/r against time (0.1–100 min, log–log). |
| Growth exponents | **Table II, p. 836**: the r = 0.011 cm sphere has inverse slope **n = 5.5**; the −17.8 °C mean is 6.9. Eq. (4), p. 835, gives x/r ~ t^(1/6.1–7.1). |
| Neck geometry | Fig. 1, p. 834: x is the neck **radius** and r is the sphere radius. This is the same ratio as Molaro's x/a. |
| Mechanism claim | p. 836: n ≈ 7 and no centre approach, so Kingery concludes **surface diffusion**. Fig. 4 (p. 837) shows time to x/r = 0.1 scaling as r⁴. Fig. 5 is an Arrhenius plot with an activation energy of 27.5 kcal/mol. |

### Thomas, Ratke & Kochan (1994), *Crushing strength of porous ice-mineral bodies*, Adv. Space Res. 14(12), (12)207–(12)216

Molaro's reference list gives the authors as Thomas, Kochan & Ratke; the paper itself prints Thomas, Ratke and Kochan.

| what | where |
|---|---|
| **The two-sphere experiment and data** | **p. (12)210 (PDF p. 4).** Two ice particles under a microscope in a cold lab at −20 °C. Fig. 2 shows photos after 3 min, 1 day and 1 week. **Fig. 3** plots x/r against time (10³–10⁶ s, log–log) at **T = 253 K, r ≈ 120 µm**, with fit **x/r ~ t^(1/3)**. The authors conclude sublimation–recondensation dominates at −20 °C. |
| Theory and its validity | p. (12)209: exponent table (n = 3 vapour, 4–5 surface diffusion, 6 grain boundary). The power laws hold only for x/r < 0.5. |
| Particle preparation | p. (12)211: water sprayed into liquid nitrogen. This is stated for the agglomerate experiments; the two-sphere grains are not described separately. |
| Restated | Summary, p. (12)215. |

The rest of the paper (strength measurements, Figs. 4–7; surface diffusion at 30 K, Fig. 8) is not part of Molaro's comparison.

---

## 2. The data

| file | source | points | window |
|---|---|---|---|
| `data/kingery1960_T-17.8_r110um.csv` | Molaro Fig. 8(b) | 4 | x/a 0.249 → 0.333 over 1.6 → 8.3 min |
| `data/thomas1994_T-20_r120um.csv` | Molaro Fig. 8(c) | 11 | x/a 0.099 → 0.443 over 3.3 → 190 h (~8 d) |

- **Kingery points.** `Molaro2019_Fig8b.csv` was empty, so I digitised the four points from `Moalro2019_Fig8b.png`. I located the marker centroids and mapped them through the axis box and tick pixels. The calibration is in the CSV header. The values agree with an eye-read to within 0.003.
- **Thomas points.** Copied from your `Molaro2019_Fig8c.csv`, with columns swapped to (time, neck).
- **Second-hand data.** Both series are Molaro's re-plot, not the original figures. Kingery's Fig. 3(c) and Thomas's Fig. 3 are in the PDFs. Digitising them directly would remove one layer of hand extraction before any manuscript use.
- **Clock zero.** Neither experiment knows its contact time. Thomas's axis starts at their first image (Fig. 2a was taken at 3 min). Compare by **anchoring our clock at the first measured neck size**, as in Molaro Fig. 12 and `fit_neck_growth.py --anchor-neck`. Do not compare raw clocks.

---

## 3. The simulations we would need

All numbers come from `analysis/size_runs.py` (output in `data/run_sizing.csv`). The common setup is the one validated in `studies/sinter_exponent`:

- Two equal spheres, axisymmetric (`-axisym 1`).
- Exact tangency (`-ic_grain_union 1`), so there is no chosen initial neck.
- α_c = 1e-3 held constant: `-alpha_pointwise 1 -alpha_model 0 -alpha_c0 1e-3`.
- Physical d₀, `--Dchannel ice`.
- `-dtmax = 2·tau_sub`.
- Sealed box at h = 1.00 as the baseline.
- Grain centres at x = 1.2R and 3.2R, with Lx = 4.4R (0.2R padding at each end) and Ly = 1.2R.

### How eps is chosen: resolve the neck just before the first data point

A neck is measurable only above `x/a ≥ sqrt(12·eps/a)`; below that, the fillet is thinner than the interface. That floor goes at 0.9× the first measured x/a:

    eps = a · (0.9 · u_first)² / 12

At α_c = 1e-3, the K&P bounds from `comp_eps.py` are microns, one to two decades looser than this. **The neck floor is what binds.** Mesh spacing is h = eps/√2.

| | **Kingery** | **Thomas (strict)** | Thomas (relaxed) |
|---|---:|---:|---:|
| T | −17.8 °C | −20 °C | −20 °C |
| grain radius a | 110 µm | 120 µm | 120 µm |
| resolve from | point 1, x/a = 0.249 | point 1, x/a = 0.099 | point 3, x/a = 0.181 |
| floor x/a | 0.224 | 0.089 | 0.163 |
| **eps** | **0.459 µm** | **0.080 µm** | 0.266 µm |
| h | 324 nm | 56 nm | 188 nm |
| domain Lx × Ly | 484 × 132 µm | 528 × 144 µm | 528 × 144 µm |
| **Nx × Ny** | **1493 × 408** | **9368 × 2555** | 2811 × 767 |
| nodes / DoF | 0.61 M / 1.8 M | 23.9 M / 71.8 M | 2.16 M / 6.5 M |
| ranks (200k DoF/rank) | 10 | 360 | 33 |
| tau_sub → dtmax | 1351 s → 2702 s | 49.8 s → 99.5 s | 555 s → 1110 s |
| d₀ (`-d0_sub0`) | 1.0064e-9 m | 1.0152e-9 m | 1.0152e-9 m |
| **experiment duration** | **~8 min** (1.6 → 8.3 min) | **~8 days** (3.3 → 190 h) | ~7 days (20 → 190 h) |
| model time to first point (est.) | 52 h | 2.5 h | 23 h |
| model time to last point (est.) | 155 h | 659 h (27 d) | 659 h |
| **model slower than experiment, over the data window** | **~920×** | **~3.5×** | ~3.7× |
| **t_final** (1.25 × last point) | **194 h (6.98e5 s)** | **824 h (2.97e6 s)** | 824 h |
| **steps** | **~360** | **~30,000** | ~2,800 |
| **wall time** (8 / 15 / 23 s per step) | **0.8 / 1.5 / 2.3 h** | **67 / 125 / 191 h** | 6 / 12 / 18 h |
| **cost** (Tier 1, $0.012/core-h) | **< $1** | **$290 / $540 / $830** | $2 / $5 / $7 |

How the estimates were made:

- **Model time.** The model trajectory t(u) is calibrated on the `mesh_pair` fine arm: α_c = 1e-3, −20 °C, tangent start, u = 0.194 → 0.303 in 15 → 79 h, so t ~ u^3.73. It is scaled to each case by the kinetic-limit law t ∝ β_sub·a²/d₀. Applying this to Thomas means extrapolating from u = 0.19–0.30 out to 0.09–0.44. **Treat the model times as ±2×.** Near u = 0.44 the exact fillet geometry slows growth further, so the Thomas end time is more likely an underestimate.
- **Seconds per step.** The low end is the Molaro axisymmetric runs (~7 s/step at ~100k DoF/rank). The high end is the phase-field strong-scaling fit, t ≈ 15.6 s + 940/P core-s per step, scaled by DoF (`project_hpc_cost_model`). Startup and CFL-limited early steps add ~100 to every step count.
- **Steps do not depend on α_c.** In the kinetic limit, both the neck time and tau_sub scale as 1/α_c, so the number of steps is ≈ (a/eps)²·u^3.7 regardless of α_c. A different α_c changes `t_final`, not the cost.

### What the table says

**Kingery is essentially free but is not a fair test of this model.**

- 10 ranks, about 360 steps, under $1.
- The model should take about 900× longer than the experiment to cross the same x/a range.
- Kingery's own exponent (n = 5.5 for this sphere, 6.9 on average) says surface diffusion. Our model has no surface diffusion, so it will not match the shape either.
- It belongs in the manuscript as the **contrast case**. That fits `project_molaro_validation`, which identifies surface diffusion as the missing mechanism.
- Since it is cheap, run both boundary arms: sealed, and Dirichlet at h = 1.00.

**Thomas is the comparison that matters, and it is the expensive one.**

- Its fitted t^(1/3) is the vapour route this model contains.
- The model is predicted to be only ~3.5× slow there, against ~150× for Molaro's own −20 °C pair.
- Resolving the first point (x/a = 0.099) needs eps = 0.080 µm: 24 M nodes, ~30k steps, roughly 3–8 days of wall time on 360 ranks. Several 24 h restart legs would be needed, at $300–800.
- The cost comes from eps scaling as u_first². Starting the comparison at **point 3 (x/a = 0.181)** instead gives eps 3.3× larger, 11× fewer nodes and steps, and **~$5 in a single 12 h job**. It still covers 9 of the 11 points and a full decade of time (20 → 190 h).
- A mirror half-domain (contact plane on the x = 0 zero-flux wall) would halve the node count again. This is untested: `neck_width.py` and the Dirichlet arm would both need checking first.

**Recommendation.** Run Kingery (both arms) and Thomas-relaxed as one cheap staged batch. Use Thomas-relaxed to confirm the predicted ~3.5× rate and the 1/3 slope. Fund the strict Thomas run only if the first two points (x/a < 0.18) are needed for the paper.

### Inputs to write (none created yet)

Geometry (`inputs/geometry/<family>/`). Grain centres: `-ice_grain_cx 1.2R,3.2R`, `-ice_grain_cy 0,0`, `-ice_grain_R/ax/ay R,R`.

| | Kingery | Thomas strict | Thomas relaxed |
|---|---|---|---|
| `-Lx`, `-Ly` | 4.84e-4, 1.32e-4 | 5.28e-4, 1.44e-4 | 5.28e-4, 1.44e-4 |
| `-ice_grain_cx` | 1.32e-4, 3.52e-4 | 1.44e-4, 3.84e-4 | 1.44e-4, 3.84e-4 |
| `-Nx`, `-Ny` | 1493, 408 | 9368, 2555 | 2811, 767 |
| `-eps` | 4.585e-7 | 7.97e-8 | 2.657e-7 |
| `-eps_valid_temp` | −17.8 | −20 | −20 |

Experiment (`inputs/experiment/<family>/`). All three use `-alpha_pointwise 1 -alpha_model 0 -alpha_c0 1e-3`, `-humidity 1.00`, `-dtCFL 1 -dtCFL_dphimax 0.2`, `-t_out_log 80`, and `-t_out_first 1`.

| | Kingery | Thomas strict | Thomas relaxed |
|---|---|---|---|
| `-temp` | −17.8 | −20 | −20 |
| `-d0_sub0` | 1.0064e-9 | 1.0152e-9 | 1.0152e-9 |
| `-dtmax` | 2702 | 99.5 | 1110 |
| `-t_final` | 6.98e5 | 2.97e6 | 2.97e6 |

The temperature guard is ±1 °C. Molaro labels Kingery "−18 °C", but use Kingery's −17.8 °C and regenerate the mesh if that changes.

### Analysis

```bash
python postprocess/neck_width.py <run_dir> --axisym          # --axisym is REQUIRED
python postprocess/fit_neck_growth.py <run>/neck_width.csv <data csv> \
    --anchor-neck-rn <x/a of first compared point> --demmenie
#   Kingery 0.2485 · Thomas strict 0.0992 · Thomas relaxed 0.1811
```

`neck_width.py --axisym` returns the neck **width** (2x); divide by 2a for x/a.
`fit_neck_growth.py` currently reads data series in the Molaro layout
(`time_min,neck_size_um,...,lg_diam_um,sm_diam_um`). Before the fit can run,
the two CSVs here need converting to that layout (neck = 2·x/a·a, both
diameters = 2a), or the loader needs a `time_h,relative_neck` branch.
