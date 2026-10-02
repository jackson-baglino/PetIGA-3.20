# Thomas et al. (1994) two-sphere sintering: data and replication plan

Molaro et al. (2019) compare their Swinkels–Ashby model against three historical
neck-growth datasets in their **Fig. 8**. This folder covers one of them, Thomas,
Ratke & Kochan (1994), Fig. 8(c). It is the one our model should match best:
- its fitted growth law is t^(1/3), the sublimation–condensation route this model contains;
- it runs for a week, so it spans a decade and a half of time.

Kingery (1960) and Hobbs & Mason (1964) were considered and dropped
(2026-10-02). Kingery's n ≈ 7 is a surface-diffusion signature our model
cannot produce, and Hobbs & Mason is paywalled. The Kingery sizing is in
commit `0dc51f16`.

![Thomas neck data](figures/thomas_neck_data.png)

*Thomas Fig. 3, digitised from Molaro et al. (2019) Fig. 8(c).*
- *(a) Relative neck size x/a, as published.*
- *(b) Absolute neck width 2x = 2a·(x/a) with a = 120 µm. This is the quantity `neck_width.py --axisym` reports and our Molaro figures show. Thomas gives only "r ≈ 120 µm", so the absolute scale carries that radius uncertainty; a 10 % error in a shifts (b) by 10 % without changing the exponent.*
- *Filled points (3–11) are the ones the proposed run resolves; open points (1–2) fall below its mesh floor (dashed). The line is the `d_fixed` fit over points 3–11, drawn only over that range.*

Regenerate (from `enceladus_DSM/`):

```bash
python studies/historical_sintering/analysis/fit_power_law.py   # data/thomas_powerlaw_fits.csv
python studies/historical_sintering/analysis/plot_data.py       # figures/thomas_neck_data.*
python studies/historical_sintering/analysis/size_runs.py       # data/run_sizing.csv
```

---

## 1. Where the experiment is described

**Molaro et al. (2019), JGR Planets 124, 243–276**

| what | where |
|---|---|
| Discussion | §3.2 *Behavior of Water Ice*, pp. 253–254. SA81 predicts growth "significantly faster" than Thomas, despite a grain size and temperature comparable to Kingery's. Molaro suggests humidity as the cause. |
| Figure | **Fig. 8(c)**, p. 256: T = −20 °C, a = 120 µm. Axes are time [h] against relative neck size. |
| Definition of relative neck size | §2.3, p. 249: neck **radius** / grain **radius**, x/a |

**Thomas, Ratke & Kochan (1994), *Crushing strength of porous ice-mineral bodies*, Adv. Space Res. 14(12), (12)207–(12)216**

Molaro's reference list gives the authors as Thomas, Kochan & Ratke; the paper prints Thomas, Ratke and Kochan.

| what | where |
|---|---|
| **The two-sphere experiment** | **p. (12)210 (PDF p. 4).** Two ice particles under a microscope in a cold lab at −20 °C. **Fig. 2** shows photos after 3 min, 1 day and 1 week. **Fig. 3** plots x/r against time (10³–10⁶ s, log–log) at T = 253 K, r ≈ 120 µm, with fit **x/r ~ t^(1/3)**. The authors conclude sublimation–recondensation dominates. |
| Theory and validity | p. (12)209: n = 3 for vapour, 4–5 for surface diffusion, 6 for grain boundary. The power laws hold only for x/r < 0.5. |
| Particle preparation | p. (12)211: water sprayed into liquid nitrogen, which was then evaporated at −20 °C. Stated for the agglomerate experiments. |
| Restated | Summary, p. (12)215. |

---

## 2. The data and its power-law fit

`data/thomas1994_T-20_r120um.csv` holds 11 points, x/a 0.099 → 0.443 over 3.3 → 190 h (~8 days). It is your digitisation of Molaro Fig. 8(c), so the data are second-hand. Thomas's own Fig. 3 could be digitised directly before any manuscript use.

The data were fitted to x/a = C·t^a (`analysis/fit_power_law.py`, 95 % CIs, results in `data/thomas_powerlaw_fits.csv`):

| window | form | a | ± 95 % | t₀ [h] | R² |
|---|---|---:|---:|---:|---:|
| all 11 points | `d_fixed`, C·t^a | **0.356** | 0.043 | — | 0.975 |
| all 11 points | `d_free`, C·(t+t₀)^a | 0.386 | 0.082 | 1.5 ± 5.8 | 0.984 |
| points 3–11 (simulated window) | `d_fixed` | **0.403** | 0.044 | — | 0.985 |
| points 3–11 | `d_free` | **0.329** | 0.096 | −9 ± 12 | 0.986 |

- **Over the full record the data agree with Thomas's t^(1/3).** a = 0.36 ± 0.04 contains 1/3.
- **Over the simulated window, a depends on the protocol.** With the clock as reported, a = 0.40 ± 0.04. Freeing the contact time gives 0.33 ± 0.10.
  - Neither experiment knows its contact time, so `d_free` is the honest number. It also has the wide CI.
  - This is the same protocol trap documented in `studies/sinter_exponent`. **Always fit the model and the data with the same form over the same x/a window.**
- **Point-to-point local slopes are unusable on this data.** They swing from −0.05 to 0.80, because digitisation scatter dominates over these short log-intervals. Use window fits only.
- **What the model should give.**
  - The kinetic-limited ideal local slope sags from 0.32 to 0.30 over x/a = 0.15–0.30, and lower at 0.44 (`sinter_exponent/PLAN.md`).
  - The resolved `mesh_pair` arm at −20 °C gave 0.283.
  - The data's 0.33–0.40 window therefore sits **at or above** the model's ideal. The model may come in slightly low on the exponent, but not by a large factor.

---

## 3. Vapour saturation: what Thomas offers

For Molaro's grains, the 3–4 % shrinkage calibrated the chamber humidity (h ≈ 0.998). **Thomas reports nothing equivalent.**

- **No humidity, enclosure or grain-size data.** The paper states no humidity, vapour pressure or enclosure for the two-sphere run, only "a cold lab at −20 °C". It reports no grain radii over time.
- **The Fig. 2 photos can't resolve shrinkage.** They are 3 min, 1 day and 1 week; I extracted them from the PDF. There is no scale bar, and the crop changes between panels. By eye, the lower grain spans about 80 % of the frame width in all three. So no shrinkage is visible at roughly the 5 % level, which is too coarse to resolve a Molaro-sized 3–4 % change.
- **The pair is not isolated.** Fig. 2 shows both grains embedded in a bed of other ice particles. That buffers the local vapour toward saturation, but it also means the pair probably has other contacts.
- **Indirect evidence points to near-saturation.**
  - The neck grows as ≈ t^(1/3) for a full week and never recedes. Molaro's −5 °C necks did recede, under net sublimation.
  - By Demmenie et al.'s (2025) argument, under-saturation pulls the exponent toward 1/7. Thomas's 0.33–0.40 is at the saturated end.
- **The last two points are not evidence of shrinkage.** They read 0.443 then 0.440, which is within digitisation scatter.

**So: run h = 1.00** in two arms, a sealed box and Dirichlet walls, to bracket the boundary treatment. Optionally add one h = 0.998 arm to show how far Molaro-level under-saturation would move the exponent. All three cost the same (below).

---

## 4. The simulation (relaxed arm: resolve from point 3)

### Why "relaxed"

A neck is measurable only above x/a ≥ sqrt(12·eps/a). That floor goes at 0.9× the first point to be compared, so eps = a·(0.9·u)²/12.

Resolving point 1 (x/a = 0.099) needs eps = 0.080 µm: 24 M nodes, ~30k steps and $300–800. Starting at **point 3 (x/a = 0.181)** gives an eps 3.3× larger and 11× fewer nodes and steps. It still covers 9 of 11 points over a decade of time (20 → 190 h).

### Setup

The setup is the one validated in `studies/sinter_exponent`:
- Two equal spheres, axisymmetric, starting from exact tangency.
- α_c = 1e-3 held constant.
- Physical d₀ and `--Dchannel ice`.
- `-dtmax = 2·tau_sub`.

| | value |
|---|---|
| T | −20 °C |
| grain radius a | 120 µm (both grains) |
| floor x/a | 0.163 (0.9 × 0.181) |
| **eps** | **0.2657 µm** (h = eps/√2 = 188 nm) |
| domain | Lx × Ly = 528 × 144 µm (4.4a × 1.2a); grain centres at x = 144 and 384 µm |
| **mesh** | **Nx × Ny = 2811 × 767**: 2.16 M nodes, 6.5 M DoF, 33 ranks at 200k DoF/rank |
| tau_sub → dtmax | 555 s → **1110 s** |
| experiment duration | 3.3 → 190 h (~8 days); compared window 20 → 190 h |
| model time to point 3 / to last point (est.) | 23 h / 660 h (27 d) |
| model slower than experiment, over the window | **~3.7×** (Molaro's own −20 °C pair: ~150×) |
| **t_final** | **824 h = 2.97e6 s** (1.25 × last point) |
| **steps** | **~2,800** |
| **wall time** (8 / 15 / 23 s per step) | **6 / 12 / 18 h**, one job with no restart legs |
| **cost per arm** | **$2 / $5 / $7** (Tier 1, $0.012/core-h) |

How the estimates were made:

- **Model time.** Calibrated on the `mesh_pair` fine arm: −20 °C, α_c = 1e-3, tangent start, u = 0.194 → 0.303 in 15 → 79 h, so t ~ u^3.73. Scaled to this case by t ∝ β_sub·a²/d₀. Reaching u = 0.44 is an extrapolation, and exact-fillet slowing there makes the end time more likely an underestimate. **Treat the model times as ±2×.**
- **Seconds per step.** The low end is the Molaro axisymmetric runs (~7 s/step at ~100k DoF/rank). The high end is the phase-field strong-scaling fit (`project_hpc_cost_model`).
- **Steps do not depend on α_c.** Both the neck time and tau_sub scale as 1/α_c, so the step count ≈ (a/eps)²·u^3.7. A different α_c changes `t_final`, not the cost.

### Why tau_sub is 555 s here

With 99 % of the τ_sub bracket kinetic, `comp_eps.py`'s formula reduces to

    tau_sub ≈ eps² · beta_sub / d0,     beta_sub ∝ 1/α_c

so tau_sub scales with **eps²** and **1/α_c**:

| run | eps | α_c | tau_sub |
|---|---:|---:|---:|
| this run (Thomas relaxed) | 0.266 µm | 1e-3 | 555 s |
| `mesh_pair` fine (−20 °C) | 0.240 µm | 1e-3 | 452 s |
| Molaro tangent set (−20 °C) | 0.240 µm | 1.34e-2 | 35.6 s |
| Thomas strict (dropped) | 0.080 µm | 1e-3 | 50 s |

Check: (2.657e-7)² × 7.941e6 / 1.0152e-9 = 552 s, against 555 s from `comp_eps.py`. The difference is the thermal and vapour terms.

### Inputs to write (none created yet)

Geometry, `inputs/geometry/<family>/`:

```
-axisym 1  -ic_grain_union 1  -ic_type multi_grains  -dim 2
-Lx 5.28e-4  -Ly 1.44e-4  -Lz 0
-ice_grain_cx 1.44e-4,3.84e-4   -ice_grain_cy 0.0,0.0
-ice_grain_R 1.2e-4,1.2e-4  -ice_grain_ax 1.2e-4,1.2e-4  -ice_grain_ay 1.2e-4,1.2e-4
-Nx 2811  -Ny 767  -Nz 1
-eps 2.657e-7  -eps_valid_temp -20
-periodic 0
```

Experiment, `inputs/experiment/<family>/`. There is one file per arm; the arms differ only in boundary flags and humidity.

```
-t_final 2.97e6  -temp -20.0  -humidity 1.00  -grad_temp0 0.0,0.0,0.0
-alpha_pointwise 1  -alpha_model 0  -alpha_c0 1.0e-3  -alpha_lo 1.0e-3  -alpha_hi 1.0e-1
-d0_sub0 1.0152e-9
-dtmax 1110  -dtCFL 1  -dtCFL_dphimax 0.2
-t_out_log 80  -t_out_first 1
# Dirichlet arm: add -flag_BC_Tfix 1 -flag_BC_rhovfix 1
# optional under-saturated arm: Dirichlet + -humidity 0.998
```

### Analysis

```bash
python postprocess/neck_width.py <run_dir> --axisym          # --axisym is REQUIRED
python postprocess/fit_neck_growth.py <run>/neck_width.csv <data csv> \
    --anchor-neck-rn 0.1811 --demmenie
```

- `neck_width.py --axisym` returns the neck **width** (2x); divide by 2a for x/a.
- `fit_neck_growth.py` reads data series in the Molaro layout (`time_min,neck_size_um,...,lg_diam_um,sm_diam_um`). Convert the Thomas CSV first (neck = 2·(x/a)·120 µm, both diameters = 240 µm), or add a `time_h,relative_neck` branch to the loader.
- Report the model's `d_fixed` and `d_free` exponents over x/a = 0.181–0.443, next to the table in §2.
