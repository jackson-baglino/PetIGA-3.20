# Thomas et al. (1994) two-sphere sintering: data and exponent comparison

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
- *The line is the `d_fixed` fit over points 3–11, drawn only over that range. The dashed line and open/filled split mark the mesh floor of a run that was planned and then withdrawn (§4).*

Regenerate (from `enceladus_DSM/`):

```bash
python studies/historical_sintering/analysis/fit_power_law.py   # data/thomas_powerlaw_fits.csv
python studies/historical_sintering/analysis/plot_data.py       # figures/thomas_neck_data.*
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
| points 3–11 (late window) | `d_fixed` | **0.403** | 0.044 | — | 0.985 |
| points 3–11 | `d_free` | **0.329** | 0.096 | −9 ± 12 | 0.986 |

- **Over the full record the data agree with Thomas's t^(1/3).** a = 0.36 ± 0.04 contains 1/3.
- **Over points 3–11, a depends on the protocol.** With the clock as reported, a = 0.40 ± 0.04. Freeing the contact time gives 0.33 ± 0.10.
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

## 4. Decision (2026-10-02): no Thomas simulation; compare exponents instead

**We are not running Thomas.** Any version at the Molaro settings costs more than $15:
- α_c = 0.1 and ε = 1.18e-7 m (from δ_ice ≤ 1) give τ_sub = 1.33 s.
- The comparison needs ~184 h of simulated time, against Molaro's 2 h.

Every under-$15 variant gives up something the Molaro runs kept: a coarser ε (δ_ice ≈ 2.25), dtmax ≥ 20·τ_sub, or an untested mirror domain. The earlier sizing was done at α_c = 1e-3 and is withdrawn. `size_runs.py` and its CSV are in git history (`a4908306`).

**The comparison uses the growth exponent a in x ~ t^a.** Each source is fitted in its own window. `d_free` is C·(t+t₀)^a with t₀ free, and is the fair form when the contact time is unknown.

| source | conditions | window | a (`d_free`) | a (`d_fixed`) |
|---|---|---|---:|---:|
| **Thomas et al. (1994)** data | −20 °C, cold lab, humidity unreported | x/a 0.10–0.44, 3–190 h | 0.39 ± 0.08 (all) · 0.33 ± 0.10 (pts 3–11) | 0.36 ± 0.04 (all) · 0.40 ± 0.04 (pts 3–11) |
| Thomas's own fit | — | — | — | **1/3** |
| **Molaro et al. (2019)** data | −20 °C cryostage, grains shrank 3–4 % | 33 → 65 µm width, 78 min | 0.204 ± 0.053 | 0.185 ± 0.019 |
| **our model**, Molaro −20 °C final run | α_c = 0.1, h = 0.99715 (shrinkage-calibrated) | 33.5 → 48.5 µm, 117 min after anchor | **0.139 ± 0.001** | — |
| **our model**, Molaro −5 °C final run | α_c = 0.1, h = 0.99674 | 34.4 → 53.3 µm, 54 min after anchor | **0.147 ± 0.002** | — |
| our model, saturated, neck resolved (`mesh_pair` fine) | α_c = 1e-3, h = 1.00, tangent start | 32.8 → 51.1 µm | 0.283 ± 0.001 | — |
| Demmenie et al. (2025) data | −3 °C, held at ice saturation | r/R 0.09–0.35 | 0.26–0.33 | — |
| Demmenie et al. (2025), under-saturated control | — | — | ≈ 1/7 | — |

- **Model exponents** come from `neck_width.csv` in `~/SimulationResults/HPC_results/enceladus_DSM/GrainPairSintering/`:
  - −20 °C: second attempt (`batch_2026-09-08__17.20.46_molaro_T-20_round2`).
  - −5 °C: first attempt (`..._T-5_h0.99674_2h_a1e-1_dirichlet`).
  - Each clock is anchored at Molaro's first measured width.
- **Data exponents** come from `analysis/fit_power_law.py` and `studies/sinter_exponent/README.md`.

**How the exponents line up:**
- **Saturation sets the exponent, as Demmenie et al. argue.** Our model gives a ≈ 0.28–0.33 at saturation. At Molaro's shrinkage-calibrated humidity (h ≈ 0.997) it drops to a ≈ 0.14–0.15, close to Demmenie's under-saturated 1/7.
- **The Molaro data sit between the two (0.20).** That fits a chamber that was slightly under-saturated, which their 3–4 % grain shrinkage already shows.
- **Thomas sits at or above the saturated value (0.33–0.40).** By the same argument, their grains were at least at saturation.
  - It is consistent with a vapour-saturated environment: the pair sits in a bed of other ice particles in a cold room, and the neck never recedes over a week.
  - It does **not** show super-saturation. Only the clock-as-reported fit over points 3–11 (0.40 ± 0.04) exceeds 1/3; with t₀ free it is 0.33 ± 0.10.
  - Calling it super-saturated would need evidence the paper doesn't give, such as grain growth, which Fig. 2 can't resolve (§3).
- **Rate.** At the α_c = 0.1 used for Fig. 2, the model reproduces about half of Molaro's −20 °C neck growth. That missing half matches the surface-diffusion share SA81 predicts. The two mechanisms are co-dominant in Molaro's own mechanism breakdown (their Fig. 10), and **our model has no surface diffusion** (`project-molaro-validation`, `docs/molaro_validation_synthesis.md`).
  - Molaro's model *does* include surface diffusion. Its weak point is the vapour term: Hertz–Knudsen with α = 1 in vacuum, ignoring gas-phase diffusion (their Appendix A).
  - We have **no** reliable model-vs-Thomas rate ratio. The "~3.5× slower" estimate was made at α_c = 1e-3 and does not carry over to 0.1.

**Caveats.**
- The model runs are pre-necked (r = 14 µm). The early fillet transient is absorbed by the free t₀, not removed, so the model exponents are best read as ±0.03.
- Exponents fitted over different x/a windows are not strictly comparable: the ideal saturated slope itself sags from 0.33 to about 0.30 between x/a 0.1 and 0.35.
- Thomas's and Molaro's points are both second-hand digitisations.
