# Experimental ice-sintering data to compare against

Compiled 2026-09-30 for the Molaro-validation discussion. **Every DOI below
was resolved on Crossref** (title and first author matched) before it went in
this file; content notes come from the abstract or from the source named, and
anything not yet read in the paper itself is marked *(check)*.

## The temperature argument, in numbers

Vapour-transport sintering scales with the vapour the gas can carry,
~ rho_vs(T) * D_v(T) (both from `preprocess/comp_eps.py`):

| T [°C] | rho_vs [kg m⁻³] | rho_vs / rho_vs(−20) | D_v / D_v(−20) |
|---|---|---|---|
| −5 | 3.311e-3 | 3.90 | 1.11 |
| −10 | 2.139e-3 | 2.52 | 1.07 |
| −20 | 8.487e-4 | 1 | 1 |
| −40 | 1.055e-4 | 0.124 | 0.86 |
| −80 | 4.501e-7 | 5.3e-4 (1/1890) | 0.61 |

Minutes after t* for the neck to reach a given width:

| width [µm] | 38 | 42 | 44 | 46 |
|---|---|---|---|---|
| model −20 °C | 13.2 | 36.8 | 54.8 | 78.3 |
| model −5 °C (round 1) | 2.9 | 8.5 | 12.9 | 18.5 |
| **model ratio −20/−5** | **4.6** | **4.3** | **4.2** | **4.2** |
| Molaro −20 °C | 5.4 | 7.6 | 8.6 | 14.5 |
| Molaro −5 °C | 15.8 | 29.9 | 48.0 | — |
| **Molaro ratio −20/−5** | **0.34** | **0.25** | **0.18** | — |

The model's 4.2–4.6x is the vapour scaling (3.90 x 1.11 = 4.3) and the
k_eff campaign's tau_sub ratio (7822 / 2063 = 3.8). Molaro's −5 °C pair is
2.9–5.6x *slower*: against the model's temperature trend by a factor
~12–25. The two grain pairs are nearly the same size (202/145 vs 197/153 µm),
so geometry does not explain it. Molaro et al. themselves note the −5 °C grains
sintered slower than their own model predicted and sublimated faster
(PMC7458059); see `docs/molaro_validation_synthesis.md` §5 for the heating
artefact that explains both.

At −80 °C a vapour-driven process runs ~1/1900 as fast as at −20 °C
(~1/3100 with D_v): a day at −20 °C is ~5–8 years at −80 °C, which is why
cold storage preserves microstructure.

## Two-grain neck growth (compare with the Molaro setup)

- **Hobbs & Mason (1964)**, Phil. Mag. 9(98), 181.
  [10.1080/14786436408229184](https://doi.org/10.1080/14786436408229184).
  Neck growth between ice spheres, 50–700 µm diameter, **−3 to −20 °C**, in
  air at atmospheric pressure. Quantitative agreement with
  evaporation–condensation *including diffusion through the gas*; volume or
  surface diffusion argued ~4 orders of magnitude slower. Molaro et al. say
  their model "most closely matches" these data. **The closest match to our
  model** (vapour transport through air, two spheres, a temperature range
  spanning both Molaro series). A correction was published:
  [10.1080/14786436408211918](https://doi.org/10.1080/14786436408211918).
  *(check: their tabulated B(T) coefficients / neck-vs-time data)*
- **Kingery (1960)**, J. Appl. Phys. 31, 833.
  [10.1063/1.1735704](https://doi.org/10.1063/1.1735704).
  Ice-particle welding down to −25 °C; attributed to surface diffusion from
  the time and size dependence. The competing-mechanism reference.
- **Kuroiwa (1961)**, Tellus 13(2), 252.
  [10.1111/j.2153-3490.1961.tb00082.x](https://doi.org/10.1111/j.2153-3490.1961.tb00082.x).
  Neck growth between ice spheres (radius < 100 µm) in **saturated air and in
  kerosene** -- kerosene suppresses vapour transport, isolating solid-state
  diffusion. Useful for the vapour vs surface-diffusion split
  (`molaro_validation_synthesis.md` §4). *(check: temperature range)*
- **Chen & Baker (2010)**, Hydrol. Process.
  [10.1002/hyp.7787](https://doi.org/10.1002/hyp.7787).
  2-D ice-sphere arrays, quasi-isothermal (−10 °C) and under a gradient;
  neck morphology differs with the thermal condition. Qualitative.

## Packings / snow under isothermal conditions (compare with k_eff campaign)

- **Kaempfer & Schneebeli (2007)**, JGR Atmospheres.
  [10.1029/2007JD009047](https://doi.org/10.1029/2007JD009047).
  Fresh snow held isothermal at **−1.6, −8.3, −19.1 and −54 °C for ~1 year**,
  micro-CT monthly: density, SSA, SMI, trabecular structure. At −54 °C SSA fell
  only 19 % in a year. **The best temperature-series benchmark for the k_eff
  campaign**: its SSA(t) at four temperatures tests the model's claim that T
  only rescales time (t / tau_sub). Also direct evidence for the cold-storage
  point.
- **Legagneux, Taillandier & Domine (2004)**, J. Appl. Phys. 95, 6175.
  [10.1063/1.1710718](https://doi.org/10.1063/1.1710718).
  Isothermal SSA decay of four snows at **−15 °C only**, closed systems;
  SSA = B − A ln(t + Δt). Tests the SHAPE of SSA(t), not the T dependence.
- **Flin et al. (2004)**, Ann. Glaciol.
  [10.3189/172756404781814942](https://doi.org/10.3189/172756404781814942).
  3-D micro-CT measurements of isothermal snow evolution. *(check: T)*
- **Willibald et al. (2019)**, Front. Earth Sci.
  [10.3389/feart.2019.00229](https://doi.org/10.3389/feart.2019.00229).
  2 mm ice beads vs rounded snow at −10 °C, 1–10 d; stiffness and strength
  from CT + FEM. Mechanical, single temperature.

## Reviews

- **Blackford (2007)**, J. Phys. D 40, R355.
  [10.1088/0022-3727/40/21/R02](https://doi.org/10.1088/0022-3727/40/21/R02).
  Review of ice sintering and microstructure; the place to find further
  temperature-series data. *(check)*

## Also in the repo

- Demmenie et al. (2025), J. Phys. Chem. Lett.
  [10.1021/acs.jpclett.5c00050](https://doi.org/10.1021/acs.jpclett.5c00050),
  "Ice Sintering by Sublimation and Condensation" -- `studies/sinter_exponent/`.
- Molaro et al. (2019) [10.1029/2018JE005773](https://doi.org/10.1029/2018JE005773).
