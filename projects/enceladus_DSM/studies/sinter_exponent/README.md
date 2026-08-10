# Stage 1 results — do our Molaro simulations fit r ~ t^(1/3)?

Regenerate everything here with:

```bash
bash studies/sinter_exponent/verification/run_exponent_fits.sh
```

No simulation required: it reads the three committed Molaro neck curves from
`$RESULTS_BASE/_neck_csv/` and the validation data in `inputs/validation/`.

**Short answer: no, and the Molaro geometry cannot answer the question either
way.** Details below. See `PLAN.md` for the benchmark being tested against and
for the stages that follow.

---

## The numbers

From `verification/molaro_fits.csv`. `a` is the fitted exponent in
`r ~ t^a`; `m = 1/a` is the Kuczynski exponent (Demmenie's `m`, **not** their
mechanism index `n`).

| series | form | a | ±95 % | m | R² | window r/R0 |
|---|---|---:|---:|---:|---:|---|
| model eps 0.35 µm (strict) | `d_free` | 0.121 | 0.002 | 8.24 | 0.9997 | 0.220–0.275 |
| model eps 0.35 µm (strict) | `d_fixed` | 0.087 | 0.002 | 11.49 | 0.9901 | 0.220–0.275 |
| model eps 0.35 µm (strict) | `kucz` | 0.188 | 0.017 | **5.32** | 0.9914 | 0.220–0.275 |
| Molaro 2019 data (−20 °C) | `d_free` | 0.204 | 0.053 | 4.89 | 0.9861 | 0.197–0.380 |
| Molaro 2019 data (−20 °C) | `d_fixed` | 0.185 | 0.019 | 5.42 | 0.9837 | 0.197–0.380 |
| Molaro 2019 data (−20 °C) | `kucz` | 0.248 | 0.075 | **4.04** | 0.9742 | 0.197–0.380 |

Benchmark: Demmenie et al. (2025), **alpha = 0.29 / 0.33 / 0.30 / 0.26 ±0.01**,
i.e. m = 3.

The eps 0.60 µm and eps 0.86 µm arms are **absent from the table on purpose**.
Their resolution floors (0.289 and 0.345) sit above their entire measured range
(0.198–0.275), so they have no fittable window at all. The script says so on
stderr and plots them faint rather than dropping them silently.

## What it means

**1. Neither the Molaro data nor our model shows 1/3 at −20 °C.** The data's
local log-slope is flat at ≈ 0.185 — squarely between 1/7 and 1/5, not near
1/3. On the Kuczynski form it gives m = 4.04, consistent with the `n ~ 5`
recorded in `docs/molaro_validation_synthesis.md`. Whatever is happening in the
Molaro cryostage, it is not the clean evaporation–condensation Demmenie see.
That is not surprising: Demmenie's whole argument is that the literature's
scatter between 1/3 and 1/7 comes from imperfect humidity control, and the
Molaro chamber's humidity was neither controlled nor measured — their grains
measurably net-shrank 3–4 %, which a saturated cell cannot do.

**2. Our model never reaches a power law in these runs.** Its local slope
climbs monotonically from 0.02 to 0.12 across the full 2 h and is *still
rising* at the end (`verification/molaro_local_slope.png`). The run spends its
entire length in the fillet-rounding transient that follows a pre-formed neck.
The three eps arms lie on top of each other, so this is a property of the
initial condition, not of the mesh.

**3. The gap is much smaller on the form that accounts for the initial neck.**
Compare `d_fixed` — m = 11.5 model vs 5.4 data, an apparent factor of two — with
`kucz`, which subtracts r0: m = 5.32 vs 4.04. Most of the "disagreement" in the
naive fit is the model and the experiment having different amounts of neck at
t = 0, not different physics.

**4. The exponent is mostly an artifact of protocol.** `verification/
molaro_resampled_*` re-fits the *same* model curve at the *same* nine times the
experiment sampled. Nothing about the physics changes; only the sampling:

| sampling | `d_fixed` a |
|---|---|
| solver's dense output | 0.087 |
| at the data's nine times | 0.078 |
| `d_free` (t0 released) | 0.110 |

So any single quoted exponent is meaningless without its fit form and window.
This is why `fit_neck_growth.py` reports all three forms and records the window
in every row of its CSV.

**5. The real obstacle is geometric, and it is fixable.** The neck fillet has
radius `rho ~ r²/(2R)` and is resolved only above `r/R >= sqrt(12*eps/R)`. For
the Molaro pair that floor is 0.22, and the entire dataset spans 0.19–0.38 —
about a quarter of a decade, sitting on the floor. **You cannot distinguish
1/3 from 1/7 over a quarter decade.** The fix is not a finer mesh (the floor
only falls as `sqrt(eps)`, and eps is bounded below by cost) but **bigger
grains**: at Demmenie's R = 500 µm the floor drops to 0.087, which is also
about where Maeno & Ebinuma (1983) put the vapour/surface-diffusion crossover.

That is the argument for stage 2/3, and it is a stronger reason to run the
1 mm geometry than "the new paper used 1 mm".

## Files

| file | what |
|---|---|
| `molaro_fits.csv` | fit table: parameters, 95 % CIs, R², window, eps, R0, floor |
| `molaro_growth_loglog.png` | r/R0 vs t, log–log, with 1/1 · 1/3 · 1/5 · 1/7 guides, the Demmenie band, and each arm's resolution floor |
| `molaro_local_slope.png` | `d ln r / d ln t` vs r/R0 — the decisive figure |
| `molaro_resampled_*` | the same strict arm read at the data's sample times |
| `run_exponent_fits.sh` | driver; regenerates all of the above |

## Caveats

- The model arms are the 2 h union-IC runs of 2026-07-28
  (`2D_molaro_axisym_T-20pair_union_eps{loose,mid,strict}`), not new runs.
- `d_free` fits with only 6–9 points carry wide CIs; the `kucz` fit on the
  6-point resampled series has ±0.69 on `a` and should not be read as a result.
- The −5 °C Molaro series is deliberately not fitted. It is artifact-dominated
  (`docs/molaro_validation_synthesis.md` §5) and its own CSV header says so.
