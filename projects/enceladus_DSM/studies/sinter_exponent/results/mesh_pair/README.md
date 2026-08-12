# mesh_pair — tangent-contact Molaro grains at α_c = 1e-3 vs the −20 °C data

Regenerate with:

```bash
bash studies/sinter_exponent/analysis/run_mesh_pair_fits.sh
```

Both arms of `batch_2026-08-11__17.24.16_mesh_pair`. Identical grains, tangent
contact, α_c = 1e-3 pointwise, 79 h; **only ε and the mesh differ**:

| arm | ε [µm] | mesh | neck width reached in 78.6 h | samples | its floor √(12ε/R) |
|---|---:|---|---|---:|---:|
| coarse | 0.867 | 631 × 198 | 13.7 → **53.5** µm | 53 | 60.1 µm |
| fine | 0.240 | 2281 × 713 | 4.7 → **25.6** µm | 57 | 31.6 µm |
| Molaro −20 °C | — | — | 32.8 → 64.8 µm (in 78 **min**) | 9 | — |

---

## Headline: the two arms do not agree, so the pair's question is answered — badly

**The coarse mesh is not adequate.** At matched neck width the arms differ by a
factor of **36–150 in rate**, and re-zeroing them at a shared neck does *not*
collapse them (`meshpair_compare.png` panel b), so this is not the coarse arm
merely being "ahead" on the same trajectory. From
`meshpair_compare_rate_ratio.csv`:

| neck width [µm] | coarse [h] | fine [h] | fine/coarse |
|---:|---:|---:|---:|
| 13.9 | 0.05 | 7.6 | 150 |
| 17.8 | 0.37 | 20.7 | 56 |
| 21.7 | 1.00 | 43.2 | 43 |
| 25.6 | 2.19 | 78.6 | 36 |

**The exponent is not converged either, once the window is honest.** Fitted over
each arm's *own* range the two look reassuringly similar (a = 0.19–0.21 coarse,
0.21–0.24 fine) — but they are covering different parts of the trajectory. Over
the width range the arms actually share (13.7–25.6 µm):

| series | `d_free` | `d_fixed` | `kucz` |
|---|---:|---:|---:|
| coarse | 0.238 | 0.157 | 0.280 |
| fine | **0.285** | **0.261** | 0.592 † |
| Demmenie 2025 band | 0.25–0.34 | | |

† parked near the fit bound; degenerate, do not quote.

Refining ε moves the exponent **toward** 1/3: the fine arm's 0.26–0.285 sits in
the Demmenie band, the coarse arm's 0.16–0.24 does not. This is the same
same-window discipline that `results/dv_sweep/` documents — there a wider window
on one arm manufactured a difference; here it manufactured an *agreement*.

**Neither arm resolves its own neck.** The coarse run ends at 53.5 µm against a
60.1 µm floor; the fine run ends at 25.6 µm against a 31.6 µm floor. Every
number on this page is therefore below √(12ε/R) — see the caveats section. The
fine arm is closer to its floor and to the Demmenie band, which is consistent
with the coarse arm's neck being substantially an ε-scale artefact (the two
arms' *initial* widths, 13.7 vs 4.7 µm, scale roughly as √ε as an unresolved
diffuse bridge should).

## The requested alignment: possible for the coarse arm only

Shifting the clock so the model's neck equals the experiment's first measured
width (32.81 µm) at t = 0 — Molaro et al.'s own Fig. 12 convention — works for
the coarse arm, which reaches that width at **7.74 h**. It is **impossible for
the fine arm**: that run tops out at 25.56 µm, 7.25 µm *below* the experiment's
starting neck, so the two have **no common neck width at all**. Extrapolating
its `d_free` fit puts the crossing at ≈ 325 h, 4× the length of the run — a
number worth knowing and not worth plotting, so `meshpair_compare.png` panel (c)
annotates it rather than drawing a curve.

---

## Coarse-arm results (superseded as a converged result; kept as the record)

Everything below is the coarse arm alone. Given the above it should be read as
"what the unconverged arm says", not as a model–experiment comparison.

From `meshpair_full_fits.csv` (`--rmin 0`, each series over its own full range):

| series | form | a | ±95 % | m | R² | window r/R₀ | n |
|---|---|---:|---:|---:|---:|---|---:|
| model coarse (tangent) | `d_free`  | 0.208 | 0.002 | 4.80 | 0.9997 | 0.079–0.309 | 53 |
| model coarse (tangent) | `d_fixed` | 0.189 | 0.004 | 5.30 | 0.9933 | 0.079–0.309 | 53 |
| model coarse (tangent) | `kucz`    | 0.211 | 0.003 | **4.74** | 0.9991 | 0.079–0.309 | 53 |
| Molaro 2019 (−20 °C)   | `d_free`  | 0.204 | 0.053 | 4.89 | 0.9861 | 0.197–0.380 | 8 |
| Molaro 2019 (−20 °C)   | `d_fixed` | 0.185 | 0.019 | 5.42 | 0.9837 | 0.197–0.380 | 8 |
| Molaro 2019 (−20 °C)   | `kucz`    | 0.248 | 0.075 | **4.04** | 0.9742 | 0.197–0.380 | 8 |

Benchmark: Demmenie et al. (2025), α = 0.29/0.33/0.30/0.26 ±0.01, i.e. m = 3.

This is a large change from the pre-necked runs in `../molaro_prenecked/`, where
the same model gave a = 0.087–0.121 and never reached a power law at all. The
difference is the initial condition, not the physics: **starting at tangent
contact gives the run a real clock zero**, so it spends 53 samples on the power
law instead of spending the whole run rounding a fillet it was handed. That also
makes `d_fixed` (t₀ ≡ 0) legitimate here for the first time, and all three forms
now agree with each other to ±0.02 instead of spanning a factor of two.

## Three caveats, in order of how much they matter

**1. The whole curve sits below the coarse arm's resolution floor.** The fillet
is resolved only above r/R₀ ≥ √(12ε/R) = **0.346** (neck width 60.1 µm); the run
ends at 0.309 (53.5 µm). Under the default protocol the model therefore has *no
fittable window at all* — that is what `meshpair_coarse_*` shows, and it is
exactly what `batches/mesh_pair.txt` predicted. The numbers above come from
overriding the floor. **The fine arm does not rescue this**: its floor is 31.6 µm
and it only reaches 25.6 µm, so it is below its own floor too, over its whole
run. The pair was expected to give one resolved arm; it gave none.

**2. The rate is off by a factor of ~155.** Growing the neck from 33.6 to
53.5 µm width takes the model 69.7 h and took Molaro 27 min. That is the
α_c = 1e-3 choice — the bottom of Braun's range, deliberately deep in the
kinetics-limited regime — not a fit failure. On the log–log figure it is the
~2-decade horizontal offset between the blue curve and the orange points; the
two run *parallel*, which is the point.

**3. The first ~15 samples are interface relaxation, not sintering.** The run
opens at width 13.7 µm because two tangent circles with a 6ε diffuse band
overlap on contact. The local slope climbs from 0.08 to its 0.20 plateau over
r/R₀ = 0.08 → 0.13 (`*_local_slope.png`). `d_free` absorbs this with t₀ = 265 s;
`d_fixed` does not, which is why it reads 0.19 rather than 0.21.

## The decisive figure

`meshpair_full_local_slope.png` — d ln r / d ln t against r/R₀, model and data on
the same axes with the exact-fillet ideal drawn:

- model plateaus at **0.20**, rising to 0.23 by the end of the run;
- data sits flat at **0.18–0.20** across its whole range;
- the exact-fillet kinetic-limited ideal runs 0.33 → 0.29;
- Demmenie's band is 0.25–0.34.

Model and experiment lie on top of each other, and both lie a clear 0.10 below
what a perfect kinetic-limited vapour model should give. So the deficit is
*shared*, which looked like the useful result.

**The fine arm withdraws that conclusion.** Its local slope over the same widths
is 0.26–0.285, not 0.20 — so the "shared deficit" was a property of the coarse
mesh, not of the physics. What survives is weaker and still worth having: the
experiment's exponent is ≈ 0.19, the best-resolved arm gives ≈ 0.27, and the gap
between them is now the open question rather than the agreement.

## Straight-line fit — `meshpair_linear.png`

A separate, blunter question: **how many microns of neck per hour, and does the
rate hold steady?** `analysis/plot_neck_linear.py` fits `w = A·t + C` to the neck
*width* on linear axes and overlays the Molaro points.

**Both clocks are anchored at a common starting neck.** The experiment opens
with a 32.81 µm neck already formed; the model opens at tangent contact and
takes **7.74 h** to reach that width. Comparing them on raw clocks compares two
different parts of the same trajectory. Re-zeroing both at w = 32.81 µm removes
the mismatch without fitting anything away, and is the convention Molaro et al.
used in their own Fig. 12. Times below are *since the anchor*.

| series | form | A [µm/h] | ±95 % | C [µm] | R² | rms | max resid | window [µm] | n |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|
| model (full run) | free | 0.494 | 0.074 | 27.9 | 0.773 | 5.40 | 10.39 | 13.7–53.5 | 53 |
| model (shared width) | free | 0.290 | 0.025 | 35.3 | 0.960 | 1.23 | 2.33 | 33.0–53.5 | 24 |
| **model (pinned @anchor)** | pinned | **0.349** | 0.027 | 32.81 | 0.890 | 2.04 | 4.04 | 33.0–53.5 | 24 |
| Molaro (full range) | free | 24.89 | 6.52 | 37.1 | 0.889 | 3.76 | 5.05 | 32.8–64.8 | 9 |
| Molaro (shared width) | free | 50.78 | 20.97 | 33.3 | 0.883 | 1.97 | 3.80 | 32.8–46.8 | 5 |
| **Molaro (pinned @anchor)** | pinned | **52.91** | 12.57 | 32.81 | 0.880 | 1.99 | 3.93 | 32.8–46.8 | 5 |

Anchoring is what makes the intercept stop being a fitted nuisance: pinned at
w = 32.81 µm the line has **one** free parameter, the rate, so the two `A`
values are the comparison with nothing else moving. It costs a little R²
(0.89 vs 0.96 free) because the concavity can no longer be partly absorbed into
an offset — that is the honest number, not a worse fit.

**The pinned fit also validates the time-scale factor.** Pinned rate ratio
A_data/A_model = **151.5×**, against the model-free elapsed-time ratio of
**151.8×** — agreement to 0.2 %. The free-slope ratio reads 175× only because a
free intercept lets the data's truncated window (32.8–46.8 µm) tilt its chord.
Two independent routes to the same number is the strongest evidence here that
the model–experiment difference really is a single rate factor.

**A line is not a description of the full run.** Over all 79 h it gives R² = 0.77
with a 10.4 µm max residual on a 39.9 µm span — 26 % — because `r ~ t^0.2` is
strongly concave and a chord cannot follow it. Restricted to the width range the
experiment actually covers, a line is a fair local approximation, and
**A = 0.349 ± 0.027 µm/h** (pinned) is the number worth quoting. The residual
panel shows the leftover arc, which is the concavity the line cannot absorb,
at ±4 µm.

**The requested `w = A·(t − t0) + C` is not fitted as written**, because it is
over-parameterized: it expands to `A·t + (C − A·t0)`, so `t0` and `C` trade off
exactly and only two of the three parameters are identifiable. Anchoring is the
principled way to spend that third parameter — it fixes the origin *from the
data* instead of asking the fit to invent one. The script reports the free line
both ways (slope-intercept `(A, C)` and zero-crossing `t0 = −C/A`) alongside the
pinned fit. The free `t0` column is in the CSV for completeness but is a large
negative number (−122 h for the model) that says only "extrapolating a chord
backwards off a concave curve misses the origin", not anything about when the
grains touched. **Quote the pinned `A`.**

**The collapse is the real result.** Stretching the Molaro clock by the single
factor **S = 152** — the elapsed-time ratio at equal neck size, model 70.9 h vs
data 28.0 min across 32.8–53.5 µm — puts the experiment's points on the model
curve within their own error bars. One scalar, no shape adjustment. Combined
with the matching exponents above, the model differs from the −20 °C experiment
by a pure rate factor over this range, and that factor is the α_c = 1e-3 choice.

Two cautions on S, both about windowing rather than physics. It is quoted from
elapsed times, not from the ratio of *free* fitted slopes (which reads 175):
with only 9 experimental points the data's fit window cannot land on the window
edge and truncates at 46.8 µm instead of 53.5, and with a free intercept a chord
over that shorter, steeper sub-range overstates the ratio. Pinning the intercept
removes that freedom and recovers 151.5×. Earlier drafts also clipped the shared
window on the model only, leaving the data's chord spanning 32.8–64.8 µm against
the model's 32.8–53.5; because both curves are concave that mismatch pulled S
down to 86, off by 1.8×. The window must be the intersection, applied to both.

## Why there is no `--anchor-neck` variant here

Anchoring re-zeroes the clock at a shared neck size, which is how you overlay a
run that opened with a mature neck onto one that did not (Molaro's own Fig. 12
convention, used throughout `../molaro_prenecked/`). These runs start at true
tangent contact, so t = 0 is physical; anchoring would discard the one property
that makes `d_fixed` meaningful, and it drops the model's entire early rise. It
was tried and gives a = 0.23 (`d_free`) / 0.066 (`d_fixed`) — the split is a
protocol artifact, and the files were removed rather than left to be misread.

## Files

| file | what |
|---|---|
| `meshpair_coarse_*` | **default protocol.** Floor enforced → the model has no fittable window. |
| `meshpair_full_*` | **headline.** `--rmin 0`, each series over its own full range. |
| `meshpair_shared_*` | both clipped to the r/R₀ = 0.197–0.309 overlap; leaves the data 3 points, and the `kucz` fit there (m = 1.34) is degenerate — do not quote it. |
| `meshpair_linear*` | straight-line fit to neck **width** on linear axes, both clocks anchored at w = 32.81 µm, Molaro overlaid on a ×152 clock, with a residual panel |
| **`meshpair_compare.png`** | **both arms + data: raw clocks, arms re-zeroed at a shared neck, and aligned on the experiment's first width** |
| `meshpair_compare_fits.csv` | power-law fits for both arms over their own *and* their shared window |
| `meshpair_compare_rate_ratio.csv` | fine/coarse elapsed time at matched neck width |
| `meshpair_both_*` | both arms + data through `fit_neck_growth.py`, each over its own range |
| `meshpair_aligned_*` | `--anchor-neck 1.6405e-5`: coarse + data on a common clock (the fine arm is reported as never reaching it) |
| `../../analysis/run_mesh_pair_fits.sh` | driver; regenerates the power-law sets |
| `../../analysis/plot_neck_linear.py` | driver; regenerates the linear set |
| `../../analysis/compare_mesh_pair.py` | driver; regenerates the arm-vs-arm comparison |
