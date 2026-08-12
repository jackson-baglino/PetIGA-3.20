# mesh_pair — tangent-contact Molaro grains at α_c = 1e-3 vs the −20 °C data

Regenerate with:

```bash
bash studies/sinter_exponent/analysis/run_mesh_pair_fits.sh
```

Run analysed: `batch_2026-08-11__17.24.16_mesh_pair/molaro_2D_L387um_eps0.87um_axisym_T-20pair_tangent__molaro_T-20_h1.00_79h_a1e-3_coarse`
(631 × 198, ε = 8.675e-7 m, R̄ = 86.75 µm, α_c = 1e-3 constant and pointwise,
79 h, 53 neck samples, neck width 13.7 → 53.5 µm).

**The fine arm (ε = 0.24 µm) is not in this folder** — only the coarse arm was
pulled off the cluster. The driver picks the fine arm up automatically once its
folder appears, and the comparison the batch was designed to make is not
answered until it does.

---

## Headline

**The tangent-start run reproduces the Molaro exponent — and it is *not* 1/3.**
Every fit form gives the model and the data the same answer to within the data's
95 % CI, at m ≈ 4.7–5.4 rather than Demmenie's m = 3.

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
overriding the floor. The fine arm's floor is 0.182 (31.6 µm), below the entire
data range, which is why it is the arm that settles this.

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
*shared*, which is the useful result: whatever suppresses the exponent in the
Molaro cryostage, this model reproduces it at α_c = 1e-3. It is not yet
established that the mechanism is the same one.

## Straight-line fit — `meshpair_linear.png`

A separate, blunter question: **how many microns of neck per hour, and does the
rate hold steady?** `analysis/plot_neck_linear.py` fits `w = A·t + C` to the neck
*width* on linear axes and overlays the Molaro points.

| series | A [µm/h] | ±95 % | C [µm] | R² | rms | max resid | window [µm] | n |
|---|---:|---:|---:|---:|---:|---:|---|---:|
| model (full run) | 0.494 | 0.074 | 24.0 | 0.773 | 5.40 | 10.39 | 13.7–53.5 | 53 |
| model (shared width) | 0.290 | 0.025 | 33.0 | 0.960 | 1.23 | 2.33 | 33.0–53.5 | 24 |
| Molaro (full range) | 24.89 | 6.52 | 37.1 | 0.889 | 3.76 | 5.05 | 32.8–64.8 | 9 |
| Molaro (shared width) | 50.78 | 20.97 | 33.3 | 0.883 | 1.97 | 3.80 | 32.8–46.8 | 5 |

**A line is not a description of the full run.** Over all 79 h it gives R² = 0.77
with a 10.4 µm max residual on a 39.9 µm span — 26 % — because `r ~ t^0.2` is
strongly concave and a chord cannot follow it. Restricted to the width range the
experiment actually covers, a line is a fair local approximation (R² = 0.96,
rms 1.2 µm), and **A = 0.290 ± 0.025 µm/h** is the number worth quoting. The
residual panel shows the leftover arc, which is the concavity the line cannot
absorb, at ±2 µm.

**The requested `w = A·(t − t0) + C` is not fitted as written**, because it is
over-parameterized: it expands to `A·t + (C − A·t0)`, so `t0` and `C` trade off
exactly and only two of the three parameters are identifiable. The script fits
the line once and reports both identifiable readings — slope-intercept `(A, C)`
and zero-crossing `t0 = −C/A`. Use `C`; the `t0` column is in the CSV for
completeness but is a large negative number (−114 h for the model) that says
only "extrapolating a chord backwards off a concave curve misses the origin",
not anything about when the grains touched.

**The collapse is the real result.** Stretching the Molaro clock by the single
factor **S = 152** — the elapsed-time ratio at equal neck size, model 70.9 h vs
data 28.0 min across 32.8–53.5 µm — puts the experiment's points on the model
curve within their own error bars. One scalar, no shape adjustment. Combined
with the matching exponents above, the model differs from the −20 °C experiment
by a pure rate factor over this range, and that factor is the α_c = 1e-3 choice.

Two cautions on S. It is quoted from elapsed times, not from the ratio of fitted
slopes (which reads 175): with only 9 experimental points the data's fit window
cannot land on the window edge and truncates at 46.8 µm instead of 53.5, and a
chord over that shorter, steeper sub-range overstates the ratio. Earlier drafts
of this comparison also clipped the shared window on the model only, leaving the
data's chord spanning 32.8–64.8 µm against the model's 32.8–53.5; because both
curves are concave that mismatch pulled S down to 86, off by 1.8×. The window
must be the intersection, applied to both.

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
| `meshpair_linear*` | straight-line fit to neck **width** on linear axes, Molaro overlaid on a ×152 clock, with a residual panel |
| `../../analysis/run_mesh_pair_fits.sh` | driver; regenerates the power-law sets |
| `../../analysis/plot_neck_linear.py` | driver; regenerates the linear set |
