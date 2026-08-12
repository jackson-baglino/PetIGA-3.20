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
| `../../analysis/run_mesh_pair_fits.sh` | driver; regenerates all of the above |
