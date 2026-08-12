# mesh_pair — tangent-contact Molaro grains at α_c = 1e-3 vs the −20 °C data

Both arms of `batch_2026-08-11__17.24.16_mesh_pair`. Identical grains, tangent
contact, α_c = 1e-3 pointwise, 79 h; **only ε and the mesh differ**.

| arm | ε [µm] | mesh | neck width in 78.6 h | reaches 32.81 µm at | its floor √(12ε/R) |
|---|---:|---|---|---:|---:|
| coarse | 0.867 | 631 × 198 | 13.7 → 53.5 µm | 7.74 h | 60.1 µm |
| fine | 0.240 | 2281 × 713 | 9.3 → 51.1 µm | 15.05 h | 31.6 µm |
| Molaro −20 °C | — | — | 32.8 → 64.8 µm (in 78 **min**) | 0 (by definition) | — |

Regenerate:

```bash
# the CSVs these all read — --axisym is REQUIRED (see below)
python postprocess/neck_width.py <run_dir> --axisym

python studies/sinter_exponent/analysis/compare_mesh_pair.py \
    --out studies/sinter_exponent/results/mesh_pair      # the headline figure
bash studies/sinter_exponent/analysis/run_mesh_pair_fits.sh   # the fit tables
```

> **`--axisym` is not optional.** These are axisymmetric r–z runs: the grid's y
> is the radius and the axis is y = 0, so the chord `neck_width.py` measures
> runs from the axis to the contour and *is* the neck radius. Without the flag
> every width is half its true value. An earlier revision of this folder had the
> fine arm extracted without it and concluded the two arms disagreed by 36–150×;
> that was the missing flag, not the physics.

---

## The figure

`meshpair_compare.png` — both arms, the Molaro points, and a dotted best fit
`w = C(t + t₀)^a` per arm. **Each simulation clock is shifted so that t = 0
falls where its own neck reaches 32.81 µm**, the experiment's first measured
width. Molaro et al. do not know when their grains touched, so their t = 0 is
already some t > 0 for a run that starts from true tangency; anchoring on a
common neck size is the only defensible shared origin, it fits nothing, and it
is what they did themselves in their Fig. 12.

## The numbers

| series | a | ±95 % | R² | n | window |
|---|---:|---:|---:|---:|---|
| coarse, post-anchor | 0.229 | 0.004 | 0.9999 | 24 | 32.8–53.5 µm |
| **fine, post-anchor** | **0.283** | 0.001 | 1.0000 | 17 | 32.8–51.1 µm |
| Molaro 2019 (−20 °C) | 0.204 | 0.053 | 0.9861 | 8 | 32.8–64.8 µm |

Benchmark: Demmenie et al. (2025), α = 0.25–0.34 (m = 3, sublimation–condensation).

**The fine arm is the one to quote.** Its resolution floor is 31.6 µm, so over
the compared range (32.8–51.1 µm) its neck is genuinely resolved — it is the
only arm here of which that is true, and the only one with a fittable window
under `fit_neck_growth.py`'s default protocol. It gives **a = 0.283**, inside
the Demmenie band. The coarse arm's entire run sits below its own 60.1 µm floor
and reads 0.229.

So the two arms do **not** give the same exponent: refining ε moves it from
0.229 toward 1/3. The rate, by contrast, is close to converged — at matched neck
width the arms differ by 1.9× at 32.8 µm, narrowing to 1.2× by 51 µm.

**Against the experiment, the exponent is the open question.** Molaro's data
gives 0.204 ± 0.053; the resolved arm gives 0.283 ± 0.001. Both are below
Demmenie's 1/3, but ours is inside their band and theirs is not, so the model
and this experiment are not measuring the same thing. Molaro's chamber humidity
was neither controlled nor measured and their grains net-shrank 3–4 %, which a
saturated cell cannot do — see `docs/molaro_validation_synthesis.md`.

## The rate offset

Both arms are ~150× slower than the experiment over the same neck-width range:

| arm | 32.81 µm → top | simulation | experiment | ratio |
|---|---|---:|---:|---:|
| coarse | 53.5 µm | 70.9 h | 28.0 min | 152× |
| fine | 51.1 µm | 63.6 h | 24.4 min | 156× |

That the two arms agree on this to 3 % is the strongest evidence the offset is
physical rather than numerical: it is the α_c = 1e-3 choice, the bottom of
Braun's range and deliberately deep in the kinetics-limited regime.

## Straight-line fit — `meshpair_linear.png`

A separate, blunter question: how many microns of neck per hour? Fitted on the
**coarse** arm only, with both clocks anchored at 32.81 µm and the intercept
pinned there, so the line has one free parameter:

| series | A [µm/h] | ±95 % | R² |
|---|---:|---:|---:|
| coarse (pinned @ anchor) | 0.349 | 0.027 | 0.890 |
| Molaro (pinned @ anchor) | 52.91 | 12.57 | 0.880 |

A line is *not* a description of the full run (R² = 0.77 over all 79 h, 10.4 µm
max residual) — `r ~ t^0.23` is concave — so A is quoted only over the width
range the experiment covers. The requested three-parameter form
`w = A(t − t₀) + C` is over-parameterized (`t₀` and `C` trade off exactly);
anchoring fixes the origin from the data instead, leaving the identifiable
one-parameter form. Pinned rate ratio 151.5× against the model-free elapsed-time
ratio 151.8× — two independent routes to the same number.

## Files

| file | what |
|---|---|
| `meshpair_compare.png` / `_fits.csv` | **the headline**: both arms + data, clocks anchored at 32.81 µm, dotted power-law fit per arm |
| `meshpair_coarse_*` | default protocol, floor enforced — the coarse arm has no fittable window, the fine arm does |
| `meshpair_full_*` | `--rmin 0`, each series over its own full range |
| `meshpair_shared_*` | both clipped to the r/R₀ = 0.197–0.309 overlap; the `kucz` fit there is degenerate, do not quote it |
| `meshpair_linear*` | straight-line fit to neck width, coarse arm only |
| `../../analysis/compare_mesh_pair.py` | driver for the headline figure |
| `../../analysis/run_mesh_pair_fits.sh` | driver for the power-law tables |
| `../../analysis/plot_neck_linear.py` | driver for the linear fit |
