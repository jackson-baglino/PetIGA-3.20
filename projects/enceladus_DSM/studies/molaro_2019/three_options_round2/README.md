# Round 2 at the corrected wall — and why the tuning campaign should stop

Batch `2026-09-08__17.20.46_molaro_T-20_round2`, three arms on dom2, all
completed. Measured with `postprocess/compare_arms_vs_molaro.py` →
[`summary.csv`](summary.csv), [`comparison.png`](comparison.png).
The cross-round result is [`tuning_ceiling.png`](tuning_ceiling.png), from
[`plot_tuning_ceiling.py`](plot_tuning_ceiling.py).

All numbers are on the anchored 78-minute window `[t*, t*+78 min]`. Molaro's
targets: neck **64.78 µm**, large grain **−2.93 %**.

## Round 2 as run

| arm | `-humidity` | 1−h | neck @78 | err | ΔR_large @78 | err | RMS | χ |
|---|---|---|---|---|---|---|---|---|
| 1 untuned | 0.99715 | 2.87e-03 | 45.98 µm | −18.80 | **−2.94 %** | −0.01 | 11.37 | 7.09 |
| 2a `D_v`×30 | 0.99928 | 7.43e-04 | 56.26 µm | −8.52 | −3.38 % | −0.46 | 4.79 | 3.13 |
| 2b `D_v`×100 | 0.99923 | 7.93e-04 | 54.55 µm | −10.23 | −3.96 % | −1.04 | 5.40 | 3.59 |

Two things to read off immediately.

**The window fix worked.** The untuned arm now reproduces the measured grain
recession essentially exactly: −2.94 % against −2.93 %. That was the whole
point of correcting the 120-min-vs-78-min bug, and the Newton step landed it
in one iteration.

**And `D_v` ×100 got worse on the neck**, 62.92 → 54.55 µm, exactly as
observed. That is not a regression — it is the trade-off finally being paid.
In round 1 that arm fit the neck *because* it was barely losing any grain
(−0.24 %); made to lose the right amount, it cannot.

## The result: transport tuning has a ceiling, and we are at it

Within each `D_v` family the two rounds differ **only** in the wall humidity.
So each family is a clean two-point sweep in the one free boundary parameter at
fixed transport, and Molaro's −2.93 % is a vertical cut through it:

| family | 1−h at the cut | neck at the cut | share of observed growth |
|---|---|---|---|
| untuned (nominal `D_v`) | 2.87e-03 | 45.98 µm | **41.2 %** |
| `D_v` ×30 | 6.45e-04 | **56.88 µm** | **75.3 %** |
| `D_v` ×100 | 5.91e-04 | **56.87 µm** | **75.3 %** |
| Molaro et al. (2019) | — | 64.78 µm | 100 % |

**×30 and ×100 land 0.01 µm apart.** A factor 3.3 in effective transport buys
nothing once the model is also required to reproduce the observed mass loss.
That is a ceiling, not a slow approach to the data.

This is the number that should be quoted, because it is the only fair one: a
model that grows a good neck while losing none of the grain is not reproducing
the experiment. Scored on both observables at once, the vapour-only model
delivers **75 % of the observed neck growth and stops**.

Two independent corroborations that this is real and not an artefact of the
two-point interpolation:

- **At nominal `D_v` the neck is completely decoupled from the wall.** The
  untuned arm's neck is 46.06 µm at 1−h = 1.86e-03 and 45.98 µm at
  2.87e-03 — and ~0.06 µm of that 0.08 µm gap is the interpolation-convention
  change below, not the wall, so a 54 % change in the wall moves the neck by
  roughly 0.02 µm. Neck growth
  is driven by the internal curvature difference, recession by the external
  wall, and at nominal transport they genuinely do not talk to each other. The
  untuned model's 18.7 µm neck deficit is therefore *not* a boundary-condition
  artefact. It is the physics.
- **The coupling turns on only at inflated `D_v`**, which is itself the tell:
  the wall reaches the neck efficiently only when transport is inflated well
  past the molecular value, so the very thing that made the neck fit in round 1
  is what makes the two observables fight in round 2.

## The exponent says the shape is wrong too, not just the amplitude

Fitting `w = C·(t' + t₀)^a` with `t₀` free (Demmenie's protocol) over the
anchored window:

| series | a |
|---|---|
| Molaro data | 0.235 ± 0.035 |
| untuned | 0.1349 ± 0.0007 |
| `D_v` ×30 | 0.1462 ± 0.0003 |
| `D_v` ×100 | 0.1222 ± 0.0009 |

The model grows too *slowly*, by ~2.5σ on their own scatter, and `D_v` does not
move it monotonically (×30 is above ×100). A mechanism that only rescales the
vapour flux cannot fix an exponent, so this is further evidence that what is
missing is a second transport channel rather than a bigger first one.

⚠️ `docs/molaro_validation_synthesis.md` §2 records "growth exponent n ~ 5
matches the data (x⁵ − x₀⁵ ~ t, R² = 0.99)". That is a different functional
form and a different protocol from the fit above, so the two are not directly
comparable — but they point different ways, and the discrepancy should be
resolved before either is published.

## Recommendation: stop tuning

Not as a concession — as the result. Three reasons:

1. **The ceiling is measured, not assumed.** Two different transport
   coefficients agree on it to 0.01 µm. Continuing to raise `D_v` is now a
   prediction we have already tested twice.
2. **The residual has a named mechanism.** Molaro's own Appendix A opens by
   calling vapour transport "one of the **two** dominant diffusion mechanisms
   driving neck growth in ice" — the other being surface diffusion, which this
   model does not carry. A vapour-only model reaching 75 % and stalling is the
   expected vapour share, not a failure to converge.
3. **Every further knob is less defensible than the last.** `alpha_c` is
   already exhausted (`alpha_c_estimate.csv`: ≤1.25× headroom above 0.1, and
   the attachment-only limit needs ≥0.139). The `M_0`×5 / `alph_sub`÷100 arm
   was dropped in round 1 for filling the neck by non-mass-conserving
   Allen–Cahn relaxation. What is left is fitting `D_v` per experiment, which
   round 2 has just shown does not even buy a better neck.

**The −5 °C series is not the next run.** It would test the same exhausted
model against the series `docs/molaro_validation_synthesis.md` §5 argues is
artefact-dominated. If the goal is a stronger validation claim, the next step
is a **surface-diffusion term in the phase-field model** — the future-work item
§4 of that document already identifies — and the number above is the
quantitative target it has to close: 7.9 µm of neck at 78 min, 25 % of the
observed growth.

## The oscillation on the neck curves was a measurement artefact

The growth curves carried a clean sinusoid riding on the power law. It is not
the simulation.

`.vts` snapshots are written on a **coarser grid than the solve** — 1080 × 541
against a 5394 × 2697 mesh — so the sample spacing is
`dy = 4.17e-7 m = 3.5·eps`. The sub-cell crossing in `chord_width()`
interpolated `phi` **linearly** across that, and linear interpolation of a
sigmoid over 3.5·eps is not a small correction: the two samples bracketing
`phi = 0.5` can sit at 0.03 and 0.97. The resulting error depends on where the
interface falls *inside* a cell, so a neck growing steadily outward sweeps
through sub-cell offsets at a steady rate — and the error comes out as a
sinusoid. It correlated with the radial sub-cell phase at |R| = 0.66–0.88.

**The fix is to interpolate in `logit(phi)`**, which is *exact* rather than
merely higher-order: this model's 1D equilibrium profile is
`phi = 1/(1 + exp(-s/eps))`, so `logit(phi) = s/eps` is **linear in distance**
and a straight line through two samples recovers the crossing at any spacing.

| | before (linear in φ) | after (logit) |
|---|---|---|
| fast-residual RMS, untuned | 0.0549 µm | 0.0232 µm (−58 %) |
| fast-residual RMS, `D_v`×30 | 0.0631 µm | 0.0371 µm (−41 %) |
| fast-residual RMS, `D_v`×100 | 0.0503 µm | 0.0335 µm (−34 %) |
| synthetic profile at 3.5·eps | 0.127 µm ptp | 0 (1.6e-16 µm) |

[`neck_interp_effect.png`](neck_interp_effect.png), from
[`plot_interp_artefact.py`](plot_interp_artefact.py); the analytic gate with no
simulation involved is `../verification/verify_neck_interp.py`.

**It changes nothing about the conclusions.** Final neck widths moved by
±0.06 µm, the ceiling moved by 0.01 µm (56.87/56.88 → 56.88/56.87 µm, 75.3 %
both), and every exponent moved by less than its own fit uncertainty. Also
tried and rejected: widening the axial vertex fit to a 5-point least-squares
parabola, which is a measured no-op (0.0232 → 0.0237 µm RMS). The residual
~0.03 µm is well below anything the conclusions rest on.

## Which round to keep: round 2

Round 1's raw batch has **already been deleted**; its numbers survive in
`../three_options/summary.csv`, which is all `plot_tuning_ceiling.py` needs.
Keep round 2's raw output. On the merits it is the better set regardless:

- its walls are calibrated on the **corrected** 78-minute window — round 1's
  were fitted against the full-run number and sit 1.54× too saturated;
- its untuned arm reproduces the measured recession (−2.94 % vs −2.93 %), so it
  is the defensible untuned-model-vs-experiment comparison; round 1's was 35 %
  short;
- round 1's *better* neck RMS was bought by barely losing any grain, which is
  not a fit to the experiment;
- only round 2 can carry the interpolation fix — round 1's snapshots are gone.

⚠️ One consequence worth recording: the round-1 values in
`../three_options/summary.csv` were measured with the **old** linear-in-φ
crossing and cannot be re-measured. The two ends of each family line in
`tuning_ceiling.png` therefore use slightly different conventions. The offset is
±0.06 µm against a 7.9 µm deficit, and recomputing the ceiling after the fix
moved it by 0.01 µm, so this does not affect any conclusion — but do not quote
round-1 neck widths to better than ~0.1 µm.

## What is worth one more run, if anything

A single **confirmation** arm at `D_v` ×30, 1−h = 6.45e-04, to replace the
two-point interpolation with a measurement. It should land at 56.9 µm and
−2.93 %. That converts the headline number from "interpolated" to "measured"
for the cost of ~6 h on 437 ranks, and it is the last vapour-only run worth
doing.
