# Three tuning options at T = −20 °C — results

Batch `2026-09-03__17.32.48_molaro_three_options_T-20`, four arms on dom2
(43.6 M DoF, 437 ranks), all completed, 121 minute-cadence snapshots each.
Measured with `postprocess/compare_arms_vs_molaro.py` →
[`summary.csv`](summary.csv), [`comparison.png`](comparison.png).

**Everything below is on the anchored 78-minute window** `[t*, t*+78 min]`,
where `t*` is when the model neck first reaches Molaro's first measurement
(32.81 µm). Their targets are neck **64.78 µm** and large grain **−2.93 %**.

| arm | `-humidity` | neck @78 min | err | ΔR_large @78 min | err | neck RMS | χ |
|---|---|---|---|---|---|---|---|
| 1 untuned | 0.99816 | 46.06 µm | **−18.72** | −1.90 % | +1.03 | 11.34 µm | 7.08 |
| 2a `D_v` ×30 | 0.99992 | 60.26 µm | −4.52 | −0.41 % | +2.52 | 3.20 µm | 2.18 |
| 2b `D_v` ×100 | 0.99996 | **62.92 µm** | **−1.86** | −0.24 % | +2.69 | **2.24 µm** | **1.88** |
| 3 `M_0`×5, `alph`÷100 | 0.99264 | 75.48 µm | **+10.70** | −4.78 % | −1.86 | 11.20 µm | 8.97 |

χ is the RMS residual in units of their own (asymmetric) error bars, so χ ≈ 1
means the model sits inside the measurement scatter.

Wall time: 186 / 342 / 216 / 754 min. Arm 3 is by far the most expensive, as
`mob_scale 5` tightens the interface-CFL limit.

## 1. The best fit is D_v ×100, on the neck alone

Arm 2b tracks all nine of their points at 2.24 µm RMS against error bars of
±1.7–2.4 µm — i.e. within the measurement scatter — and lands 1.86 µm short at
78 min. Arm 2a is close behind. Arms 1 and 3 have nearly the same RMS (11.3
vs 11.2 µm) but in opposite directions: arm 1 is 18.7 µm *under*, arm 3 is 10.7
µm *over*, and arm 3 is already over by 10.6 µm at their second data point, so
its **shape** is wrong and not just its amplitude.

## 2. Transport tuning asymptotes almost exactly at their neck value

This was the reason for running ×30 and ×100 rather than one arm, and the
prediction held. ×30 → ×100 is 3.3× more `D_v` for only **+2.66 µm** of neck,
so the attachment-limited ceiling is near 63–65 µm — within ~1 µm of Molaro's
64.78. Extrapolating from the ×10 arm alone (`~D^0.3`) had predicted a match at
≈×28; the actual ×30 result is 60.26 µm, so that extrapolation was mildly
optimistic and saturation was already setting in, exactly as `L* = D_v·β_HK`
(41.8 µm at ×30, 139 µm at ×100) said it would.

**The consequence is a real physical statement:** at `alpha_c = 0.1` the
untuned model's neck shortfall is entirely a *transport* shortfall. Remove the
transport limitation and the neck lands on the data without touching the
kinetics.

## 3. But no arm fits both observables, and that is structural

Panel C of the figure is the result: nothing sits at the origin. The two
observables are fed by the same vapour and pull in opposite directions —

- neck growth is driven by the **internal** curvature difference between the
  neck and the grain surfaces;
- grain recession is driven by the **external** wall undersaturation.

Raising `D_v` strengthens the neck as a vapour sink relative to the wall, so
the arms that fit the neck lose the shrinkage: 2a and 2b recede only 0.41 % and
0.24 % against the −2.93 % target, a factor 7 and 12 short.

## 4. Part of the shrinkage miss was a measurement-window bug

`run_batch_measure.sh` anchored the neck at `t*` but reported `dR_large` over
the **whole run**. `R_large(t)` is linear to R² = 1.0000 at fixed humidity, so
that overstates the shrinkage by exactly 120/78 = 1.54×, and a humidity fitted
against it lands 1.54× too saturated.

Arm 1 shows it plainly: **−2.92 % full-run** (looks like a bullseye against
their −2.93 %) versus **−1.90 % over the anchored window** (35 % short). The
h = 0.99816 in that arm came from a Newton step on the earlier h = 0.99797 run's
full-run −3.218 %, so it inherited the same 1.54× error.

Fixed in `run_batch_measure.sh` (2026-09-08): `dR_large_pct` / `dR_small_pct`
are now the anchored window, with `dR_large_fullrun_pct` carried alongside so
the two can never be confused again.

## 5. Corrected humidities

`R_large(t)` being linear makes the Newton step exact:
`(1−h)_next = (1−h)_now · (−2.93) / dR_measured`, offset from
`h_eq = 1.0000234`.

| arm | 1−h now | ΔR got | factor | 1−h next | **h next** |
|---|---|---|---|---|---|
| 1 untuned | 1.863e-03 | −1.90 % | 1.54 | 2.875e-03 | **0.99715** |
| 2a `D_v` ×30 | 1.034e-04 | −0.41 % | 7.17 | 7.414e-04 | **0.99928** |
| 2b `D_v` ×100 | 6.340e-05 | −0.24 % | 12.44 | 7.884e-04 | **0.99923** |
| 3 `M_0`×5 | 7.383e-03 | −4.78 % | 0.61 | 4.522e-03 | 0.99550 |

**2a and 2b converge on nearly the same corrected wall** (0.99928 vs 0.99923).
That matters: it says the high-`D_v` regime has one consistent wall value rather
than a per-arm fudge, so the effective-transport framing does not collapse into
two independent fits.

Note the derived series-resistance formula was off by 7–12× at high `D_v` —
much more than the 10 % it missed by at nominal `D_v`. It assumes the wall is
the only sink; once `D_v` is large the neck competes for the same vapour and
the grains recede far less than the formula predicts. The formula is a starting
point at high `D_v`, not a calibration.

## 6. Arm 3's neck is not vapour-fed

Arms 1 and 3 share the same nominal `D_v`. Arm 3 has 4× the wall
undersaturation and 2.5× the grain recession, yet its neck reaches 75.5 µm
against arm 1's 46.1. Under vapour-only transport a more undersaturated wall
*competes* with the neck and should slow it — which is exactly what was
measured before (45.71 µm undersaturated vs 56.26 µm saturated at nominal
`D_v`). Arm 3 inverts that ordering.

That is the signature of neck filling by Allen–Cahn interfacial relaxation
rather than by vapour, which is what `mob_scale 5` with `alph_scale 0.01` sets
up: AC curvature motion is deliberately not coupled to vapour
(`docs/model_description.md` §3.4), so the ice arriving in the neck does not
have to come from anywhere. Combined with the worst shape of the four, this arm
should be dropped rather than re-tuned.

## 7. What to do next

**One more −20 °C round before −5 °C.** Two arms, `D_v` ×30 at h = 0.99928 and
`D_v` ×100 at h = 0.99923, to answer the only question left: does the neck
survive a 7–12× stronger wall undersaturation?

It will cost some neck — an undersaturated wall demonstrably slows it. Scaling
the measured sensitivity at nominal `D_v` (≈ −1840 µm per unit `1−h`, from the
h = 0.99797 / h = 0.99816 pair) by the 3.65× higher recession-per-`(1−h)`
these arms actually achieved gives roughly **−5 µm**, landing near 58 µm.
That estimate carries wide error bars — it is a two-point slope extrapolated
across a factor 12 — which is precisely why it is worth one run rather than
more arithmetic.

- If the neck holds near 58–63 µm at the corrected wall, one parameter set
  reproduces both observables to ~10 %, and −5 °C is worth running.
- If it collapses back toward ~50 µm, the vapour-only model cannot do both at
  once, and *that* is the honest result — it sharpens
  `docs/molaro_validation_synthesis.md` §4 from "reproduces ~50 % of the rate"
  into a specific, quantified trade-off.

Either way −5 °C needs a **new geometry file**: `eps` is temperature-dependent
and baked into the element count, and every Molaro geometry carries
`-eps_valid_temp -20`, which will refuse a `-temp -5.0` run.
