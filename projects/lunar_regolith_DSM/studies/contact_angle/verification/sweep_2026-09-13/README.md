# Contact-angle sweep — batch 2026-09-13 (corrected)

Rerun of 2026-09-12 with the thin-interface correction ON, φ clamped in the
wall term, and `dtmax = τ_sub/2`. All 11 runs completed, **61/61 distinct
snapshots each**, zero bounds violations, zero divergences.

## Verdict

**The angle sweep passes. The ε ladder is complete for the first time. The
sessile runs are still not converged and should not be quoted.**

## 1. Angle sweep — ✅ PASS

Channel, ε/R = 1/50. θ_inf from a free three-parameter fit.

| θ_Young | θ(t_end) | θ_inf | error | τ [d] | drift across fit windows |
|---|---|---|---|---|---|
| 30° | 33.160 | 30.127 ± 0.020 | **+0.127** | 18.7 | 0.397 |
| 60° | 61.186 | 60.046 ± 0.004 | **+0.046** | 15.3 | 0.068 |
| 90° | 90.617 | 90.007 ± 0.001 | **+0.007** | 14.7 | 0.006 |
| 120° | 120.268 | 120.000 ± 0.000 | **+0.000** | 16.2 | 0.007 |
| 150° | 149.004 | 150.024 ± 0.001 | **+0.024** | 21.1 | 0.023 |

**Max 0.127°, RMS 0.062°** against a 2° criterion. The four per-snapshot
estimates agree below printed precision, and θ_inf is stable across fit windows
(≤ 0.4°), so the extrapolation is sound.

### The correction slightly increased the error

| θ_Young | 30° | 60° | 90° | 120° | 150° | RMS |
|---|---|---|---|---|---|---|
| correction OFF (09-12) | +0.079 | +0.035 | +0.007 | +0.001 | +0.020 | 0.040 |
| correction ON (09-13) | +0.127 | +0.046 | +0.007 | +0.000 | +0.024 | **0.062** |

Systematic but tiny — 0.06° RMS is 30× inside tolerance, and the case for the
correction was never that it improves *this* number. It addresses a spurious
surface-diffusion operator that, by `gt_deficit.tex`'s own argument, this
measurement cannot see. Worth recording rather than hiding.

## 2. ε convergence — ✅ consistent with convergence

| ε/R | θ_inf | error | extrapolation spread |
|---|---|---|---|
| 1/25 | 60.080 | +0.080 | ±0.115 |
| 1/50 | 60.042 | +0.042 | ±0.068 |
| 1/100 | 60.032 | +0.032 | ±0.071 |

The 1/100 point ran cleanly — 61/61 distinct snapshots, φ_min = −1.5e-22 —
after stalling twice in the previous batch. The clamp fixed it.

**Correction to an earlier reading of this table.** It first looked like the
error was flattening toward a ~0.035° floor (1/25→1/50 falling 2.4×, 1/50→1/100
only 1.24×). That over-read three points. Quantifying each θ_inf by refitting
over seven fit windows, the 1/50 and 1/100 errors differ by **0.011°** while
each carries an extrapolation spread of **~0.07°** — the difference is smaller
than the uncertainty in either point. There is no evidence for a floor; the data
are consistent with the expected convergence, the measurement simply cannot
resolve below ~0.05° at this run length.

Resolving it would need longer runs (to tighten the extrapolation), not finer ε.

## 3. Sessile drop — ❌ still not converged, do not quote

| θ_Young | θ(t_end) | θ_inf | error | τ [d] | drift | slope at t_end |
|---|---|---|---|---|---|---|
| 30° | 39.401 | 34.0 ± 2.1 | +4.0 | 101.3 | 2.40 | +0.028 °/d |
| 60° | 64.194 | 62.4 ± 0.6 | +2.4 | 76.0 | 0.81 | +0.000 °/d |
| 120° | 115.861 | 117.6 ± 0.7 | −2.4 | 79.6 | 1.03 | +0.044 °/d |
| 150° | 137.994 | 143.0 ± 1.5 | −7.0 | 85.8 | 2.96 | −0.028 °/d |

**τ is 76–101 days, not the 30–45 previously estimated.** That earlier estimate
came from fitting an exponential to a trajectory only 1.4–2.2 τ long, which
biases τ low — so the 174-day `t_final` chosen from it is again only 1.7–2.3 τ.
Window drift of 0.8–3.0° (channel: ≤0.4°) confirms the fit has not settled.

Reaching 4.5 τ needs `t_final ≈ 3.9e7 s` (450 days), and since the current τ is
itself measured from ~2 τ of data it may still be biased low. See "What to do"
below before spending that.

## 4. Mass conservation — ✅

Exactly constant to printed precision in every channel run; ≤ 2.1e-4 % in the
sessile runs.

## What to do about the sessile case

The sessile test exists to check the angle is not an artifact of confinement.
Three options, in order of cost:

1. **Drop it.** The channel result stands on its own: five angles within
   0.127°, a converging ε ladder, and four independent estimates per snapshot.
   Confinement would have to conspire across all five angles to fake that.
2. **Shrink the box.** Relaxation is vapour-diffusion-limited over the domain,
   and the sessile box is 3× taller than the channel (300 vs 100 µm). Reducing
   Ly to ~2R should cut τ substantially at the cost of some ceiling influence —
   cheaper than a 450-day run and probably decisive enough.
3. **Run it to 450 days** (`-t_final 3.9e7`). ~4× the current cost, and with the
   risk that τ is still underestimated.

Option 2 is the recommendation: it tests the same thing for far less.
