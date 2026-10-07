# Neck growth at saturation against t^(1/3)

The Demmenie-conditions run: two equal grains (R = 84.4 um, initial neck radius
14 um) at -20 C in vapour saturated over the grains' own curvature
(h = 1 + 2 d0/R), alpha_c = 1e-3, 100 h, one-grain mirror domain.
Batch file `../batches/demmenie_mirror_T-20.txt`; run on 2026-10-06 (job 4066835).
Data: `LOC/GrainPairSintering/batch_2026-10-06__09.21.04_demmenie_mirror/`.

This run is not compared with Molaro et al. (2019): their chamber was
unsaturated and their grains shrank. It is compared with the
evaporation-condensation law w ~ t^(1/3) and with the exponents Demmenie,
Woutersen & Bonn (2025) measured at saturation, 0.26-0.33.

## Result (2026-10-07)

Fits use Demmenie's form w = C (t + t0)^a, t0 free, over the samples after the
relaxation period (11 tau_sub = 20 min; tau_sub = 108.9 s).

| fit | exponent | rms | largest misfit |
|---|---|---|---|
| free | a = 0.212 +- 0.001 (t0 = 4.8 h) | 0.28 % | 2.0 % |
| one-third, a fixed | 1/3 (t0 = 23.9 h) | 2.2 % | 12.7 % |
| Kuczynski, w^m - w_r^m = K dt | m = 4.84, 1/m = 0.207 | 0.37 % | 1.3 % |

- The free exponent rises with the start of the fit window: 0.212 from 1 h,
  0.222 from 8 h, 0.228 from 20 h, 0.232 from 38 h, 0.234 from 53 h. It levels
  off near 0.23-0.24 and does not approach 1/3.
- The grain radius changes by +0.06 % over 100 h, so the saturation condition
  holds: nothing but the neck drives mass transfer.
- The neck grows from 28 to 63 um, i.e. w/2R from 0.17 to 0.37.

So the model does NOT reproduce t^(1/3) here: a = 0.21-0.23, below Demmenie's
0.26-0.33 and below the envelope-ODE prediction of 0.29-0.30
(`../envelope_ode/README.md`). The cause is open; see `../TODO.md`.

## Files

- `analyze_demmenie.py <run>`: the fits and four diagnostic plots, into
  `<run>/plots/demmenie/`. `postprocess/run_batch_measure.sh` calls it for
  any saturated run and skips the Molaro comparison there.
- `fig_demmenie.py <run>`: the manuscript figure, built from the Molaro figure's
  helpers so the two match (colour bars, sections, instants, panels),
  `../manuscript/Figure4_saturated_neck_growth.{pdf,png}`.
