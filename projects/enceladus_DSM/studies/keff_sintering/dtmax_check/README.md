# dtmax 2·τ_sub vs 1.09·τ_sub on a packing — PASS

The lunar ladder that moved dtmax from 1.09·τ_sub to 2·τ_sub (2026-10-02) had
no packing. This compares the gated seed301 packing (φ 0.325, −20 °C) run at
both settings, with all other options the same apart from the k_eff cadence.

- 1.09·τ_sub: the first 3a shakedown,
  `GRP/batch_2026-09-28__13.50.25_keff_3a_shakedown/…seed301…T-20…__kf5`
  (local `LOC/GPS/keff_b3a_shakedown/`). 367 steps.
- 2·τ_sub: the 3a rerun, `CAMP/…seed301…T-20__snow_T-20_h1.00_30d`
  (job 3774294). 240 steps.

Differences for t ≥ 11 τ_sub (1 d), against a pass limit of 0.5% (seed
scatter is ~3%):

| quantity | max \|diff\| | mean |
|---|---|---|
| SSA(t) | 0.21% | −0.08% |
| k_iso(t) | 0.13% | +0.04% |
| k_iso at matched SSA | 0.08% | −0.03% |

So 2·τ_sub costs 35% fewer steps and changes the results by about 0.1%,
far inside the seed scatter. Every 2·τ_sub run stands.

`compare_dtmax.py` (driver), `dtmax_check.png`, `result.txt`.
