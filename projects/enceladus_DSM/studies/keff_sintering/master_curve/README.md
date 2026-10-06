# The closure k_iso = F(θ, connectivity) — first pass, 2026-10-06

θ = t / τ_sub(T) is the sintering age. `fit_master.py` fits every production
run (25 packings, 70 runs, all five temperatures) for θ ≥ 30, i.e. after the
early transient that depends on the interface width
(`../eps_sensitivity/README.md`). Outputs: `master_curve.png`,
`master_packings.csv`, `master_fit.txt`.

## What the data give

- **Form.** A power law in age fits better than a logarithm (rms 0.74% against
  0.96% of k/k_ref per packing; the curves bend upward on a log axis):

      k_iso(θ) = k_ref · (θ / 30)^n

- **Growth exponent.** n = 0.066 ± 0.005 (sd over the 15 packings with
  φ ≤ 0.375); the equivalent log-law rate is a = 0.166 ± 0.014 per decade.
  Per porosity: n = 0.068, 0.068, 0.064, 0.059, 0.053 for φ = 0.275 … 0.475.
- **Connectivity does not set the growth rate while the solid is well
  connected.** Within φ ≤ 0.375 the correlation of the rate with contacts per
  grain at the band, z, is +0.17. The rate falls only for the poorly connected
  packings (z ≲ 3, φ ≥ 0.425), where the scatter also doubles.
- **Level.** k_ref/k_ice = 1.37·exp(−5.0·φ) for φ ≤ 0.375, with 6.6%
  packing-to-packing scatter. z explains little of that scatter at fixed φ
  (correlation +0.34).

## F, first pass (θ ≥ 30, φ ≤ 0.375, 2D, air-filled pores)

    k_iso = 1.37 · k_ice · exp(−5.0 φ) · (θ/30)^0.066

It predicts the absolute k_iso of every sample of every run to 7.4% rms, which
is the packing scatter. With a packing's own k_ref the error is under 1%.

## Still to do

- Replace the level's φ by a connectivity measure that carries to 3D. z at the
  band does not do it; test the spanning-cluster fraction, or distance from
  the connectivity limit.
- Add the seeds 4–5 runs (3f) and the six redone 3a runs.
- Refit against SSA as the state variable (k ∝ SSA^−0.8) and decide which form
  the paper leads with.
- Give τ_sub(T, R, α_c) in closed form next to F, with its range of validity
  (attachment-limited: R ≪ L*).
- Repeat on the vacuum-pore k_eff once it exists.
