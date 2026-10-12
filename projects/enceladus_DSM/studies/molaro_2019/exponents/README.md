# Neck-growth exponents of the two Molaro runs

`fit_molaro_exponents.py` fits w = C (t + t0)^a, t0 free, to the two final
manuscript runs (alpha_c = 0.1, nominal D_v, thin-interface term on) over the
samples after the relaxation period (11 tau_sub), exactly as
`../demmenie/analyze_demmenie.py` does for the saturated run. Run 2026-10-12.

| | exponent a | by fit window |
|---|---|---|
| model, -20 C (h = 0.99715) | 0.135 +- 0.001 | 0.135-0.147 |
| model, -5 C (h = 0.99674) | 0.147 +- 0.001 | 0.151 falling to 0.107 for windows starting after 40 min |
| model, saturated pair (-20 C) | 0.212 +- 0.001 | 0.212-0.234 |
| Molaro's data, -20 C (9 points) | 0.23 +- 0.03 | |
| Molaro's data, -5 C (3 points) | not determined | |

- The relaxation period is 15 s (-20 C) and 4 s (-5 C) here, so it removes
  only the first two samples; restricting the fit to the part of the run that
  overlaps the experiment (from the model time at Molaro's first width, for
  78 or 48 min) gives the same exponents, 0.135 and 0.147.
- Molaro's -5 C series has three points, as many as the fit has parameters,
  so it has no exponent. Their -20 C exponent moves with the fit protocol
  (0.20-0.24 in earlier fits); +- 0.03 is the formal error of this one.
- Both unsaturated model runs sit at about 1/7, the exponent Demmenie et al.
  (2025) report for undersaturated conditions, and the saturated run is
  higher: the model follows the saturation trend.

- At -5 C the exponent falls as the window starts later (0.15 to 0.11): the
  grains are shrinking fastest there, which slows the neck. The red line in
  that panel of `exponents.png` is the three-point fit and means nothing.

Files: `exponents.txt` (full numbers), `exponents.png` (exponent against the
start of the fit window).
