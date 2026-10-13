# Surface diffusion against vapor transport in Molaro et al. (2019)

`sa81_mechanisms.py`, 2026-10-13. Question: does surface diffusion account for
the part of Molaro's -20 C neck growth that our vapor-only model does not
reproduce (59 %)? The script implements the surface-diffusion term (V1) and
the vapor term (V3) of the Swinkels & Ashby model exactly as printed in
Molaro et al. (Table 3, Eqs. 1-4, 6, 17-20, Table 5) and evaluates them.

## Result

- **At the necks they measured (radius 16-32 um, grains of about 87 um,
  -20 C), surface diffusion is 0.7 % falling to 0.04 % of the vapor term**,
  using their own equations and the top of their coefficient range.
- **Their vapor term alone gives the measured rate**: dx/dt = 7.8e-9 falling
  to 1.5e-9 m/s against a measured mean of 3.4e-9 m/s. That term is
  Hertz-Knudsen with alpha = 1 and no gas diffusion (their Appendix A).
- Surface diffusion does dominate, in their model, for small necks: below a
  neck radius of about 4 um at this grain size and temperature. That is the
  early part of their Figure 10, where V1 is the largest current. The two
  largest currents in that figure are V1 and V3, but not at the same time.

## What this means for our comparison

The 59 % our model lacks at -20 C is not surface diffusion in Molaro's own
framework. Their model matches the data because its vapor route is the
fastest conceivable (alpha = 1, vacuum), not because of a second mechanism.
Ours has alpha_c = 0.1 and diffusion through air.

## How far the check goes

- Figure 10 is reproduced in shape, not in value. The vapor current agrees
  (2e-21 rising to 3e-20 m3/s against about 1e-21 rising to 2e-20 in the
  figure). The surface current here starts about 100 times higher than
  theirs and the crossover is later (x/a = 0.2 against an earlier time in the
  figure). The difference is in the coupling to grain-boundary diffusion:
  with their parameters A = db Db/(ds Ds) = 3.5e8, and the approximate
  solution (Eq. 20) raises the driving curvature K2 by a factor of 4000.
  They solve the four equations iteratively.
- So the numbers above OVERSTATE surface diffusion relative to their figure.
  With the coupling off the ratio at the measured necks is 8e-7 to 4e-8.
- An earlier estimate with Kuczynski's simple law gave 1e-5 to 1e-7; it sits
  between the two.
- S (their smoothing function) is set to 1; it multiplies both terms.

Files: `sa81_mechanisms.txt`, `sa81_mechanisms.png`.
