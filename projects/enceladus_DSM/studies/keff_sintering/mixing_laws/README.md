# Volume-weighted geometric mean against the simulated k_eff

`geometric_mean.py`, 2026-10-13. Prompted by a paper on freezing in sand
molds whose Eq. 2 the user described as "the product of the thermal
conductivities weighted by their volume". The PDF is in iCloud and could not
be opened from the session, so the rule tested is the one that description
names:

    k_geo = k_ice^(1 - phi) * k_air^phi ,   k_geo/k_ice = exp(-4.74 phi)

## Result

- **Porosity dependence: reproduced for porosity <= 0.375.** At every
  sintering age the simulated k_eff/k_ice falls as exp(-5.0 phi)
  (-4.99, -5.08, -4.99 at theta = 30, 331, 1256) against -4.74 for the
  geometric mean. The ratio k_eff/k_geo is the same at the three porosities
  to within 3 %.
- **Level: not reproduced, because it changes with sintering.** k_eff/k_geo
  is 0.80 at the opening frame, 1.26 at theta = 30, 1.47 at 331 and 1.66 at
  1256. A rule that knows only volume fractions cannot follow that.
- **Porosity 0.425 and 0.475 fall away** (ratio 1.10 and 1.03 at theta = 30),
  the anisotropic packings of `../percolation/`.
- So for phi <= 0.375, k_eff = C * k_ice^(1-phi) * k_air^phi with C a
  function of the sintering state only; with the Figure 7 law,
  C = 1.26 (SSA/SSA_r)^-0.80.

Caveat: one conductivity contrast (k_air/k_ice = 0.0087) and two dimensions.
The match of the slope to ln(k_air/k_ice) may be particular to this contrast;
it would take conductivity replays of stored snapshots at another k_air to
tell.

Files: `geometric_mean.txt`, `geometric_mean.png`.
