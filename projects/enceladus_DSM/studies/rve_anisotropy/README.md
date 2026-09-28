# Is the domain an RVE, and where does the anisotropy come from?

Asked 2026-09-28, before the porosity campaign. The worry is a reviewer's:
k_xx and k_yy differ by ~10%, and k_xy reaches 10% of k_iso, so maybe the
domain is too small, or maybe the packing generator builds something strange.

**Existing data only**: no new simulations and no new solves.

- **Geometry.** Every packing on disk (`grains.dat`): 36 with both axes
  periodic, plus 6 x-periodic or open controls.
- **PetIGA k_eff.** L/R = 40:
  - 4 packings: the pilot runs (arith law) and the warm-end batch (tensor law).
  - Seed 2 from its tensor replay, from 16 d on.
- **PetIGA k_eff.** L/R = 64: 1 packing, rev64 (arith law).
- **Finite-volume solver.** k_yy/k_xx at t = 0 for 9 packings at L/R = 10, 20
  and 40, from `studies/packing_design/rev_bias.csv`.

```bash
venv_enceladus/bin/python studies/rve_anisotropy/geometry.py   # geometry.csv, seam_profile.csv
venv_enceladus/bin/python studies/rve_anisotropy/figures.py    # rve_anisotropy.png, k_data.csv
```

![rve](rve_anisotropy.png)

## Method

- **Contact fabric.** F = ⟨n nᵀ⟩ over contacts, where n is the unit vector
  between the centres of two touching grains.
  - A pair counts as a contact when its gap is below the diffuse band,
    9.2·eps with eps = R_ave/50. The solver joins every such pair with ice.
  - In a granular solid, heat crosses from grain to grain through contacts. So
    F predicts the direction of the conductivity anisotropy from geometry
    alone, with no solver involved.
  - F_xx/F_yy > 1 means more contacts lie horizontally, so k_xx > k_yy is
    expected. F_xy ≠ 0 tilts the principal axes.
- **Seam density.** For 256 horizontal lines, count the contacts whose
  centre-to-centre segment crosses the line, per unit length. y = 0 is the
  periodic join. Vertical lines are the control, because the x boundary is not
  a stitched seam.
- **RVE criteria**, applied to each quantity separately:
  - the mean over seeds must not depend on L;
  - the scatter between seeds must shrink with L, as 1/L in 2D;
  - any quantity that symmetry forbids must average to zero over seeds.

## Findings

### 1. The horizontal bias is real and comes from the deposition rule. It is not a small-domain effect.

- F_xx/F_yy is above 1 at every size:

  | L/R | 10 | 20 | 40 | 64 |
  |---|---|---|---|---|
  | F_xx/F_yy, mean | 1.11 | 1.10 | 1.09 | 1.09 |

  The mean is flat with size and the scatter shrinks (sd 0.08–0.10 at L/R
  10–20, 0.027 at 40, 0.023 at 64): textbook RVE behaviour.
- **It persists with no y seam.** The x-periodic and open packings give 1.11
  (sd 0.04), so it isn't produced by the stitching.
- **The mechanism is the rolling step.** A grain lands on the apex of a
  support and rolls toward its side until it finds a second contact. Both
  motions tilt contacts toward the horizontal.
- **k agrees on average.**
  - k_xx/k_yy at t = 0 averages 1.08 at L/R ≥ 40 (sd 0.045).
  - At 30 d the four L/R = 40 packings give 1.10 ± 0.03 (SE), under the tensor
    law.
- **Not designed in.** The generator was built to mimic settling, and
  `studies/packing_design/README.md` expected vertical load paths, i.e.
  k_yy > k_xx. The data say the opposite, and that note should be corrected.
- **Whether real snow does the same** needs checking against measured snow
  anisotropy before any claim is made. Cite only after verifying.

### 2. k_xy is a finite-sample fluctuation, not a material property and not a solver error.

- **Mirror symmetry.** Deposition is mirror-symmetric in x: grains roll in
  whichever direction they lean, and exact ties are broken at random. So the
  average of k_xy over seeds must be zero.
  - At 30 d, L/R = 40: k_xy/k_iso = +0.097, −0.083, +0.015, −0.070. The mean
    is −0.010 ± 0.042 (SE).
- **It shrinks with size.** F_xy scatter falls from 0.040 (L/R 10) to 0.019
  (20), 0.006 (40) and 0.0045 (64). The 1/L prediction from L/R = 40 is 0.0039.
- **Each packing's k_xy has the sign of its own F_xy**, in 5 of 5 PetIGA
  packings. Seed 2's F_xy is essentially zero, so it counts for little. The
  magnitudes do not track, so the fabric predicts the sign only.
- **k_iso is unaffected by the tilt.** k_iso = trace/2, which does not change
  when the axes rotate.

### 3. The y seam is a real generator artifact: a layer with fewer contacts.

- **Density.** Mean contact density across the seam is 0.65 of the interior
  (range 0.33–1.1). It is below the interior in 34 of 36 xy packings.
  - Pilot seeds 1 and 4 are at 0.33 and 0.36. Seeds 2 and 3 are at 0.95 and
    0.90.
- **Cause.** y-periodicity is stitched: the bottom of the cropped window is
  joined to its top, overlaps along the join are pushed apart, and gaps are
  left for the void filler.
- **Upper bound on its effect**, treating the seam as a layer ~2R thick in
  series with vertical heat flow:

  | L/R | seam density | k_yy low by at most | k_iso low by at most |
  |---|---|---|---|
  | 40 | typical (0.6) | 3% | 1.6% |
  | 40 | worst (0.33) | 9% | 4.6% |
  | 64 | 0.68 | 1.4% | 0.7% |

  The effect shrinks as R/L.
- **No clear trace in k at t = 0.** Seam density correlates with k_xx/k_yy at
  only r = −0.21. Seed 4 (seam 0.36) is nearly isotropic. So the fabric
  dominates the anisotropy; the seam is a bounded extra effect on k_yy.
- **Quantifying the seam exactly needs one solve.** k_yy with and without the
  seam, on the same packing. It has not been done, deliberately (existing data
  only).

### 4. Per-seed differences at L/R = 40 are not explained by the fabric.

- The fabric predicts the per-seed anisotropy well at L/R ≤ 20 (r = 0.91), but
  not at L/R ≥ 40 (r = 0.02).
- At production size, seeds differ by about ±5% in k_xx/k_yy for reasons the
  contact count doesn't capture, such as contact quality and long chains.
- That is scatter, which averaging over seeds handles. It is not a bias.

## Decision

- **Domain size: L/R = 40 is adequate for the reported quantities:** k_iso,
  its rise, SSA, and k_eff(SSA).
  - k_iso seed CV is 3.6% at 30 d.
  - The single L/R = 64 run sits inside the L/R = 40 spread, with the same
    relative response.
  - Fabric statistics converge as an RVE should.
- **Report anisotropy as a mean over seeds, never per seed.** k_xx/k_yy is
  1.10 ± 0.03 at 30 d from 4 seeds.
  - Report k_xy as a zero-mean fluctuation that shrinks as 1/L, not as a
    property.
- **Fix the seam before generating the porosity packings.** It's the one
  genuine generator artifact. Cheapest fix: a gate on seam contact density
  (e.g. ≥ 0.8 of the interior), rejecting and re-seeding like the other gates.
  - Pilot seeds 1 and 4 would fail it.
  - The temperature comparison is paired on the same packings, so it is
    unaffected.
- **Open, closable from existing data:** L/R = 64 has one sintered packing.
  If the 2026-09-21 rev64 seeds 2–4 (jobs 3264515–3264517) finished, their
  `k_eff.csv` + `SSA_evo.dat` would test the sintered-state RVE directly.
