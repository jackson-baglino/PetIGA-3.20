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

## Follow-up (2026-09-28): fixing the seam

**The proposed 3×3 superdomain was not adopted.** Copying every grain into the
8 surrounding tiles makes the deposition periodic in y by construction. But
gravity deposition has a direction, so the stitch doesn't disappear; it moves:

- The images of the FIRST grains deposited sit at y + Ly from the start, as a
  ceiling.
- The LAST grains in the centre tile must fill the gap between the rising bed
  and that fixed ceiling. Gravity cannot place them there cleanly, and the gap
  they leave is the same contact-poor layer, at the join between the last and
  the first layers.
- The first layer also needs a floor, which puts a flat layer right at that
  join.

x already works this way: deposition wraps sideways, and that is why the
x boundary shows no defect.

**Adopted: a y-seam acceptance gate** (`--min-seam-contact`, default 0.76, in
`preprocess/generate_packing_gravity.py`, measured by
`packing_lib.seam_contact_ratio`).

- **Threshold.** 0.76 is the 10th percentile of the same measure over ordinary
  interior bands (2059 bands, 29 packings). A seam passes when it is no worse
  connected than 9 interior bands in 10.
- **Test.** 12 packings: φ 0.30 / 0.325 / 0.35 / 0.40 × seeds 1–3, L/R 40.
  - 11 accepted, with seams 0.76–0.99 (previously 0.32–1.06, median 0.72).
  - Bulk fabric is unchanged: F_xx/F_yy 1.088 ± 0.028 against 1.099 ± 0.053
    for the earlier packings. The gate selects on the join, not the interior.
  - Attempts per accepted packing: 2–57. Most rejections also fail the
    pre-existing void gate; a poorly meshed seam leaves a void along the join,
    so the two gates see the same defect.
  - `--max-tries` raised from 24 to 128. A full attempt costs ~3.5 s at L/R 40.

**Seeds are now salted with the porosity** (`_stream`).

- Seeding from the number alone gave seed N the same drop positions and radii
  at every porosity. The old porosity sweep compared one realization under
  different rolling budgets: F_xy was positive for seeds 1–3 at every
  φ ≥ 0.40.
- With the salt, F_xy of the 11 test packings is −0.0001 ± 0.008 (5 of 11
  positive).
- Packings made before this change are not reproduced by the new code; their
  `grains.dat` is the record. `metadata.json` now carries `rng_stream`.

## 2026-10-04 — the anisotropy reverses at high porosity (production runs)

The 25 production packings at −20 °C (3a + 3b) give seed-mean k_xx/k_yy at 30 d
of **1.02, 1.02, 0.91, 0.79, 0.64** for φ = 0.275 … 0.475 (sd 0.13, 0.03, 0.10,
0.16, 0.17). Above φ ≈ 0.35 the deposition direction (y) conducts better,
although the contact fabric leans the other way (F_xx/F_yy = 1.04–1.12). The
"k_xx > k_yy from the rolling rule" finding above holds only for φ ≤ 0.325.

`chord_anisotropy.py` tested the obvious explanation: longer continuous ice
along y. It isn't that. The mean ice chord ratio L_x/L_y is 1.006–1.016 at
every φ (`chord_anisotropy.csv`, `.png`), essentially isotropic and slightly x-leaning.
Per packing it correlates with k_xx/k_yy (r = 0.81 at t = 0), but mostly because
both trend with φ. A ~1% chord change cannot carry a 36% k change.

What the data do support: the reversal and its seed scatter both grow as the
solid approaches 2D percolation (≈ 0.40–0.45), where conduction runs through a
few backbone paths and small directional biases in connectivity are strongly
amplified. The φ 0.475 snapshots show vertical ice columns. Sintering
strengthens the effect: 0.71 → 0.64 at φ 0.475 from t = 0 to 30 d.

**Open:** a directional-connectivity measure of the backbone, e.g. the
spanning-cluster mass, or a max-flow along x vs y on the raster, would test the
amplification idea directly. Until then, report the measured k_xx/k_yy per φ
as a seed mean with its scatter, and do not attribute it to the contact fabric.
