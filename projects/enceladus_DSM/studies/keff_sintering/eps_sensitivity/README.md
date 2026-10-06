# Interface-width sensitivity (eps ×2) — result 2026-10-06

Does the k_eff trajectory depend on the diffuse-interface width? φ 0.325,
−20 °C, seeds 1701–1703, run at eps = 1 µm (production) and eps = 2 µm
(`batch_eps2.txt`, jobs in `CAMP/…eps2.00um…`; $0.45–0.56 each) with the same
experiment file, so dtmax is identical in seconds. Driver `compare_eps.py`;
outputs `eps_sensitivity.png`, `eps_sensitivity.csv`.

## Result: the late-time rate is eps-independent; the early transient is not

Rise k(30 d)/k(t_b) − 1, paired per packing (eps 1 µm / eps 2 µm):

| baseline t_b | seed 1701 | seed 1702 | seed 1703 | mean difference |
|---|---|---|---|---|
| 1 d | 25.1 / 31.9 | 28.4 / 39.5 | 25.6 / 35.7 | **+9.3 points** |
| 2 d | 18.3 / 20.0 | 20.5 / 24.9 | 18.9 / 21.7 | +2.9 |
| 3 d | 14.6 / 15.8 | 16.9 / 19.1 | 15.3 / 16.6 | +1.6 |
| 5 d | 11.0 / 11.5 | 12.7 / 13.8 | 11.7 / 12.0 | +0.6 |
| 8 d | 8.1 / 8.0 | 9.4 / 9.6 | 8.8 / 8.7 | **0.0** |
| 12 d | 5.7 / 5.5 | 6.7 / 6.2 | 6.3 / 6.1 | −0.3 |

- From day 8 on the two widths grow k at the same rate, and the absolute k_iso
  differs by a constant +1.5% (two packings) to +3.6% (one).
- Before that, the eps = 2 µm run lags and then overshoots. At t = 0 it starts
  11% more conductive (the wider band bridges more of each tangent contact)
  with 8% less SSA (more throats sealed), and its relaxation is slower:
  τ_sub ∝ eps², so 11 τ_sub is 4 d at 2 µm against 1 d at 1 µm.
- The seed scatter of the rise is ~2.3 points, so the day-1 difference is
  4σ and the day-8 difference is zero.

## What it means

The eps-dependence is confined to the early transient. That is the regime the
campaign already flags as unresolved: a neck is only represented above
r/R = √(12·eps/R) (0.49 at 1 µm, 0.69 at 2 µm), so early growth, which happens
in necks below that floor, depends on eps until the necks outgrow it.

**For the manuscript:**
- State the late-time rate and the absolute k_iso (within 2–4%) as robust to
  the interface width.
- Do **not** quote "rise from day 1" as an eps-independent number. Either
  measure rises from a later baseline, or say explicitly that the first days
  include a width-dependent relaxation.
- k(SSA) is also eps-sensitive at early times, because SSA itself is (sealed
  throats).

**Open:** this test brackets the production width from above only. Whether the
1 µm early curve is itself converged needs one run at eps = 0.5 µm (about $40;
96 M unknowns). The scaling suggests the 1 µm transient ends near 2 d
(8 d ÷ 4), i.e. after the 11 τ_sub baseline.
