# k_eff domain-size convergence (supplement) — result 2026-10-04

φ 0.325, −20 °C, ungated packings (seam and percolation gates only):
L/R 20 / 30 / 40 / 56 / 80 with 8 / 6 / 5 / 4 / 3 seeds. The L/R 40 point is
the production set (1701–1705). Plus the 5 gated pre-B L/R 40 packings.
Runs: `CAMP/` (batch_rve, 3a, 3b). Driver `analyze_rve.py`; outputs
`rve_convergence.png`, `rve_convergence.csv`.

## Absolute k_iso: L/R 40 is within the standard error of L/R 80

| state | L/R 40 gap to L/R 80 | standard error |
|---|---|---|
| 11 τ_sub | +1.3% | ±4.0% |
| 100 τ_sub | +0.6% | ±4.6% |
| 330 τ_sub | +0.8% | ±5.0% |

With only 3–4 seeds at the largest sizes, this rules out gaps larger than ~8%,
not smaller ones. L/R 30 sits +13% high (4 of its 6 packings start at
k_iso 0.47–0.54 against ~0.385 elsewhere, and it has the largest voids,
1.71 R). Nothing in its setup differs, so this is read as realization
statistics at a small size. The seed CV falls from 11.6% (L/R 20) to ~5%
(L/R 56–80), roughly as 1/L.

## The sintering rise (the manuscript's quantity) is size-independent

k(330 τ)/k(11 τ) − 1, seed mean ± sd (SE):

| L/R | n | rise |
|---|---|---|
| 20 | 8 | 27.7 ± 3.6% (1.3) |
| 30 | 6 | 26.9 ± 2.4% (1.0) |
| **40** | 5 | **27.2 ± 2.3% (1.0)** |
| 56 | 4 | 30.7 ± 2.3% (1.1) |
| 80 | 3 | 27.9 ± 2.2% (1.3) |
| gated 40 | 5 | 28.8 ± 2.3% (1.0) |

L/R 40 vs 80: −0.7 ± 1.6 points. Even the L/R 30 offset in absolute k
disappears in the rise. **Gate effect:** the gated packings rise +1.6 ± 1.4
points more than the ungated ones at the same size, which is not
significant. Option B was a methods choice and did not change the result.
