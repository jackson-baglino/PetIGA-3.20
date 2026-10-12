# Does the ice span the periodic cell? (production packings)

`percolation.py`, run 2026-10-12 on the 25 production packings: the four
snapshots stored locally for each -20 C run (0, 10, 20, 30 d) and the final
snapshot of each -5 C run (the most sintered state, t/tau_sub = 1256).
No simulation; about three minutes.

A phase "wraps" in x when one of its clusters is connected to its own periodic
image one cell to the right, so heat (ice) or vapor (pore) can cross the cell
without leaving the phase. Threshold phi = 0.5; 8-neighbour ice, 4-neighbour
pore.

## Result

- **The ice wraps in both directions in every packing at every porosity**, from
  10 d on, and at 30 d at -5 C. At t = 0 the one exception is porosity 0.475,
  seed 2002, which does not wrap in x until the first necks form (it does by
  10 d). 92-99 % of the ice belongs to the wrapping cluster.
- **The pore never wraps** in both directions, at any porosity or time (one
  packing in one direction at t = 0). The pore space is a set of closed
  pockets.

So the statement that the ice "remains disconnected in the horizontal
direction" at porosity 0.425 and 0.475 is NOT true of these packings. The
anisotropy there (k_xx/k_yy = 0.79 and 0.64) comes from a spanning ice
network that conducts less well across than along the deposition direction,
not from a broken one. The test of 2026-09-15 that found the ice failing to
span by porosity 0.45 was on an earlier, differently built set of packings
(`studies/packing_design/`); it does not carry over.

What the test does not measure: how good the connection is. At t = 0 grains
that only touch are joined at phi = 0.5 by the diffuse band, so "wraps" can
mean "through point contacts".

## What does change at high porosity: the horizontal direction (2026-10-12)

`directional_exponents.py` fits the Figure 7 power law to k_xx and k_yy
separately (y is the direction of deposition).

| porosity | p of k_xx | p of k_yy | p of the mean | k_xx/k_yy |
|---|---|---|---|---|
| 0.275 | -0.79 +- 0.02 | -0.79 +- 0.04 | -0.79 | 1.02 +- 0.13 |
| 0.325 | -0.80 +- 0.07 | -0.84 +- 0.05 | -0.82 | 1.03 +- 0.03 |
| 0.375 | -0.70 +- 0.10 | -0.86 +- 0.11 | -0.78 | 0.93 +- 0.09 |
| 0.425 | -0.60 +- 0.12 | -0.83 +- 0.13 | -0.73 | 0.82 +- 0.15 |
| 0.475 | -0.48 +- 0.11 | -0.73 +- 0.08 | -0.63 | 0.66 +- 0.18 |

- Along the deposition direction the law holds to porosity 0.425 (p between
  -0.79 and -0.86) and weakens only at 0.475. Across it the exponent falls
  steadily from 0.375 on. The weaker exponent of the isotropic mean in
  Figure 7(d) is the horizontal direction pulling it down.
- The packings become anisotropic in the same order: k_xx/k_yy 1.02, 1.03,
  0.93, 0.82, 0.66, with k_xx < k_yy in 4 of 5 packings at 0.425 and 5 of 5
  at 0.475. A systematic direction is a property of how the packing is
  built (grains dropped under gravity rest on those below), not of sampling.
- `weakest_section.py`: the least ice on any cross-section falls with
  porosity about equally in both directions (0.52 to 0.30 in x, 0.53 to 0.32
  in y), and the number of ice bridges on a cut hardly changes. So the
  anisotropy is not visible in how much ice crosses a cut; it is in how the
  contacts are arranged.
- Domain size: the convergence study (`../rve_convergence/`) was run at
  porosity 0.325 only. Whether 2 mm is a representative volume at 0.425 and
  0.475 is untested. The spread of k_xx/k_yy between packings is larger
  there (0.15, 0.18) but it is 0.13 at 0.275 too, so the data neither show
  nor rule out a porosity-dependent representative volume.

Files: `directional_exponents.txt`, `weakest_section.txt`/`.csv`, `percolation.csv` (every run and snapshot), `percolation.txt` (counts
by porosity and instant), `percolation.png`.
