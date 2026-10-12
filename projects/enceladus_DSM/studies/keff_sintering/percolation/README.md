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
  seed 2005, which does not wrap in x until the first necks form (it does by
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

Files: `percolation.csv` (every run and snapshot), `percolation.txt` (counts
by porosity and instant), `percolation.png`.
