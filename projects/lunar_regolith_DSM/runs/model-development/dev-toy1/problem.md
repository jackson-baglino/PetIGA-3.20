# Problem brief (DRY RUN ONLY): 1D heat conduction with a source

This toy problem exists only to test the plumbing of the model line (S3:Build). It is not a result and
must not be cited. Keep every stage short.

## The system

A rod 0 ≤ x ≤ L with constant conductivity k, density ρ and heat capacity c_p, a uniform volumetric heat
source q, and fixed temperatures T(0) = T(L) = T_0. Find T(x, t) from T(x, 0) = T_0.

## Decided by the user

- Notation: T temperature, κ = k/(ρ c_p) the diffusivity, q the source.
- The implementation is a small Python module, not the C solver: `toy/heat1d.py` in the feature
  worktree, function `solve(n, dt, t_end, params) -> (x, T)`, second-order finite differences in space,
  backward Euler in time, NumPy only.
- Comparison test: the steady state must equal the analytic parabola T_0 + q x (L − x) / (2 k).

## What the line must produce (keep each short)

M1: the scope (mechanisms in: conduction and source; out: everything else), the scales and the one
dimensionless group. M2: a choice between finite differences and finite elements (one paragraph each).
M3: the strong form in a ledger (no papers need to be cited; derived steps only), nondimensionalized.
M4: the discrete equations, the residual and the Jacobian (tridiagonal), and a manufactured solution.
M5: tests (steady state against the parabola; manufactured solution with order 2 in h and 1 in dt).
M6: the module. M7: mutations. M8: skip (record "not applicable to the toy"). M9: a one-page note.
