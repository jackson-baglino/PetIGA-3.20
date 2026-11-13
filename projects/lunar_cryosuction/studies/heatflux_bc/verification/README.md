# Thermal face conditions: verification

Checks for the two opt-in thermal boundary conditions added in commit
`4d7e0a61` (branch `heatflux-bc`, FREEZE.md ledger):

- `-qflux_faces x0 -qflux q [-qflux_t0 t0]`: prescribed heat flux,
  k dT/dn = q on the named faces (q > 0 heats the domain), on for t >= t0.
- `-Twall_faces y0,y1 -Twall_profile erfc|erfc_steady|linear -Twall_axis a
  -Twall_s0 s0 ...`: prescribed wall temperature T = T_w(s,t),
  s = |x_a - s0| the distance along the wall from the heated end, imposed
  weakly by symmetric Nitsche (beta = C k(phi)/h_n, `-Twall_nitsche`, C = 10).
  - `erfc`: T_init + dT erfc(s / (2 sqrt(kappa (t - t0)))), T_init before t0
    (`-Twall_kappa`, `-Twall_dT`, `-Twall_Tinit` default `-temp`, `-Twall_t0`)
  - `erfc_steady`: the same shape frozen at t - t0 = `-Twall_tref`
  - `linear`: T_init + dT (1 - s/L) for s < L, T_init beyond (`-Twall_length`)

Run (from the project root, with `PETIGA_DIR` pointing at the tree that holds
this folder):

```bash
bash studies/heatflux_bc/verification/verify_heatflux.sh
```

Every simulation goes through `scripts/Studio/run_cryosuction.sh`, one at a
time (the largest, the 155 x 52 channel of W1, is 25k DOF). Inputs are in `inputs/geometry/heatflux/` and
`inputs/experiment/heatflux/`. `check_heatflux.py` reads each run folder
(spline fields evaluated with igakit through `postprocess/pplib.py`; stored
heat integrated exactly from the control values, since
int N_i dx = (U[i+p+1] - U[i])/(p+1)) and prints one CSV row per check.
Results: `verify_heatflux.csv`, full log `verify_heatflux.log`.

Every field in these tests is a uniform phase (phi = 0 or 1 at every control
point), so there is no interface and no phase change: the vapor equation is
uncoupled, and the latent term in the energy budget is identically zero.
The latent part of the budget is therefore not exercised here (see
"Not covered").

## What each check asks, and the result

Numbers from the run of 2026-11-13 on commit `4d7e0a61` (all PASS; the full
rows are in `verify_heatflux.csv`).

| Test | What it checks | Quantity | Measured | Expected | Tolerance (why) |
|---|---|---|---|---|---|
| W1 | Runs without a thermal face are unchanged WITH `-wall_faces` (the contact-angle branch was restructured); merge-base with main, theta-60 channel, 32 steps | `SSA_evo.dat` | identical | identical | exact |
| G1 | Same without `-wall_faces` (`verify_wall_bc.sh`, also G2–G4) | `SSA_evo.dat` | identical | identical | exact |
| J | Thermal boundary block of J (Nitsche T and phi columns) vs central differences, isolated by switching the terms off; interface crosses the walls, dT/dn != 0 | rel. error | 2.8e-8 | 0 | 1e-6 (FD truncation at h = 1e-7, as G2) |
| J | `-snes_test_jacobian`, one step, full system | ‖J−J_fd‖/‖J‖ | 1.3e-4 | — | INFO only: same size as the full-system column of the J check (1.5e-5 to 1.3e-4), in which the isolated boundary block agrees to 2.8e-8, so the difference is attributed to the interior form (inferred, not tested separately) |
| a | Steady flux q = 100 W/m² on x0, T fixed on x1, air: T = T0 + (q/k)(L − x); two meshes (h = 5, 2.5 µm) | max\|T − T_lin\|/(qL/k) | 1.1e-14, 1.1e-14 | 0 | 1e-8. The linear profile lies in the p = 2 space, so the discrete steady state is exact on every mesh; refinement cannot reduce a round-off error (convergence under refinement is shown in c) |
| a | Heat leaving through x1 at steady state | \|Q_out − q Ly\|/(q Ly) | 0 | 0 | 1e-8 |
| b | Flux into a 600 µm slab, air (q = 100) and ice (q = 1000): T(0,t) − T0 = (2q/k) sqrt(αt/π) at 4 times; whole profile with the far face by images | \|ΔT(0) − exact\|/exact; max_x error / ΔT(0) | ≤ 3.4e-5; ≤ 3.4e-5 | 0 | 1e-3 ((h/sqrt(αt))² ≈ 1e-3 at the earliest time) |
| c | erfc wall temperature (κ = α_air) on x0, x1, y0, y1; interior must be the 1D erfc, uniform in y; t = t_f/64 … t_f | max_xy\|T − T_erfc\|/dT | 4.0e-4 (h = 8 µm), 3.2e-4 (h = 4 µm), at t_f/64; ≤ 3.1e-6 from t_f/4 | 0 | 1e-3 |
| c | Order of convergence h → h/2 at the time where a dt/4 run changes the error by < 5 % (t_f/4: 0.2 %) | order | 1.98 | 2 | > 1.5. At t_f/64 and t_f/16 the dt/4 run halves the error: the time integration of the jump at t = 0 dominates, so the order there (0.3–0.4) is reported as INFO |
| c2 | Linear wall temperature on all faces, steady: exact in the space, so Nitsche must be consistent | max_xy\|T − T_lin\|/dT | 7.1e-15 | 0 | 1e-8 (a pure penalty would leave an O(k/(βh)) error) |
| d | Energy, flux runs: ΔE = ρc ∫ΔT dA (exact from control values) vs ξ_T q Ly t | max_t\|ΔE − Q_in(t)\|/Q_in(t_f) | 1.0e-5 (air and ice) | 0 | 1e-4 |
| d | First step of the flux runs | offset/(q Ly dt0) | 0.20 | — | INFO: TSALPHA starts with U_t = 0 while the flux switches on at t = 0; afterwards the offset stays within 3.2e-11 to 3.7e-11 J/m (it moves only when dt changes), against a final input of 4e-6 |
| d | Energy, erfc wall (h = 4 µm): ΔE vs time integral of −k dT/dn on the wall faces (`bflux.dat`, trapezoid) | ΔE / Q_in | 1.3287e-4 / 1.3275e-4 J/m (9e-4) | — | 2e-3 (bflux omits the Nitsche β(T − T_w) exchange; trapezoid on a 1/sqrt(t) flux) |


## Not covered

- Latent heat in the budget: with no interface the term
  xi_T rho_i L d(int phi)/dt is zero. Ice dynamics are many orders slower than
  heat diffusion at these scales, so a test with a moving interface would need
  a run long enough for the thermal problem to be trivially steady.
- `erfc_steady` is not run on its own: it is the `erfc` formula evaluated at
  a fixed time, sharing all code with `erfc` (test c) except the argument.
- `-geom_file` (mapped) meshes: h_n uses dx/du from the geometry map; only
  Cartesian patches are tested.
- Axisymmetric runs: the thermal terms carry the same r weight as the wall
  term; not tested.
