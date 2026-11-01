# F2 (pore-space formulation, 90 degrees, conservative form): equations, checks done, and the implementation plan

> **Provenance (added 2026-11-01).** Copied unchanged below this note from the
> practice-round archive, `scratch/2026-10-21_practice_round/research_archive/Planning/`
> (git-ignored, local only). Other numbered docs it cites (docs 14 and 20, and `Review/memos/`) are still there.
> Scripts under `agentic/math/` and the `agent/math-lane` branch are in tag
> `archive/2026-10-21/agent-math-lane`.

Status: 2026-10-18. The equations below were re-derived by the math lane (two models, algebra-checked) and two properties were confirmed on toy numerical solves. **Nothing has been implemented in the C solver.** Sources for the design: the R8 gate memo (`Review/memos/gate_R8_diffuse_domain.md`, 22 papers) and doc 14 section 8.

## 1. The model
Fields: phi_s (sediment indicator, frozen) gives the pore fraction w = 1 - phi_s. Unknowns: psi (ice fraction within the pore space, 0 to 1), temperature T, vapor density rho_v. Derived: phi_i = w psi, phi_a = w (1 - psi). The sum phi_i + phi_a + phi_s = 1 holds identically.

Energy F[psi] = integral of w [ (3/eps) W(psi) + (3 eps/2) |grad psi|^2 ], with W(psi) = psi^2 (1 - psi)^2 / 2. Conservative weighted gradient flow (the pore measure w is the metric), plus the phase-change source:

  w psi_t = 3 M0 [ eps div(w grad psi) - (w/eps) psi (1 - psi)(1 - 2 psi) ] + w (alpha/rho_i) psi^2 (1 - psi)^2 (rho_v - rho_vs(T))

  vapor:  d/dt (w (1 - psi) rho_v) = div( w (1 - psi) D grad rho_v ) - rho_i w psi_t
  heat:   rho c_p T_t = div( K grad T ) + rho_i L w psi_t   (rho, c_p, K from the mixtures of phi_i = w psi, phi_a, phi_s as in the current code)

Why the conservative form: in the weak form the psi equation has no division by w (the 1/w of the divided form is what makes it singular inside grains); all terms simply carry w. Where w is tiny the equation row vanishes, so a floor w_tilde = w + upsilon (as the literature does, 1e-6 to 1e-8) is used consistently in every term.

At 90 degrees there is no wall term: the weighted gradient energy alone gives the no-flux condition grad(w).grad(psi) = 0 at the grain surface, i.e. a 90 degree contact angle (Yu et al. 2012, Eq. 20 with cos(theta) = 0). The ice-air tension and kinetics are exactly those of the two-phase model (Sigma_i = Sigma_a drops out; M0, eps, alpha are the same numbers).

## 2. What is established (and how)
| Property | How shown | Status |
|---|---|---|
| The equation is the weighted gradient flow of F; dF/dt = - M0 integral of (deltaF/deltapsi)^2 / w <= 0 without the source | lane task c3a, 16/16 checks, Opus replicate agrees | established where w > 0; the adversary notes the identity is undefined where w = 0 (grain interior), handled only through the floor |
| psi = 1 and psi = 0 are exact steady states for any w, linearly stable (well curvature 1) | c3a | established |
| phi_i + phi_a + phi_s = 1 and phi_a >= 0 identically for psi in [0,1], 0 <= w <= 1 | c3a | algebraic fact; that the flow keeps psi in [0,1] is a maximum-principle argument (adversary), not proved in the lane |
| w = 1 gives the two-phase equation | c3a | established |
| Front in a pore of slowly varying cross-section moves at V = -3 M0 eps d(ln w)/dx, i.e. toward smaller w (area minimisation) | c3b (two models agree); 1D solve gives 0.998 of the prediction for four slopes | established to first order |
| F1's overshoot of phi_a on a curved grain scales as 0.18 / R_s (min phi_a = -0.015, -0.022, -0.036 at R_s = 12, 8, 5, in units of eps), Lambda = 0 or 1 alike; F2 gives phi_a >= 0 to round-off | toy 2D solve `agentic/math/numerics/f1_vs_f2.py` (explicit, periodic, no vapor) | direct support for the curvature-residual mechanism (H1) and for F2's construction; this toy shows the **negative** overshoot, not a positive film (no vapor or phase change in it) |
Not established: coupling to vapor and heat near the wall (first integral of the vapor equation with varying w), thin-interface corrections near the wall, behaviour in narrow throats, the effect of the floor on convergence, angles other than 90 degrees.

## 3. Implementation plan (needs your go; nothing is changed until then)
**Base.** The three-phase kernel lives on branch `feat/diag-air-film` (`src/assembly_sed.c`: `Residual_A2`, `Jacobian_A2`; dof = 4 with the frozen sediment as the fourth field). New branch `feat/f2-pore-space` from it, never merged without your approval, behind a runtime option (`-sed_formulation f2`, default off so every existing run is unchanged).
**Residual (new `Residual_A3`, replacing the ice, heat and vapor rows when the option is on).** With w = max(1 - phi_s, 0), w_tilde = w + upsilon (new option `-w_floor`, default 1e-6, same spirit as `-air_lim`), phi_i = w psi, phi_a = w (1 - psi):
- ice row: N w_tilde psi_t + 3 M0 eps w_tilde grad(psi).grad(N) + (3 M0/eps) w_tilde psi(1-psi)(1-2psi) N - pc (alpha/rho_i) w_tilde psi^2 (1-psi)^2 (rho_v - rho_vs) N
- heat row: unchanged form with phi_i = w psi in the mixture properties and the latent term rho_i L w_tilde psi_t
- vapor row: storage and diffusion weight phi_a (floored as now), coupling (xi_v rho_i - rho_v) w_tilde psi_t
- sediment row: unchanged (phi_s_t = 0).
**Jacobian.** New `Jacobian_A3`, derived block by block; verified with `-snes_test_jacobian` (the standard in this project: error around 1e-10), including the columns for the frozen phi_s.
**Bounds.** The psi unknown gets the plain [0,1] bounds already used for the two-phase ice variable; the sediment-aware node-wise bound (`-vi_sed_bound`) is not needed, because the constraint is built in.
**Initial conditions.** psi is given by the ice-air geometry alone (disc, slab, pair of grains), with phi_i = w psi derived; the existing sediment-core initial conditions get an F2 variant that sets psi from the ice shape without carving out the sediment.
**Output and post-processing.** Write psi and also derived phi_i and phi_a; `to_vtk.py` and the film/bubble scripts read the derived fields; the dof-4 output format is unchanged.
**Tests (ladder of doc 14, with the memo's T1-T8):**
0. Unit test: with phi_s = 0 the F2 residual equals the two-phase residual (`Residual_A1`) to round-off.
1. 1D ice against a flat diffuse wall: max phi_a in the ice-sediment band is zero to round-off for all times.
2. 2D ice-air interface meeting a flat wall: angle 90 degrees to better than 2 degrees as the sediment width and eps are varied independently (resolve every interface with at least 3 to 5 elements).
3. Ice cap on a circular grain against the analytic carved-circle equilibrium; no overshoot at R_s down to the smallest radius that failed in F1.
4. Convergence in eps (expect about first order), floor sweep (upsilon 1e-8 to 1e-4, watch the time step), channel-drift test (c3b numbers) with the actual code.
5. The grain pair with sediment cores (your bubble configuration), with a day-long run and a ParaView snapshot saved through `save_review.sh`.
Each rung is run by you; I prepare the commands and the acceptance numbers.
**Risks to watch (from the memo and the derivations):** stiffness and conditioning from the weights (the floor trades accuracy for stability); first-order accuracy of the wall condition when the sediment band is as wide as eps; the vapor equation near the wall (varying w breaks the exact first integral, so a thin-interface correction near the wall is untested); ice-air tension appearing reduced in narrow throats (to be tested, expected to be only the geometric drift); order parameters where w < 0.5 are not physical, so contour psi and mask diagnostics to w >= 0.5.
**Open choices (my default in brackets):** floor form [additive upsilon]; whether psi is also solved inside the grains [yes, with the floor]; weight of the source [w, which leaves pore-space kinetics identical to the two-phase model].

## 7. Revision after the six-iteration audit (supersedes the conflicting lines above)
1. **Floors.** The floor w_tilde = w + upsilon is applied only in the psi row (mass and stiffness terms). The vapor storage coupling and the latent-heat term use the unfloored w, so that d(phi_i)/dt = w psi_t exactly and ice, vapor and heat balance exactly; inside grains the coupling vanishes with w.
2. **Source weight.** The psi-row source keeps the weight w (it is a requirement for the exact 90 degree contact under driving, doc 20 iteration 3).
3. **No clamp** in w = 1 - phi_s; phi_s is in [0,1] by construction of the initial condition, so the Jacobian has no kink.
4. **Bounds.** The [0,1] variational-inequality bounds on psi are required (no discrete maximum principle is guaranteed).
5. **Initial conditions.** psi = 1 where w < 0.01 (grains); restart from an F1 state uses psi = clip(phi_i/w, 0, 1) with a count of nodes needing clipping.
6. **Output.** The formulation flag is written to the output metadata; post-processing and restart abort on a mismatch; derived phi_i = w psi and phi_a = w(1 - psi) are written.
7. **Tests added or changed.** Unit test equality at upsilon = 0 and (1 + upsilon) scaling otherwise; closed-box mass drift (xi_v = 1); PetIGA dt sweep of the discrete energy and the wall slope on the tilted-interface case; h-refinement at fixed eps and eps_s; condition number and row scaling for upsilon in {1e-8, 1e-6, 1e-4}; the ice-cap test is replaced by the curved-wall orthogonality metric (|n_psi . n_w| below 0.03 at long times for eps_s = eps) and a coupled closed-mass test; floor sweep 1e-12 to 1e-4; adaptive-step error contribution by region; a coupled psi-rho_v tapered-pore test against the sharp-interface vapor-limited front.
8. **Documented limitations.** 90 degrees only; conduction at point contacts through the diffuse band (shared with F1); an ice-air front in the wall band is expelled into the grain; wall accuracy O(eps_s).

## 8. Second revision (after the second red-team round; supersedes items 1, 2 and 5 of section 7 where they conflict)
- psi row: the floored weight w_tilde multiplies all four terms (mass, stiffness, well, source). At phi_s = 0 the ice row is exactly (1 + upsilon) times `Residual_A1` (unit test), and the forcing carries the weight of the mass term (exact 90 degree contact under driving).
- Vapor and heat coupling use the unfloored w.
- Initial condition: psi in the band and the grain is the **normal-constant extension** of the adjacent pore value (not 1). A restart from F1 extends psi the same way and converts or zeroes the rate vector.
- Add: control-point bound assertion for phi_s; per-region residual scaling and tolerance weighting of the psi rows by 1/w_tilde; option matrix with an abort on `-vi_sed_bound` + F2; a coupled test of an ice-air front in the wall band (vapor and heat on, production rho_vs/rho_i and xi_v); acceptance by dt-convergence, not by toy numbers; rung 0 is a PetIGA reproduction of the 1D drift and tilted-interface toys.
