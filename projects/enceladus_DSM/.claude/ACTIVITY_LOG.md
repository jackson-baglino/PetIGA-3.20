## 2026-09-02 (final, 2) — 10 in width is now a hard ceiling

- Jackson: nothing may be wider than 10 in, these go on PowerPoint slides.
- The six panels were already exactly 10.00 x 5.80 in with correct pHYs DPI
  metadata (checked the chunk directly, not the figsize). The violator was
  fig0_overview.png at 30 x 11.6 in -- my contact sheet, three panel-widths
  across.
- Overview restacked into ONE column, 10 x 34.8 in at 110 dpi. Two columns
  would give 5 in per panel, which the four-group legend cannot fit. It is a
  review sheet, not a slide asset, and the README now says so.
- Added `--width` (default 10.0) and `--dpi` (default 200), and a `_save()`
  helper that ASSERTS the width ceiling on every figure and prints the physical
  size. A figsize edit can no longer silently produce something unusable.
- Noted in the README: PowerPoint that ignores the PNG DPI tag inserts these at
  width*dpi/96 = 20.8 in. Setting the width box back to 10 in loses nothing
  (the image is then 200 dpi at 10 in); `--dpi 96` gives drop-in sizing if
  preferred.
- Overview suptitle split across two lines -- as one line it was wider than
  10 in and overflowed both ends.

---

## 2026-09-02 (final) — Grouped legends, 10 pt floor, 10 x 5.8 panels

- Panels regenerated at 10 x 5.8 in / 200 dpi (10 x 5 of plot plus a legend
  strip), with no type below 10 pt: FS_TITLE 15, FS_LABEL 13, FS_TICK 11.5,
  FS_LEG 10, FS_NOTE 10. Verified by resolving every fontsize expression in
  the source -- min is exactly 10.
- The flat 8-entry legend strip was a lookup problem. `_legend()` now takes
  [(heading, [labels])] and exploits matplotlib's COLUMN-MAJOR fill: pad each
  group to equal length and every group lands in its own column under its own
  bold heading. Groups are consistent across panels -- "Which bound binds",
  "Libbrecht's chamber", "Our conditions", "Reference".
- Because the heading now carries the provenance, the entries shrank to the
  number alone: "sigma = 1e-1" instead of "Libbrecht chamber, sigma = 1e-1",
  and "sigma = 2.7e-3  wall BC" instead of "Molaro wall BC, sigma = 2.7e-3".
- Panel 5's "<- shaded" suffix moved into the title, which now says outright
  that the shading follows sigma = 2.7e-3.
- Trimmed dead space: VIEW_BETA bottom 2e2 -> 1e3 (nothing runs below 1.5e3),
  VIEW_NX bottom 3e2 -> 4e2 (lowest curve 6e2).

---

## 2026-09-02 (end of day) — Legends outside; IguanaTeX snippets actually work

- **The IguanaTeX snippets did not generate.** I had validated every block
  wrapped in \[ ... \], assuming that is what IguanaTeX does. It is not:
  IguanaTeX supplies \documentclass ... \begin{document} and drops the paste
  into the body in TEXT mode. `aligned` needs an enclosing math mode, so every
  block failed. Fixed: `aligned` -> `align*`, and every bare-math block now
  carries its own \[ ... \]. Also dropped \gtrsim (block 13) -- it is amssymb
  and IguanaTeX's preamble carries only amsmath. All 13 blocks now compile
  against the EXACT wrapper Jackson pasted, and the test in docs/tex/README.md
  uses that wrapper rather than my earlier assumption.
- Legends moved OUT from under the curves to beneath the axes, 2-3 columns,
  every panel. The reserved-band trick from this morning was worse than the
  problem: it cost 1.5-2 decades of axis on panels 2/3/5/6 and still crowded
  the region labels. Windows reverted to the tight zoom (alpha_c to 2.5,
  Nx to [3e2, 1e6]).
- Legend labels shortened: "Libbrecht chamber, sigma = 1e-1" -> "sigma = 1e-1
  (Libbrecht)", the bound formulas dropped to "B-KINETIC binds (K&P 45)" (the
  formulas live in docs/tex now), and the reference lines to "production run",
  "Nx = 1e4 practical ceiling" and so on.
- Layout switched to constrained_layout, which accounts for legends outside the
  axes; tight_layout does not. The overview is now 3 x FIGSIZE wide (30 x 15)
  so each panel gets the standalone's area -- at 24 in the outside legends
  overflowed into the neighbouring column.

---

## 2026-09-02 (later still) — Slide-scale type; legends off the curves

- Jackson: legends in panels 5 and 6 were covering curves, panel 4's equation
  needed more air, and the whole thing has to read from the back of a room.
- Type scale is now one scheme for both the standalone panels and the
  overview: title 15, labels 13, ticks 11.5, legend 10.5, notes 10. Standalone
  figures went 8.2x5.8 -> 10x7 in, the overview 20x11.4 -> 24x13.6, so the
  per-panel area stays comparable and one font scheme serves both.
- Every crowded panel now reserves a band that NO DATA CAN REACH and parks its
  legend there, instead of floating it over the curves:
    * panels 2, 3: above alpha_c = 1 (physical ceiling, drawn and labelled)
    * panel 5: above the largest resolvable eps
    * panel 6: below the smallest Nx
  Panels 5 and 6 stay exact reciprocals because Nx's reserved band is at the
  bottom and eps's is at the top -- same band, seen from the two ends.
- Panel 4's equation gets a blank line between it and the heading, and sits in
  the title at 13 pt.
- Panel 1: moved the factor-27 note to the upper left (the only corner with no
  curve, legend or arrow) -- it was colliding with the kink callout in the
  overview -- and trimmed the dead headroom to ylim 0.45.

---

## 2026-09-02 (later) — Zoomed the Libbrecht panels to the usable window

- Jackson: the beta_sub equation box was sitting on panel 4's legend, and the
  full-range axes wasted the detail on values we could never run. Fixed both.
- Panel 4's defining equation moved into the axis title, out of the data area
  entirely -- it cannot collide with a legend or a curve there.
- Every axis now clipped to the LITERATURE band opened up two decades each
  side: alpha_c [1e-5, 1], beta_sub [2e2, 2e8] around M&F's [2e4, 2e6], Nx
  [3e2, 1e6] around the Nx = 1e4 practical ceiling. The eps window is the
  exact reciprocal of the Nx window so panels 5 and 6 read against each other.
  The usable structure is now legible: at sigma = 1e-1 Libbrecht's own alpha_c
  sits right on the 1e4 ceiling, and the crossovers are visible instead of
  being compressed into a hairline.
- `_mark_offscale()` puts a triangle on the frame where a curve exits and
  states the extreme it reaches (3e+32, 1e-30, ...), so clipping never reads
  as missing data. Markers stagger by curve index and carry a white backing.
- Dropped the molecular-limit line (a = 3.19 A, Lx/a = 1.4e6) from panels 5-6:
  it belongs to the regime we just zoomed out of. It stays in the README, the
  .tex and the console summary.
- Added `docs/tex/constraints_iguanatex.txt` -- Jackson wants to paste single
  equations into PowerPoint via IguanaTeX, not \input fragments. 13 numbered
  snippets, amsmath/amssymb only, each regression-tested to compile inside
  \[ ... \] (which is what IguanaTeX wraps them in). Header lines are
  %-commented so a sloppy copy is harmless. The recheck script is in
  docs/tex/README.md.

---

## 2026-09-02 — Regenerated the Libbrecht intractability case; constraints in LaTeX

- Added `preprocess/plot_libbrecht_constraints.py`: six slide-sized panels plus
  a combined overview, tracing Libbrecht (2017) sigma0(T) -> alpha_c -> beta_sub
  -> eps -> Nx. Every eps and mesh number goes through `comp_eps.compute_eps()`
  itself, so the figures cannot drift from the sizer. Output in
  `studies/libbrecht_kinetics/` with a README stating the argument.
- The load-bearing input is v_n. Eq. (45) only bites if the front moves, so the
  default is Molaro's own MEASURED integrated neck rate 3.416e-9 m/s from
  `studies/molaro_2019/alpha_c_estimate.csv`, not comp_eps's generic 1e-9. The
  README says explicitly that `--vn_feature` makes Eq. (45) collapse to
  eps <= R_feat and is therefore circular, not a rebuttal.
- Four sigmas plotted: Libbrecht's own chamber conditions (1e-1, 1e-2) against
  ours (2.7e-3 = the Molaro wall undersaturation; 4.5e-4 = d0/rho_fillet). At
  -20 C the wall-BC value already demands Nx = 1.4e7 (9.2e13 nodes, eps = 0.47
  angstrom) and the fillet value 3.4e31 (5.7e62 nodes, eps 25 orders below a
  water molecule). Production is Nx = 5394 / 14.5 M nodes at eps = 0.118 um.
- B-VAPOR and B-CURV never bind for this geometry: D*_ia < D_v makes B-HEAT the
  tighter of the two K&P Eq. 43 channels, and 0.05*R_ave = 4.3 um is a decade
  loose. Only B-HEAT and B-KINETIC appear in the shading.
- Added `docs/tex/`: `constraints_kinetics.tex`, `constraints_eps.tex`,
  `constraints_derived.tex`, `constraints_libbrecht.tex` (beamer-ready
  fragments, no preamble) and a compilable `constraints_standalone.tex`.
  Transcribed from comp_eps.py, including the scaled/unscaled beta trap, the
  D_heat channel choice, and the alpha_sub/M_sub identity.
- .gitignore: narrow exception so `studies/libbrecht_kinetics/*.png` is tracked,
  matching the existing carve-outs for verification and docs figures.

---

## 2026-08-26 — Vector contour export for Figure 1

- Added `postprocess/contour_svg.py`: extracts a level set from a `solV_*.vts`
  snapshot and writes it as a clean, editable SVG (one `<path>` per loop, um
  coordinates, mm page size, pt stroke weight, named groups for Inkscape).
  Motivation: ParaView/matplotlib SVG exports are clip-path soup and cannot be
  edited by hand.
- `--mirror` reflects about r = 0 and rejoins the open axisymmetric contour into
  a single closed outline; the half-plane form closes along the axis instead.
- Exported the phi_i = 0.5 outline at the final step (288, t = 7210.28 s) of the
  Dv10 curvature-BC Molaro T = -20 C r = 14 um run, mirrored and half-plane, to
  the manuscript Figure1__SinteringMechanisms folder in iCloud, each with a PNG
  preview. Final neck half-width ~30.5 um (from 14 um).

---

# Activity log — enceladus_DSM

Newest entries first.

## 2026-08-18 (later) — Rung 0 reviewed; four output/consistency bugs fixed

Rung 0 (10-step IC check, job1059799) ran clean. The IC is confirmed
quantitatively, not just visually: ice volume +0.0073 % vs analytic (matching
the diffuse-tail correction (4pi^3/3)R eps^2/V exactly, and the same figure
verify_grain_shrinkage.py predicts), interface area +2.40 %, 5-95 % band 5.95
eps / 8.42 elements vs 5.889 / 8.33, 1-99 % band 9.31 eps / 13.17 vs 9.190 /
13.00, and phi on the axis at the contact = 0.5748 -- a pure cusp sampled
7.8e-8 m off the exact point, under half a grid spacing. No SNES/KSP failures.

Scope caveat: the run covered 5.03 ms of 7200 s (7e-5 %), so tot_ice is
unchanged to 7 s.f. and the CFL limiter was never exercised (0 caps, 0
rollbacks, dt reached only 1.8e-3 s). That check moves to rung 1.

Not problems, explained and left alone:
- "OUTPUT DATA ERROR" is a bare PetscPrintf, gated on n_out_log <= 0, and only
  fired because rung 0 passed -t_out_log 0. Production sets -t_out_log 80.
- Temperature outrunning vapour is xi_v = 1e-3: over 5.03 ms heat penetrates
  268 um and vapour 10.5 um. Expected.
- 3 unused PETSc options (-ksp_gmres_restart under bcgs, -Lz/-Nz in 2D).

Fixed:
- **The banner printed kinetics the solver does not use.** Under
  -alpha_pointwise 1 it showed tau_sub = 6.8e1 s from the IGNORED -beta_sub0
  (alpha_c = 8.0e-3) against a true 7.8e0 s -- an 8.7x error I read when sizing
  -dtmax. Now evaluates AlphaCondensation + SubKinetics at the IC state and
  prints those, flagging when alpha_c sits ON the -alpha_hi clamp (which
  -alpha_c0 1.0e-1 against -alpha_hi 1.0e-1 does exactly).
- D_v banner line printed the raw STP constant, 13 % above what VaporDiffus
  evaluates at -20 C. Now prints D_v(T0) with D_v0 in parentheses.
- enceladus_main.c built the startup tau_sub from the same uncorrected
  dif_vap while SubKinetics uses VaporDiffus -- a 13 % disagreement between
  two paths that should match. Both now use D_v(T0).
- **comp_eps.py's xi_v check had the inequality backwards**, warning when xi_v
  EXCEEDED rho_vs/rho_i. xi_v scales diffusion and the source but not storage,
  so the storage error is rho_vs/(xi_v rho_i) and larger xi_v is MORE
  quasi-steady. Now warns below 10*rho_vs/rho_i and reports the storage error
  (0.092 % at the default) and tau_vap instead of a pass/fail.
- estimate_vapor_bc_molaro.py printed the vapour transient as L^2/D_v, omitting
  xi_v and understating tau_vap by 1000x. Now prints tau_vap with the scaling
  and its margin against eps/v_n, which scales as L^2 (dom2 37x, dom3 17x,
  dom4 9.3x -- all fine, dom4 worth watching).
- Deleted a material_properties.c comment referencing a mob_scale field that
  does not exist.

Domain size: Jackson raised it as unsettled. Rung 1 now runs dom2 AND dom3 at
h = 1.0000 to measure the effect rather than rely on the (a/L)^3 estimate.
Nothing submitted yet.

---

## 2026-08-18 — Interface-width conventions, then the Molaro 2019 campaign setup

**Part 1: eps is Moure & Fu's, not Kaempfer & Plapp's W.**

- Re-derived the K&P (2009) Eqs. 43/45 bounds against this solver's own
  Allen-Cahn equation, from the paper text rather than from the existing
  docstrings. The solver's well is (1/2)phi^2(1-phi)^2 with an eps^2 gradient
  term -- M&F's -- so the equilibrium profile is logistic and `eps` IS M&F's
  eps. K&P's Eq.(33) uses a +-1 well with profile tanh(x/(sqrt(2) W)), so their
  W = sqrt(2)*eps. Confirmed structurally against the M&F PDF (their Eq. 5
  triple well is the same 1/2-normalised form, and they write the constraint in
  eps, not W).
- Nothing numerical was wrong: a1 = 5 and a2 = 0.1581 are the constants for
  M&F's well and already absorb the sqrt(2). Verified exactly --
  5 == (5 sqrt(2)/8)*4, and a1*a2 = 0.79 in eps-units vs 0.78 in W-units.
  Every committed opts case reproduces its eps and mesh unchanged.
- Renamed the mislabelled "beta_eff/beta_target > 0.9" diagnostic to the
  thin-interface expansion parameter delta = a1 a2 eps/(D beta_HK). K&P Eq.(40)
  quantifies the error you make when you DON'T compensate, but tau_sub carries
  the a2 terms and does compensate, so the residual is O(delta^2). Now also
  reports delta_ice - delta(D*_ia), the honest uncertainty on beta from
  compensating at the mean while violating the ice channel.
- B-CURV was R_ave with the safety factor on top; now eps <= eps_over_R * R_ave
  (--eps_over_R, default 0.05) with no safety factor. chi means curvature only.
  Fixed the "~7.5 elements" claim (it is 8.33) and the "2 sqrt(2) eps Karma
  width" (the asymptotic W is sqrt(2) eps, so h = eps/sqrt(2) is W/2, i.e. 2
  elements per width -- mid-range for Karma-Rappel practice, not 4 per a width
  that does not exist). Docs and comments cut ~40 %.

**Part 2: Molaro et al. (2019) -20 C campaign setup, all parameters derived.**

- `preprocess/estimate_alpha_c_molaro.py` inverts Gibbs-Thomson + series
  resistance at the neck against their measured dr/dt. The attachment-only
  limit needs no geometric assumption and gives alpha_c >= 0.139; the alpha_c=1
  transport ceiling shows the observed rate exceeds what vapour alone delivers.
  The saturation sweep shows only 1.25x headroom above alpha_c = 0.1, because
  Lambda = beta_HK D_v/rho is already 0.29 there. **alpha_c cannot be tuned to
  close the gap** -- it tops out near 50 % of their rate at any value, which is
  the surface-diffusion share docs/molaro_validation_synthesis.md sec. 4 finds
  independently. Raising alpha_c also tightens the K&P bound (L* ~ 1/alpha_c).
- `preprocess/estimate_vapor_bc_molaro.py` does the same for the Dirichlet
  humidity: h = 0.99728 at L/a_eff = 3, nearly alpha_c-independent but strongly
  L-dependent, so the previously fitted h = 0.998 does not transfer from the
  old 121 um domain.
- Two findings about their data. (1) Fit the LARGE grain: its 78-min trend is
  4.4 sigma and its least-squares slope (-2.93 %) reproduces the -3 % their
  caption quotes; the caption's -4 % for the small grain is not in the same
  table (endpoint -0.69 %, min -2.40 %, fit -0.89 %, all inside its own +-2 %).
  (2) The differential is NOT a ripening signal -- the imposed wall
  undersaturation outweighs d0*(2/R_sm - 2/R_lg) by 343x.
- eps = 2.5852e-07 set by the NECK FILLET (rho >= 6 eps at their FIRST data
  point), which resolves their entire r/R = 0.19-0.37 range. Every K&P ceiling
  is looser at alpha_c = 0.1, so safety never binds (0.999).
- Three domain arms at L/a_eff = 2.01/3.03/4.01 (the old geometry was 1.08, with
  the wall 20 um from the ice and a predicted ~58 % neck error), five humidity
  arms, and the batch files for both.
- `postprocess/grain_shrinkage.py` + an analytic tangent-pair gate under
  `studies/molaro_2019/verification/`. The gate caught two real defects: the
  default split plane was the midpoint between grain CENTRES, 14 um from the
  contact for unequal radii, moving +1.3 % onto the small grain (half the
  signal being fitted) -- now the radical plane; and splitting the trapezoid
  rule at an array index dropped the segment straddling the split.
- `postprocess/run_batch_measure.sh` measures on the cluster and emits one
  summary.csv, since a dom3 arm's snapshots cannot be moved.
- Nothing has been submitted. Next step is the 10-step IC check on HPC.

---

## 2026-08-14 — plot_growth_rate.py: interface velocity v_n as a ParaView field

- New macro `scripts/paraview_macros/plot_growth_rate.py` renders the interface
  normal velocity derived in docs/curvature_driven_growth.md:
      v_n = -3*M*eps*kappa + (eps*alph_sub/(5*rho_ice))*(rho_v - rho_vs)
  v_n > 0 = ice growing. Adds V_normal, V_curvature, V_phasechange and
  V_normal_um_per_day, masked to the diffuse band.
- **Chains onto plot_rhovsI rather than recomputing kappa.** It consumes that
  macro's Curvature / Supersaturation / RhoVS_I arrays, so the curvature
  implementation (which took three rounds to get right) exists once.
- Emits the decomposition, not just the total, because on the wedge run the two
  terms are ~1e-11 each and **cancel by a median of 88.5%**, leaving v_n ~5e-13.
  Plotting v_n alone would hide that the net is a small residual of two large
  opposing rates.
- mob_sub / alph_sub / rho_ice / eps are read from the run's own outp.txt
  banner (they are DERIVED at startup from d0_sub0/beta_sub0 and appear in no
  .opts), with PV_* env and constant overrides; provenance is printed per value
  and the macro refuses to guess if it cannot find them.
- Cross-checked against postprocess/gt_balance.py, which reaches the same
  quantity by the independent route (sigma - d0*kappa)/beta. On the phi=0.5
  contour: concave +5.264e-13 vs +5.20e-13, convex -4.955e-13 vs -4.97e-13.
- Band medians differ from contour values (+6.28e-13 vs +5.26e-13) because
  kappa varies across the band; the docstring says to read the contour, and
  prints the instruction after running.
- Hit the ParaView builtins-shadowing trap a third time, and here it is a real
  user-path bug rather than a test artifact: clicking plot_rhovsI first runs a
  Programmable Filter that star-imports over max/min/sum in __main__, which
  macros share. Guarded with builtins bindings.

---

## 2026-08-13 — Where curvature actually drives growth: docs/curvature_driven_growth.md

Jackson asked why a growing interface reads UNDERSATURATED when the only
explicit phase-change term is proportional to (rhov - rhoIvs), and where
curvature enters at all given d0_GT was deleted. Written up in full.

- **Answer.** Capillarity is carried by the Allen-Cahn term, not by rho_vs.
  Substituting the equilibrium profile kills the AC bracket and leaves
  -3*M*eps*|grad phi|*kappa: motion by mean curvature. Projecting the residual
  onto the translation zero mode f' (Fredholm solvability, since L f' = 0 by
  translation invariance) gives
      v_n = -3 M eps kappa + (eps*alph/(5 rho_ice))(rhov - rhovs)
  i.e. the Gibbs-Thomson condition  s = d0*kappa + beta*v_n  with
      d0   = 15 M rho_ice/(alph_sub rho_vs)
      beta = 5 rho_ice/(eps alph_sub rho_vs)
- **d0 checks out to five decimals**: 1.0166e-9 m from the run's own M and
  alph_sub, vs d0_sub0 = 1.0166e-9 printed by the run, vs the physical
  gamma*Vm/(R*T) = 1.0168e-9. The hard-coded a1 = 5.0 IS the
  int f'^2 / int loc*f' ratio, (1/6eps)/(1/30).
- **beta is 22% high at leading order**, and the gap is fully accounted for:
  it equals 5*a2*eps*(1/diff_sub + 1/dif_vap)*rho_ice/rho_vs, the thin-interface
  correction the code subtracts when building tau. Implied diff_sub = 7.7785e-6
  matches enceladus_main.c:717's 0.5*(k_a/rho_a cp_a + k_i/rho_i cp_i) to four
  digits. a2 enters beta only, never d0.
- **The paradox resolves**: at the concave surface d0*kappa = -5.05e-6 and the
  measured s = -4.71e-6, so relative to its own (depressed) equilibrium the
  interface is SUPERsaturated by +3.4e-7 and grows. It is undersaturated only
  against a flat surface.
- **New `postprocess/gt_balance.py`** tests s = d0*kappa + beta*v_n on any run:
  reads M/alph_sub/rho_ice/eps from the run's outp.txt, samples the phi=0.5
  contour, reuses plot_rhovsI.py's own filter for kappa. On the wedge run,
  5702 points over 8 timesteps: median |s - d0*kappa|/|d0*kappa| falls from
  1.000 at t=0 (vapour initialized uniform, far off equilibrium) to 0.103 at
  90 d. It RELAXES onto the GT line, which is the evidence the condition is
  emergent rather than imposed. Concave grows (+5.2e-13 m/s), convex
  sublimates (-5.0e-13 m/s). Figure + CSV in docs/figures/.
- **Practical consequence recorded**: do NOT add an explicit d0_GT*kappa term
  to RhoVS_I — it would double-count capillarity the AC term already supplies.
  The lever is -d0_sub0. Also noted that Supersaturation_GT ~ beta*v_n, i.e.
  it is an interface-velocity field (positive = growing).

---

## 2026-08-12 (later still) — "Bunny ears": a cancellation, not a denominator

Jackson reported a symmetric cusp in the curvature profile across the band.
Real, and a third distinct defect after the eps_reg residue.

- **Not the denominator.** Measured the equipartition ratio binned by phi on
  the wedge: 1.073 at the band edges, 0.990 mid-band — an 8% swing, against a
  36% peak-to-trough in kappa. The denominator could not explain it.
- **Cancellation in the numerator.** `-(L - n.H.n)/G` forms a difference whose
  operands are `|1 - 2 phi|/(eps*kappa)` times larger than the result: 238x at
  phi = 0.03 for the wedge, exactly 0x at phi = 0.5. Discretization error in
  the second derivatives is amplified by exactly that profile — none mid-band,
  most at both edges. That IS the bunny-ear shape.
- **Fix: normalize first, differentiate second.** kappa = -div(n) with n built
  pointwise and handed back to VTK's gradient as a point array
  (`dsa.VTKArray` + DataSet/Association). The large f'' terms then cancel
  exactly at the divide instead of approximately inside a difference of
  discretized derivatives, and what gets differenced is an order-1 unit vector.
  No denominator at all.
- **Resolution sweep added to the gate**, which is what makes the case:
  MAX band error on an analytic disc, divergence vs bracket —
  dx/eps 0.33: 0.87% vs 4.53%; 0.70: 3.48% vs 20.44%; 1.00: 6.45% vs 38.70%.
  Real runs sit near 0.7. The old gate at dx/eps = 0.33 was too fine to see it.
- Main gate worst case 4.73% -> 0.88%. Wedge profile spread across the band
  36% -> 8%, which is about the true level-set variation plus numerics.
- Equipartition is no longer a denominator (there is none) but survives as a
  printed diagnostic of how far the profile is from equilibrium; it remains the
  denominator for METHOD = "bracket", kept only so the gate can score both.

Separately: confirmed no absolute value anywhere in the kappa path. Curvature
does NOT change sign across phi = 0.5 within one interface, and should not —
the level sets of a curved surface are nested and all curve the same way
(analytically, an R = 200 um disc gives +4926 at phi = 0.03 rising monotonically
to +5076 at phi = 0.97). The field that flips sign at phi = 0.5 is
laplacian(phi) / d2phi/dn2, the inflection of the tanh profile. The old
eps_reg residue was proportional to phi'', which is why the earlier, broken
version appeared to flip.

---

## 2026-08-12 (later still) — Curvature was dominated by a regularization artifact

Jackson flagged the Curvature field on the wedge run
(`batch_2026-08-07__07.26.38_wedge_bc/2D_wedge_band_90deg__...`) as wrong: the
sign did not mirror between the ice block's left and right surfaces. He was
right, and the cause was the regularization, not the normal direction.

- **Diagnosis.** The two-term form `-L/G + (g.H.g)/G^3` with Tikhonov
  `G^2 = |grad phi|^2 + eps_reg^2` does not cancel between its terms. On a FLAT
  interface, where kappa must be 0, it leaves `-eps_reg^2 * phi'' / G^3` —
  zero at phi = 0.5, sign-flipped either side, divergent at the band edges.
  Measured on the wedge: -1.0e5 /m at phi = 0.03 against a real curvature of
  ~4.5e3. The residue is a function of phi, so it reads as "negative in air,
  positive in ice" at BOTH surfaces — which is exactly the failure to mirror
  that was reported.
- **Fix.** Factor G out: `kappa = -(L - n.H.n)/G`. A flat interface has
  L = n.H.n identically, so the bracket vanishes for any denominator. The unit
  normal `n = grad phi/|grad phi|` already points into the ice by construction;
  no sign convention needed.
- **Denominator from equipartition** (Jackson's suggestion), `|grad phi| =
  phi(1-phi)/eps`: analytic in phi, so it does not amplify noise where the
  measured gradient decays at the band edges. `PV_DENOM=gradient` switches to
  the measured gradient. The macro now prints the measured equipartition ratio
  (1.041 median on the wedge, so the assumption holds there).
- **eps is now a real input,** not a regularization knob — it sets the
  denominator scale, so kappa scales as 1/eps if it is wrong. Resolution order
  PV_EPS > EPS constant > `-eps` in the staged .opts > estimate off the data,
  and the source is printed every run.
- **The gate missed this** because it only sampled |phi-0.5| < 0.05, where the
  residue is zero by construction. Rewritten to score the WHOLE band against
  each level set's own radius, plus a flat-interface case that pins kappa = 0
  at every phi. Before: max error 283%. After: 4.73% (flat case exactly 0.00%).
- Diverging colormap now centred on zero; auto-rescaling to an asymmetric range
  had been putting white off-zero, compounding the misreading.

Wedge results after the fix: left surface uniformly -4583..-6246 /m (concave),
right surface uniformly +3036..+4243 /m (convex), no sign flip within either
band — the annulus geometry, as expected.

---

## 2026-08-12 (later still) — ParaView macros: phase isovolumes, GT curvature fields

- Added `scripts/paraview_macros/split_phases.py`: ice IsoVolume
  (`IcePhase` in [0.5, 1.01]) + air IsoVolume ([-0.01, 0.5]) on the active
  source. The isovolume half of `setup_movie_view.py`, standalone.
- Extended `scripts/paraview_macros/plot_rhovsI.py` with the Gibbs-Thomson
  correction removed from the solver on 2026-07-21, as a *post-processing*
  diagnostic: new `Curvature`, `RhoVS_I_GT` and `Supersaturation_GT` point
  arrays, with kappa = -div(grad phi/|grad phi|) rebuilt from the deleted
  `Curvature()` and evaluated via two passes of VTK's gradient operator.
  d0 = gamma*V_m/(R*T) pointwise. Nothing in the solver reads these.
- eps (for the regularization) and `-axisym` are auto-detected from the .opts
  the run script stages next to the data; the axisym mode is always printed,
  since the .vts itself cannot say, and getting it wrong halves every kappa.
- Added `scripts/paraview_macros/verify_curvature.py`, an analytic gate that
  runs the macro's own filter on synthetic logistic spheres: planar +/-1/R and
  axisymmetric +/-2/R. All four within 1.01% (near-axis 1.24%). The concave
  cases pin the sign — kappa must be POSITIVE on a convex grain or ripening
  runs backwards.
- Two real bugs the gate caught:
  - **Boundary ring.** kappa needs a second derivative, so VTK's one-sided
    stencil corrupts the outer TWO point layers, not one: 39% and 26% low on an
    analytic sphere, in planar runs as much as axisymmetric. Now repaired by
    replicating the first clean layer outward.
  - **Axis limit.** n_y/y is 0/0 on the axis; flooring the denominator was 13%
    off. Replaced with the exact limit H_yy/G (g_y vanishes there by symmetry).
- Also noted, for anyone writing a Programmable Filter: it execs in `__main__`,
  so its preamble star-imports `numpy_interface.algorithms` over the builtins
  (`max(a, b)` becomes `algs.max(a, axis=b)`) and rebinds `vtk`. Both files
  document it; the gate dodges it with underscore-prefixed globals.
- On the `phi0.325_0.5mm_T-20_s05` epsconv run the flat supersaturation is
  ~2e-15 (uniform T, rho_v initialized at rho_vs — no driving force at all)
  while `Supersaturation_GT` spans -4.4e-4 to +2.3e-2. That is the LSW
  ripening driver the 2026-07-21 removal took out, visible again.

---

## 2026-08-12 (later still) — ParaView macro to split ice/air isovolumes

- Added `scripts/paraview_macros/split_phases.py`: builds an ice IsoVolume
  (`IcePhase` in [0.5, 1.01]) and an air IsoVolume ([-0.01, 0.5]) on the active
  source, colours them solid, and hides the input beneath them.
- This is the isovolume half of `setup_movie_view.py` on its own, for when the
  two phases are wanted as pipeline objects rather than as a movie scene.
- Auto-detects the phase array name, and leaves the air volume invisible on 3D
  inputs where it would enclose the ice.
- Verified with pvpython on the `phi0.325_0.5mm_T-20_s05` epsconv run: 348k ice
  cells + 175k air cells vs 503k input cells — the overlap is the cells
  IsoVolume clips at phi = 0.5, so the two tile the domain.

---

## 2026-08-12 (later) — `--axisym` correction; the arms agree after all

- The fine arm's `neck_width.csv` had been extracted **without `--axisym`**.
  These are axisymmetric r–z runs, so the measured chord is the neck *radius*
  and the flag doubles it; every fine-arm width was half its true value.
- Corrected, the fine arm spans 9.3 → 51.1 µm (not 4.7 → 25.6). This withdraws
  the entry below: the arms differ by **1.9× at 32.8 µm, 1.2× by 51 µm**, not
  36–150×, and the fine arm *does* reach the experiment's starting width
  (15.05 h), so nothing needs extrapolating.
- The fine arm is above its own 31.6 µm floor over the compared range — the
  resolved arm the pair was built to produce. a = **0.283 ± 0.001**, inside the
  Demmenie band; the coarse arm (below its 60.1 µm floor throughout) gives
  0.229. Molaro's data gives 0.204 ± 0.053, so model and experiment do *not*
  agree on the exponent once the resolved arm is used.
- Rate offset vs Molaro: 152× (coarse), 156× (fine) — agreeing to 3 % across a
  3.6× change in ε, which is what marks it as the α_c choice, not a mesh effect.
- Rewrote `compare_mesh_pair.py` down to a single axes as asked: both arms, the
  Molaro points, a dotted power-law fit per arm, each simulation clock shifted
  so t = 0 lands where its own neck reaches 32.81 µm. No mechanism guide lines.
  Dropped the now-subsumed `meshpair_both_*` and `meshpair_aligned_*` sets.

---

## 2026-08-12 — Fine arm arrives: the mesh pair does not converge (WITHDRAWN — see above)

- Extracted the fine arm's neck curve (`neck_width.py`, 57 snapshots) and
  compared both arms against the Molaro −20 °C data.
- **The coarse mesh is not adequate.** At matched neck width the arms differ by
  36–150× in elapsed time, and re-zeroing at a shared neck does not collapse
  them — so it is not the coarse arm merely being ahead on one trajectory.
- The exponent is not converged either: equal over each arm's *own* range
  (0.19–0.21 vs 0.21–0.24), but 0.157–0.238 (coarse) vs 0.261–0.285 (fine) over
  the range they share. Refining ε moves the exponent toward 1/3; the fine arm
  lands inside the Demmenie band, the coarse arm does not.
- This supersedes 2026-08-11's headline: the coarse arm's agreement with the
  Molaro exponent, and the ~152× rate factor, are mesh properties. README
  updated where it previously said the opposite.
- Neither arm resolves its own neck (coarse 53.5 µm vs a 60.1 µm floor; fine
  25.6 µm vs 31.6 µm). The pair was expected to give one resolved arm, gave none.
- The requested alignment on the experiment's first neck (32.81 µm) works for
  the coarse arm (−7.74 h) and is impossible for the fine one — it tops out
  7.25 µm below that width, extrapolating to ~325 h (4× its run).
- New `analysis/compare_mesh_pair.py` (3-panel comparison + rate-ratio CSV).
- Fixed `--anchor-neck` in `fit_neck_growth.py`: it interpolated the anchor
  after the `t > 0` filter had dropped the t = 0 sample, which for a pre-necked
  experiment is the sample carrying the starting width — so it silently dropped
  the very dataset it was written to align. `Series` now keeps raw arrays.

---

## 2026-08-11 — mesh_pair coarse arm: exponent fits vs the Molaro −20 °C data

- Fitted the tangent-contact coarse run (79 h, α_c = 1e-3, ε = 0.87 µm) from
  `batch_2026-08-11__17.24.16_mesh_pair` against `molaro2019_fig11_T-20.csv`.
  a = 0.19–0.21 on all three fit forms, vs a = 0.18–0.25 for the data — inside
  the data's CI, at m ≈ 5, not Demmenie's m = 3.
- The tangent IC is the reason: it gives the run a physical clock zero, so the
  three forms agree instead of spanning 2x as they did on the pre-necked runs
  in `results/molaro_prenecked/` (a = 0.09–0.12 there).
- Recorded the caveats with the numbers: the whole curve sits below this arm's
  resolution floor √(12ε/R) = 0.346 (default protocol → no fittable window), so
  the fine arm settles it; the absolute rate is ~155x slow at α_c = 1e-3; the
  first ~15 samples are diffuse-interface relaxation off tangent contact.
- New `analysis/run_mesh_pair_fits.sh` (picks up the fine arm automatically) and
  `results/mesh_pair/` with three fit windows, both figures, and a README.
- Added `analysis/plot_neck_linear.py`: neck **width** vs time on linear axes,
  `w = A·t + C`, Molaro overlaid. A = 0.290 ± 0.025 µm/h over the shared width
  range (R² = 0.96); a line does *not* fit the full 79 h (R² = 0.77).
- Best result of the day: stretching the Molaro clock by one factor S = 152
  (elapsed time at equal neck size) puts the data on the model curve within its
  error bars — model vs experiment is a pure rate offset over 32.8–53.5 µm.
- The requested `w = A(t − t0) + C` is over-parameterized (t0 and C trade off
  exactly). Resolved it by anchoring: both clocks re-zeroed at a common neck
  width (32.81 µm; the model needs 7.74 h to get there from tangent contact),
  which fixes the origin from the data and leaves the one-parameter pinned form
  `w = A·t' + w_anchor`. A = 0.349 ± 0.027 µm/h model vs 52.91 ± 12.57 Molaro.
- That cross-validates the rate factor: pinned ratio 151.5× vs model-free
  elapsed-time ratio 151.8×, two independent routes to the same number. The
  free-slope 175× was an artifact of the data's truncated fit window.

---

## 2026-08-10 — Sintering growth exponent vs Demmenie et al. (2025)

- New target: Demmenie, Woutersen & Bonn, *J. Phys. Chem. Lett.* 16(8)
  2104–2109 (2025), doi:10.1021/acs.jpclett.5c00050 — two ~1 mm ice spheres at
  −3 °C in a box at ice-saturation, 2.5 h, giving r ~ t^alpha with alpha =
  0.29/0.33/0.30/0.26 (±0.01), i.e. Kuczynski m = 3, evaporation–condensation.
  That is this model's own transport route, and alpha is a *shape* observable,
  so it is independent of the ~50 % rate deficit against Molaro.
- Pulled `neck_width.py`, `plot_neck_convergence.py`,
  `calibrate_neck_geometry.py` and the Molaro validation CSVs out of the
  scratch quarantine.
- Added `postprocess/fit_neck_growth.py`. Reports three fit forms side by side
  (`d_free` = Demmenie's `C(t+t0)^a` with t0 free, `d_fixed`, Kuczynski
  `r^m − r0^m = Kt`) plus the local log-slope, because the exponent turns out
  to be mostly a property of the protocol: the same model curve gives a = 0.087
  over the solver's dense output, 0.078 at the data's sample times, 0.110 with
  t0 released.
- **Found `neck_width.py` could never read a snapshot.** Its local VTS parser
  called a `decode` helper it never imported; a bare `except Exception` then
  reported every NameError as "empty/corrupt file, likely an incomplete
  transfer" and the script exited 0 with a header-only CSV. Broke when pplib
  absorbed the decoder (a09f873) while the script sat in scratch/. Now
  delegates to `pplib.read_vts`, narrows the except, refuses to write an empty
  CSV.
- Added `-t_out_log N` (log-spaced snapshot cadence). Neither existing cadence
  spreads samples across decades, and one snapshot of the 9792×2671 sintering
  mesh is ~630 MB, so the budget is ~40 files. t0 defaults to `delt_t`, not
  `dtmin`: entries below the first accepted step are silently subtracted from
  the budget (`-t_out_log 8` gave 7 files before the fix).
- Stage-1 reanalysis (`studies/sinter_exponent/`, no compute): **the Molaro
  geometry cannot resolve an exponent at all.** The neck fillet `rho ~ r²/(2R)`
  is only resolved above `r/R >= sqrt(12·eps/R)` = 0.22 there, and the whole
  dataset spans 0.19–0.38 — a quarter decade sitting on the floor, over which
  1/3 and 1/7 are indistinguishable. The eps 0.60 and 0.86 µm arms have no
  fittable window whatsoever. Fix is bigger grains, not finer mesh: the floor
  falls only as sqrt(eps).
- Also: the Molaro *data* sits at a flat local slope of 0.185, between 1/7 and
  1/5, so it does not show Demmenie's exponent either — consistent with their
  claim that the 1/3-to-1/7 literature scatter tracks humidity control. And the
  model/data gap is far smaller on the r0-subtracting Kuczynski form (m = 5.32
  vs 4.04) than on the naive one (11.5 vs 5.4).
- Added five opts files: strict/coarse Demmenie arms (1 mm spheres, exactly
  tangent, `-ic_grain_union 1` — the additive IC would open with a spurious
  eps-dependent bridge), their shared −3 °C experiment file, and the Molaro
  pair restarted from tangent contact with an 8 h experiment. IC verified at
  the full pilot mesh: TOT_ICE = 1.047e-9 = the analytic two-sphere volume.
- Not yet run. Pilot batch is `studies/sinter_exponent/pilot_batch.txt`
  (~$5); production strict arm ~$60–110, gated on the pilot showing the neck
  climbs past r/R0 ≈ 0.15.

---

## 2026-07-31 — Merge effective_thermal_cond in; k_eff becomes an in-line diagnostic

- Merged `projects/effective_thermal_cond` into this solver. k_eff is now
  computed during the time loop (`-keff 1`) on a cadence independent of field
  output (`-keff_freq` / `-keff_t_interv`), because a snapshot is hundreds of MB
  while a k_eff sample is ~200 bytes — the old coupling capped k_eff(t) at the
  snapshot count. `-keff_replay <dir>` covers finished runs and replaces the
  standalone binary, which is now marked superseded.
- New: `include/keff.h`, `src/keff{,_field,_cell,_solve,_sample}.c`. Guards
  hard-error on non-periodic meshes, `-geom_file`, `-axisym` and `dim < 2`.
- Renamed `NASA_types.h`/`NASA_main.h` → `enceladus_types.h`/`enceladus_main.h`;
  deleted dead `env_helper.c/h`; raised `TARGET_DOFS_PER_CORE` 50k → 80k.
- Deleted the write-only `alph[]`/`mob[]` Gauss-point arrays and standardised on
  one indexing convention, `point->index + point->count*point->parent->index`,
  documented at the top of `src/keff_field.c`. The read side is a PetIGA form
  callback, so a sequential counter was never an option there.
- Fixed `-ic_type ice_slab`, which built a circular blob of radius `Ly` rather
  than a slab and self-overlapped on a periodic cell.
- **The "iterative solvers are unreliable for this problem" claim was a bug**,
  not numerics: `IGACreateMat` overrides `MATOP_CREATE_VECS`/`MATOP_DUPLICATE`
  with IGA-aware versions, AMG's coarse operators inherit them without the
  composed IGA, and `PCGAMG` died in setup on every mesh. CG+GAMG now agrees
  with LU to all ten printed digits and is 17× faster at 322k unknowns.
- Verification suite in `studies/snow_thermal/verification/` — analytic gate,
  solver comparison, IC resolution study, figures, raw logs, README.
- Two measurement notes that matter for the study: on the production mesh rule
  the IC's volume-fraction error is `≈118·eps`, so the discretisation floor on
  the safety 0.5 vs 1.0 comparison is **0.024%**, not the ~0.2% first estimated;
  and the diffuse band (9.2·eps = 9.2 µm at safety 0.5, 18.4 µm at safety 1.0)
  exceeds the median pore throat (4.0–12.6 µm) in most shakedown packings, so
  throat closure is the mechanism to watch, biasing k_eff high at low porosity.

---

## 2026-07-27 — Project created from lunar_regolith_DSM

- Created `enceladus_DSM` as a verbatim copy of `lunar_regolith_DSM` (itself
  the renamed `sublimation_pf`), so both projects start from the same
  extensively-tested two-phase solver and diverge only as the two papers
  require.
- The previous `dry_snow_metamorphism` solver was **discarded**, not merged:
  it had a wrong latent-heat density (mixture `rho` instead of `rho_ice`),
  stale `xi_v`/`xi_T` (1e-5/1e-4 vs 1e-3/1.0), Libbrecht kinetics on by
  default, `eps` sized at 273.15 K on a 2.8x-too-coarse mesh, and `p=1/C=0`
  which is incompatible with `effective_thermal_cond`. Its final state is
  preserved at tag `archive/dry_snow_metamorphism-legacy`.
- The three-phase (ice/sediment/air) formulation is gone from both projects.
- `studies/snow_thermal/` moved here from the lunar project; it is this
  project's effective-thermal-conductivity study.

---
