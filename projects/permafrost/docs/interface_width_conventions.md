# Interface Width: Conventions, Measurement, and Mesh Sizing

*Written 2026-07-09, after the d0phys two-grain runs. Companion to the
`eps` section in CLAUDE.md and to `preprocess/comp_eps.py`.*

This note settles three recurring questions:

1. What **is** the interface width of this model, and why do different ways of
   measuring it give different numbers?
2. How should we **measure** it (e.g. in ParaView) and what value counts as
   "on spec"?
3. How should **eps**, the **mesh size h**, and the **initial condition**
   be chosen so everything is mutually consistent?

---

## 1. The equilibrium profile has no edge — every "width" is a chosen cutoff

The 1D equilibrium solution of this model's double-well free energy is the
**logistic profile**

```
phi(s) = 1 / (1 + exp(-s/eps))
```

where `s` is the signed distance from the interface midpoint and `eps` is the
**decay length** — the single true length scale of the interface. The profile
approaches 0 and 1 *asymptotically*: there is no point where the interface
"ends." Consequently, any number quoted as "the interface width" is really
"the width between two chosen cutoff contours," and all such widths are fixed
multiples of `eps`:

| convention            | definition                              | width      | elements at h = eps/√2 |
|-----------------------|-----------------------------------------|------------|------------------------|
| Karma width           | 2√2·eps (tanh-argument convention, K&P) | 2.83·eps   | **4.0**                |
| slope width           | 1 / max\|phi'\| (tangent ramp at mid)   | 4.00·eps   | 5.7                    |
| 10%–90% band          | 2·ln(9)·eps                             | 4.39·eps   | 6.2                    |
| 5%–95% band           | 2·ln(19)·eps                            | 5.89·eps   | 8.3                    |
| **1%–99% band**       | 2·ln(99)·eps                            | 9.19·eps   | **13.0**               |
| 0.1%–99.9% band       | 2·ln(999)·eps                           | 13.8·eps   | 19.5                   |

**Every row is the same interface.** A ParaView measurement between the
phi = 0.01 and phi = 0.99 contours is a perfectly valid measurement — it is
the 1%–99% row — but it must be compared against the 1%–99% prediction
(9.19·eps ≈ 13 elements on the standard mesh), never against a rule-of-thumb
number that was calibrated in a narrower convention. Comparing "13 measured"
against "7–10 expected" mixes rulers and produces false alarm.

Verified empirically (2026-07-09 d0phys run, eps = 4.8534e-7, h = 3.41e-7):
the relaxed interface measures 13.0 ± 0.1 elements across 1%–99% for 1400+
consecutive snapshots — exactly the analytic 2·ln(99)·eps = 13.1.

## 2. Why the literature's "4–10 elements" rule uses the narrow conventions

The rule of thumb exists to control **discretization error**, which lives
where the solution has large derivatives — the steep **core** of the profile.
The tails (between phi = 0.01 and 0.05, say) are nearly flat and cost the
basis nothing to represent. So published resolution guidance (Karma & Rappel,
Provatas & Elder, ...) is stated in core-width units: the Karma/tanh width or
the slope width. Translated to the 1%–99% ruler, "4–10 per Karma width" reads
"13–32 per 1%–99% band."

Our mesh rule `h = eps/√2` (built into `comp_eps.py`) delivers **exactly 4
elements per Karma width** — the standard floor — equivalently ~6 across
10%–90% and ~13 across 1%–99%.

## 3. How to measure interface width (and resolution adequacy)

- **Keep using the ParaView contour method** (contour phi = 0.01 and 0.99,
  count elements between them, e.g. with Surface With Edges). Just compare
  against **9.19·eps/h ≈ 13** as the on-spec value for the standard mesh.
- Equivalently, contour 0.1/0.9 and compare against ~6.2, which maps directly
  onto the literature's 4–10 range.
- A count *below* the expectation means the profile is sharper than the model
  equilibrium (e.g. an unrelaxed IC — see §5) or the mesh is coarser than
  h = eps/√2. A count *above* it (persistently, after relaxation) would mean
  something is genuinely fattening the profile and is worth investigating.
- **Never** hand-tune eps or h to make the visible band hit a remembered
  element count; the band width is an *output* (see §4).

## 4. The one-way logic chain: physics → eps → profile → h → IC

The correct/consistent procedure for this double-well system flows strictly
one way:

1. **Physics bounds eps from above** (K&P 2009 Eqs. 42–46, implemented in
   `comp_eps.py`): eps must be well below the smallest grain radius, the
   kinetic length d0/(beta_sub·v_n), and the heat/vapor diffusion lengths.
   The script picks eps = safety · min(bounds). Note the bounds use the
   PHYSICAL capillary length d0 = gamma·V_m/(R·T) (Eq. 13).
2. **eps fixes the equilibrium profile analytically** — logistic with decay
   length eps, all widths in the §1 table. This is a property of the energy
   functional; it is not independently adjustable. The kinetic coefficients
   (lambda, tau_sub, mob_sub, alph_sub via M&F SI Eq. 9) are *derived from
   this same eps* — which is why the assembly must use `user->eps` verbatim.
3. **h resolves the profile core**: h = eps/√2 gives the standard 4 elements
   per Karma width. If a mesh-convergence study shows more accuracy is
   needed, reduce h (e.g. eps/2 or eps/(2√2)); the 1%–99% count then rises
   (18, 26, ...) as an automatic consequence. Refining h changes nothing
   about the width in meters — it only re-tiles the same interface.
4. **The IC reproduces the equilibrium profile**:
   `phi = 1/(1 + exp(-(R - d)/eps))` = `0.5 - 0.5·tanh(0.5·(d - R)/eps)`.
   Any other steepness just adds a one-time relaxation transient (the model
   restores its own width within ~tens of steps) plus a spurious kick to the
   coupled vapor/thermal fields.

## 5. Historical pitfalls (do not repeat)

- **The `eps_model = 0.75·eps` residual/Jacobian scaling** (removed): scaling
  eps inside the assembly while the kinetic coefficients are derived from the
  unscaled `user->eps` breaks the Kaempfer–Plapp asymptotic mapping — the
  effective capillary length and kinetic coefficient silently shift ~25%,
  corrupting the physical d0 calibration. If a narrower interface is ever
  genuinely required, change eps itself (opts) and re-run `comp_eps.py` so
  mesh and kinetics move together.
- **Stale IC steepness** (fixed 2026-07-09): the grain ICs long carried
  `tc = 1/(√2·0.75·eps)` — 1.89× steeper than equilibrium, a leftover of the
  removed scaling. Every run began at ~7 elements (1%–99%) and relaxed to 13
  over ~60 steps. ICs now use the equilibrium `tc = 0.5/eps`.
- **Judging eps from the visible band** (see CLAUDE.md): the on-screen
  diffuse band is ~6–9× eps depending on cutoff, so "the interface looks 2×
  too wide, halve eps" reasoning is off by that ruler factor and also
  invalidates the K&P bound derivation. Always recompute eps with
  `comp_eps.py`.

## 6. Open item

Whether 4 elements/Karma width is sufficient *for this problem's accuracy
targets* is an empirical question: run the same GT-physics case at h = eps/√2
and h = eps/2 (eps unchanged!) and compare grain-mass/interface evolution
curves. If they overlap, the standard mesh is converged; if not, the delta is
the coarse-mesh error bar.
