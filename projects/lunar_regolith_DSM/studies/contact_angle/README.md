# Prescribed contact angle at the regolith wall

## Why

The substrate in this two-phase model is the domain boundary, not a phase
field. Until this study it had no say in the physics: `assembly.c` returned
immediately at every boundary quadrature point and no boundary form was ever
enabled, so every wall was natural Neumann `dphi/dn = 0`. The solver's own
banner said as much — *"dphi/dn = 0 for the phase field (a 90° contact angle
where ice meets the wall)"*.

A 90° wall is an accident of omission, not a modelling choice. Ice meeting
regolith has a real contact angle set by the three surface energies, and until
it can be prescribed, every wall-bounded geometry here — pore channels, wedges,
the icy-regolith work — runs the wrong substrate physics.

## What was added

A wall free-energy surface integral,

```
F_wall = ∫_Γ f_w(φ) dΓ,   f_w(φ) = γ_as − γ_ia cos θ · h(φ),   h(φ) = φ²(3−2φ)
```

so `f_w(0) = γ_as` (wall|vapour) and `f_w(1) = γ_is` (wall|ice), with Young's
equation `γ_ia cos θ = γ_as − γ_is` fixing θ. Because `h'(0) = h'(1) = 0` the
term is inert in both bulk phases and acts only where the interface meets the
wall, so it cannot shift the bulk equilibria.

Its variation gives the natural condition `dphi/dn = cos θ · φ(1−φ)/ε`, i.e.
`m·n = cos θ` with `m` the interface normal into the ice and `n` the outward
wall normal.

### The term has no ε and no γ in it

This surprised us enough to be worth recording. The interior form is the
gradient flow of a *dimensionless* functional whose interfacial excess is
exactly `ε/6` per unit area, so `γ_ia` cancels against that normalisation, and
the residual's `(3M/ε)` prefactor cancels the remaining ε:

```
R_bnd[a][0] = −3·M·cos θ·φ(1−φ)·N0[a]
```

`cos θ` is the entire physical content. `demo/Metamorph.c` writes the same term
as `−3·ε·mob·|∇φ|·cos θ`, which agrees identically because `|∇φ| = φ(1−φ)/ε` on
the equilibrium profile — but the `h(φ)` form needs no `1/|∇φ|` regularisation
(Metamorph clamps it at 1e-5) and its Jacobian is an exact mass-matrix block.
The full derivation is above `WallH()` in `src/assembly.c`.

## Usage

Surface energies are the input; the solver derives θ and prints it.

```
-wall_faces y0,y1        # which domain faces are regolith
-gamma_ia 0.109          # ice-vapour   (default: -Sigma_i)
-gamma_as 0.300          # air-regolith
-gamma_is 0.2455         # ice-regolith  ->  theta = 60 deg
```

Defaults are a strict no-op: `γ_is = γ_as` gives `cos θ = 0`, and without
`-wall_faces` no boundary form is enabled at all. `-contact_angle_deg` exists
as a debug override that bypasses Young's equation and says so in the banner.

Only constraint: `|γ_as − γ_is| ≤ γ_ia`. The three-phase triple-well constraint
`γ_is > (γ_ia/2)(1−cos θ)` does **not** apply — that came from requiring a
positive `Sigma_s`, and there is no sediment phase or triple well here.

## Verification

### Unit gates — `verification/verify_wall_bc.sh`  ✅ all passing

| Gate | What | Result |
|---|---|---|
| G1 | Without `-wall_faces`, bit-identical to the merge-base solver | `SSA_evo.dat` identical |
| G2 | Analytic wall Jacobian vs finite differences, 8 angles | ≤ 1.1e-9 (the FD floor) |
| G3 | Boundary surface measure, 3 meshes × 4 face sets × 3 angles | ≤ 7.5e-15 |
| G4 | Sign: wetting advances the contact line | correct at 5 angles |

G3 is the one that could not be settled by reading the source: on a plain
Cartesian patch PetIGA takes its `detS = 1.0` branch instead of computing a
geometric surface Jacobian, and nothing says whether the resulting face measure
is right. It is, to 2e-16.

G2 deliberately does **not** use `-snes_test_jacobian`. That only runs inside a
TS step, where TSALPHA's restart solve re-enters it repeatedly and a debug PETSc
build takes minutes on even a tiny mesh. It also tests the wrong thing: the full
system is dominated by the interior form's `3M/ε` terms and its FD error is
identical to seven digits with the wall term on and off. Differencing both `J`
and `F` between `costhet` on and off isolates the wall block, which is what the
gate actually checks — while `||Jbnd·v||` is reported so it cannot pass
vacuously.

θ = 90 is gated the other way round: `cos θ = 0` must make the block *exactly*
zero, not merely small. It is, at 2e-28.

### Measurement gate — `verification/verify_contact_angle_measure.py`  ✅ passing

Synthesises bridges whose angle is known in closed form (the model's own tanh
profile in the *exact* signed distance to a circular arc) and pushes them
through the production `measure()` rather than a reimplementation. All eleven
angles from 15° to 165° come back exact to 0.000°.

This caught a real bug before any physics run. `pplib.circle_radius` uses the
affine form `2ax + 2by + c = x²+y²`, which cannot represent a straight
interface — it is the `R → ∞` limit — and θ = 90 in a channel is exactly
straight. It returned **78 ± 25°**. `contact_angle.py` therefore fits an
implicit conic `A(x²+y²) + Bx + Cy + D = 0`, where `A = 0` *is* the line and the
normal comes from `∇F`.

The measurement also agrees with the analytic clipped-disc angle on real solver
output: the t=0 channel IC measures **131.805°** against
`arccos(−(Ly/2)/R) = 131.810°`, with `R_arc` recovered to 0.01%.

### Pilot — ✅ θ_inf = 60.36 ± 0.04° vs Young's 60.000°

One run, θ = 60°, on the coarsest mesh (ε/R = 1/25). Details and the
fit-window table in `verification/pilot/`. The relaxation is a clean single
exponential with τ = 14.4 days; a **free** three-parameter fit (asymptote
fitted, not assumed) gives θ_inf = 60.36 ± 0.04°, and the best-conditioned
window gives 59.93 ± 0.03° — **0.07° from Young**.

Strong evidence, not the validation: one angle, one mesh, and an extrapolation
rather than a directly observed equilibrium.

### Physics sweep — ⏳ not yet run

`verification/sweep_tests.txt` defines 11 runs: the five-angle sweep, ε
convergence at fixed θ, and the sessile drop. **These have not been run.**

## Running it

```bash
make
bash studies/contact_angle/verification/verify_wall_bc.sh              # G1-G4
venv_lunar/bin/python studies/contact_angle/verification/verify_contact_angle_measure.py
```

`t_final` is now sized from the pilot's measured τ = 14.4 days: 65 days, i.e.
4.5 τ, which lands within ~1° of equilibrium. `contact_angle.py` additionally
reports θ_inf from the free exponential fit, so the equilibrium angle is
recovered even if a run stops short.

The sweep (11 runs — one batch, not chained submits):

```bash
git push
./scripts/HPC/submit_batch.sh --tag contactangle \
    --tests-file studies/contact_angle/verification/sweep_tests.txt
```

`contact_angle.py` runs automatically in `run_postprocess.sh` for any run that
declares `-wall_faces`, writing `contact_angle.csv` and `plots/contact_angle.png`.

## Acceptance

|θ_measured − θ_Young| < 2° across the sweep, with the four per-snapshot
estimates (two menisci × two walls) agreeing within their spread, and the error
decreasing with ε/R.

## Known caveats

`docs/gt_deficit/` records a 1.22× gap between requested and realised β, from
`tau_sub`'s thin-interface corrections being added and never subtracted back.
It affects interface *kinetics*, not the equilibrium angle, so it should not
contaminate these results — but it is the first thing to suspect if θ comes out
systematically off while the unit gates stay clean.

`-axisym` combined with a wetting wall is untested; the `rw` weight is carried
into the boundary form for consistency but nothing exercises it.

`run_postprocess.sh` previously used whatever `python3` was on PATH, which on
this machine has no numpy. Every plotting step failed and the only symptom was
an empty `plots/` directory after a run that reported success. It now prefers
`venv_lunar` and warns loudly if it cannot find a working interpreter.
