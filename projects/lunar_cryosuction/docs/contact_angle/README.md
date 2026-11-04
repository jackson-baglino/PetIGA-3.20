# Prescribed ice–regolith contact angle

Synthesis of the contact-angle documents as of 2026-10-21. This feature is in
the current solver on `main` and is carried forward. The derivation is
`contact_angle.tex` (PDF beside it); the evidence lives in
`studies/contact_angle/`. Every number below is copied from the file cited
next to it.

## What it does

The regolith is the domain wall, not a phase field, so it acts through a
boundary condition. Cahn's wall free energy is added to the bulk functional:

```
F_wall = ∮_Γ f_w(φ) dΓ,   f_w(φ) = γ_as + (γ_is − γ_as)·h(φ),   h(φ) = φ²(3 − 2φ)
```

`h'(0) = h'(1) = 0`, so the term is inert in both bulk phases and acts only
where the ice–vapor interface meets the wall. Its first variation gives the
natural condition

```
∂_n φ = cos θ · φ(1 − φ)/ε,   cos θ = (γ_as − γ_is)/γ_ia     (m·n = cos θ)
```

with θ measured through the ice. Young's equation is an output, not an
assumption: the wall supplies `γ_as − γ_is`, and the bulk calibration supplies
`γ_ia` (`contact_angle.tex` §2.1).

Residual and Jacobian on wall faces (`src/assembly.c`, `WallH()` and the
`FACE_WALL` branch):

```
R_a = −3·M·cos θ·(h'(φ_c)/6)·N_a        J_ab = −3·M·cos θ·(h''(φ_c)/6)·N_a·N_b·dφ_c/dφ
```

No ε and no surface energy appear in the coefficient, because
`M = mob_sub = 2·γ_ia·L` already contains `γ_ia` (`contact_angle.tex` §3).
`φ_c` is φ clamped to [0, 1]. The clamp is required: unclamped, `h'` changes
sign outside [0, 1] and the wall term drives φ further out of range. This
happened once, at ε/R = 1/100: a −2e-6 excursion grew to −0.05, the residual's
domain guard returned zero, SNES reported convergence, and the run froze while
reporting success (`contact_angle.tex` §4).

## Using it

| Option | Meaning |
|---|---|
| `-wall_faces "y0,y1"` | which domain faces are regolith (none by default) |
| `-gamma_is`, `-gamma_as` | ice–regolith and air–regolith energies [J/m²] |
| `-gamma_ia` | ice–air energy seen by the wall; default `-gamma_ia_bulk` (0.109) |
| `-contact_angle_deg` | debug only: sets θ directly, bypassing Young's equation |
| `-wall_bot_y0`, `-wall_bot_slope`, `-wall_top_*` | affine walls for the wedge geometry |
| `-test_wall_measure`, `-test_wall_jacobian` | debug gates that run once and exit |

The solver forms cos θ, prints θ, and refuses the run if `|γ_as − γ_is| > γ_ia`.
That is the only admissibility constraint in this two-phase model. The
three-phase constraint `γ_is > (γ_ia/2)(1 − cos θ)` does **not** apply here.
It came from requiring Σ_s > 0, and this model has no sediment phase
(`contact_angle.tex` §5). Defaults are a strict no-op: without `-wall_faces`
no boundary form is enabled, and `γ_is = γ_as` gives θ = 90°.

## Evidence

| What | Result | Source |
|---|---|---|
| Unit gates G1–G4: bit-identical without walls; wall Jacobian vs FD; boundary measure; sign | all pass; Jacobian ≤ 1.1e-9 (FD floor); measure ≤ 7.5e-15; θ = 90° block exactly 0 | `studies/contact_angle/README.md`, `verification/verify_wall_bc.sh` |
| Angle measurement on synthetic bridges, 15°–165° | exact to 0.000° | `verification/verify_contact_angle_measure.py` |
| Channel sweep, 30°–150°, batch 2026-09-12 | max 0.079°, RMS 0.040°; error halves from ε/R 1/25 to 1/50 | `verification/sweep_2026-09-12/` |
| Channel sweep, batch 2026-09-13 (the numbers quoted in the .tex) | max 0.127°, RMS 0.062° | `verification/sweep_2026-09-13/` |
| Larger channel and unconfined droplet, batch 2026-09-14 | RMS 0.191° (channel), 0.540° (droplet); runs stopped at 2.4–3.8 τ, short of convergence | `campaign_2026-09-14/` |
| Wedge on a curvilinear mesh, sealed BCs, batch 2026-09-19 | +0.039° (60°), −0.047° (120°) | `wedge_2026-09-19/` |
| Wedge, reservoir BCs | ice +21.3 % (60°) / −21.8 % (120°); angle off Young by −1.06° / −0.15° | `wedge_2026-09-19/` |
| Interface velocity vs series-resistance theory, 44 runs, α_c = 1e-3 and 1e-2 | channel measured/theory 0.98–1.03 for most runs; θ = 30°, 150° at α_c = 1e-3 run 8–15 % high (unexplained) | `velocity_plan_2026-10-01/README.md` |

All angle results are against a 2° acceptance criterion.

Two physical lessons came out of the campaigns. Growth needs both a vapor
reservoir and a thermal bath. With an insulated box (`-flag_BC_Tfix 0`) the
latent heat of a 1e-5 % deposit warms the domain enough to cancel the driving
supersaturation (`campaign_2026-09-14/`, group B). And a moving interface holds
a dynamic angle, not Young's angle, because `σ = d0·χ + β·v_n`
(`wedge_2026-09-19/`).

## Open items and inconsistencies

- **Not fully validated.** The term recovers the static (Young) angle, but
  the implementation does not apply something like a free-flow boundary
  condition at the moving contact line (user, 2026-10-22). To be discussed
  before the dynamic results above are relied on.
- `studies/contact_angle/README.md` says `-gamma_ia` defaults to `-Sigma_i`.
  That is stale: Σ_k were removed on 2026-09-17, and the code default is
  `-gamma_ia_bulk` (`src/lunar_main.c`).
- The thin-interface β offset was measured as 8.7e5 s/m on the August wedge
  batch (`docs/enceladus_carryover.md`) and 1.3e5 s/m on the October runs. The
  cause of the 6.8× difference has not been traced.
- The excess velocity of the strongly curved menisci (θ = 30°, 150°,
  α_c = 1e-3) is not explained. The flat-front approximation is the first
  suspect.
- The droplet runs and the 2026-09-14 channel runs did not reach equilibrium.
  The README estimates 153 d (channel) and 321 d (droplet) to match the earlier
  precision.
- `-axisym` combined with a wetting wall is untested.
- The `.tex` quotes the residual as a single code line. The code now computes
  the same term through a `flux` variable shared with the open-face branch.
  The two are equivalent, but the quote is no longer verbatim.
