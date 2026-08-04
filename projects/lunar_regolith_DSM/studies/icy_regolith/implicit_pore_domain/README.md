# Implicit pore-domain approach (Effort 1)

**Question.** How does ice adhere to and metamorphose within lunar-regolith pore
space, as a function of pore geometry? Here the regolith is modelled
*implicitly* — it is the deformed top/bottom boundary of the domain, not a
simulated field. Ice / vapour / temperature evolve in the pore space between the
regolith walls (3 DOF, no solver changes). Interior regolith grains are not
representable in a single tensor-product patch; that is the explicit approach
(`../explicit_sediment_phase/`, Effort 2).

## Geometry

A two-sided pore channel whose walls are lined with regolith grains drawn from a
grain-size distribution (median radius ~50 µm, `docs/material_parameters.md`
§2.3), with a central **throat pinch** (the defining pore constriction). Built by:

    preprocess/build_geometry_regolith_pore.py

Strategies that only place ice into the pore share one mesh
`inputs/geometry/regolith_pore.dat` (the walls) and differ only in their
`.opts` (the ice); a strategy that *reshapes* the walls (`wall_divots`) gets
its own mesh. Deterministic; asserts every geometric constraint (throat gap,
wall-slope budget, no buried/overlapping ice, bump count vs `MAX_SED_GRAINS`,
`eps/R` resolution) before writing. See the header of that script for all knobs
(`R_MED`, `N_BUMPS`, throat multipliers, wall seeds, divot sizing).

## Ice placement is the study variable

How ice sits in the pore is exactly what we're probing, so the initial ice
configuration is a swappable `--ice-placement` strategy, each testing a
different adhesion hypothesis:

| Strategy (`.opts`) | grains | Hypothesis it probes |
|---|---|---|
| `flank_caps` → `2D_regolith_pore.opts` | 6 | Ice caps on grain tops; throat bare. Baseline. |
| `throat_bridge` → `2D_regolith_pore_throat.opts` | 7 | Ice suspended in the throat — does it persist/adhere at the constriction? |
| `pore_lining` → `2D_regolith_pore_lining.opts` | 15 | Small ice grains seated in wall troughs (reentrant pore corners) — distributed adhesion. |
| `wall_divots` → `2D_regolith_pore_divots.opts` | 12 | Grains adhering in **divots carved into the walls**, each seat sized to its grain. The setup for the temperature-gradient aggregation run — see below. |

Preview each config (walls + ice) at `preprocess/regolith_pore{,_throat,_lining,_divots}.png`.
Add a strategy by registering a function in `ICE_STRATEGIES`.

### `wall_divots` — grains seated in the wall

The other three strategies place ice *on* the wall. `wall_divots` cuts each
grain a **seat**: a smooth depression in the regolith surface, half-width
1.75·R_ice and depth 0.5·R_ice, with the grain centred at the bottom of it.
The wall clips the disc on both sides, so the exposed ice is a little less than
a semicircle nestled in the depression — ice adhering in a hollow of the
regolith rather than perched on a bare surface.

Mechanically a divot is just another C-infinity bump with **negative height**
appended to `-sed_grain_*` / `-top_grain_*`. `SedimentBumpField()`
(`src/initial_conditions.c`) sums them with no sign restriction, so the wall
stays C-infinity and the solver needs no special handling.

**Mesh quality is what limits how deep a divot can be.** The mesh interpolates
linearly in `v` between floor and ceiling, so element skew near a wall tracks
that wall's slope; the validated stiff-wall geometry peaks at 46°. A divot's
own slope *adds* to whatever the wall already has, so three things keep the
total in budget (`SLOPE_BUDGET = 1.10`, ~48°):

1. **Seats go at wall extrema only** — trough bottoms and grain crests, where
   the local wall slope is zero. Seating in troughs alone (the `pore_lining`
   rule) measured 55–59°, because a trough bottom is not flat enough and the
   divot's flank lands on the neighbouring grain's rise.
2. **Softer regolith relief for this strategy only** (`h/R` 0.22–0.36 vs
   0.35–0.60), with the central amplification raised to 2.8 so the throat
   survives the softening. `R_MED` is unchanged — the grains sit more embedded
   in the substrate, they are not smaller.
3. **The budget is checked per seat over its own support**, not domain-wide.
   A global check lets the first divot that reaches the limit veto every later
   seat anywhere on the wall, however flat that spot is; supports never
   overlap, so per-seat checking is equivalent, and the finished wall is still
   asserted domain-wide.

Because the walls differ, this strategy has its **own mesh**
(`regolith_pore_divots.dat`) rather than sharing `regolith_pore.dat` with the
other three. Strategies that reshape walls are listed in `CARVES_WALLS`.

Measured on the generated mesh: Jacobian strictly positive (min/max 0.465, vs
0.403 for the validated `regolith_pore.dat`), max element shear 47° vs 46°,
throat 48% Ly, and the mesh boundary agrees with the C-side bump formula to
8.7e-9 m — 1% of `eps`.

### The temperature-gradient aggregation run

`wall_divots` exists to answer: with ice distributed along both pore walls as
comparable wall-adhering deposits, **where does it collect** under a
warm→cold gradient? Monotone cold-ward migration, a pile-up at the throat
constriction, or do the seats pin it in place?

Deliberately *not* the older `2D_pore_channel_icecap*` framing (a slab of ice
at one end, one big attractor grain at the other), which fixes the end state
in the initial condition.

```bash
./scripts/HPC/submit_lunar.sh 2D_regolith_pore_divots Tgrad_T-20_G50_3d divot_seed
```

`Tgrad_T-20_G50_3d` is −20 °C, saturated, `dT/dx = −50 K/m` (left warm), with
the left/right faces pinned and top/bottom insulating so the gradient survives
at the walls. It needs no edits — and its −20 °C matches this geometry's
`-eps_valid_temp -20` lock.

Driving ΔT across `Lx = 8.14e-4` m is **4.07e-2 K**, but only ~3.5e-3 K between
neighbouring grains ~70 µm apart — weaker per-neighbour forcing than the
two-grain sanity test (1.25e-2 K over 2.5e-4 m). Little movement in 3 days is
therefore expected rather than a fault; the response is a longer `t_final`
(`Tgrad_T-20_G50_90d.opts`), not a geometry change.

## The wedge: does taper alone move ice?

A separate, deliberately minimal geometry —
`preprocess/build_geometry_wedge.py` → `2D_wedge_band.opts` — that removes
every variable the geometries above deliberately include. Two **perfectly flat**
walls tapering from a 50 µm throat to a 200 µm mouth (half-angle 14°), **no
temperature gradient**, and a **single** ice band bridging the channel. If the
ice moves, the pore geometry moved it.

This is the analytic clean-room version of the `throat_bridge` hypothesis: the
rough-walled geometries ask whether ice persists at a constriction, and this one
isolates the mechanism that would make it do so.

### Prediction: it migrates toward the NARROW end

The walls are two rays from a virtual apex off the left edge of the domain, and
a circular arc centred on that apex is perpendicular to every ray from it — so
at the model's natural 90° (Neumann) contact angle, **both menisci are
apex-centred arcs whose radius of curvature is their distance from the apex.**
The inner meniscus has the ice on the *outside* of its circle, so it is
**concave** (κ = −1/r₁); the outer is **convex** (κ = +1/r₂). The concave face
has the lower equilibrium vapour density, so vapour flows inward.

The counterintuitive part is that the narrow end has the *larger* |κ|. That is
true, but κ carries a sign, and the sign is what sets which way material goes.
The energy accounting agrees and needs no flux argument at all. At 90° Young's
equation gives γ_wall-ice = γ_wall-vap, so sliding the band costs no wall energy
and only the two menisci count:

```
interfacial length  L = 2α(r₁ + r₂)
ice area            A = α(r₂² − r₁²) = α(r₂−r₁)(r₂+r₁)
```

With `S = r₁+r₂` the constraint gives `r₂−r₁ = A/(αS)` and `L = 2αS`, so at
fixed area **L falls monotonically as S shrinks** — toward the apex. Curvature
rises toward the apex while interfacial *area* falls, and area is what carries
the energy. Same sign as two familiar cases: sintering necks grow rather than
pinch off, and capillary condensation fills the smallest pores first.

**The run is the test of this, not a demonstration of it.** A drift toward the
wide end would be a real finding, or a sign error worth hunting.

### Why a band and not a circle

A circle meets a wall at 90° only if its centre lies *on* that wall, which
cannot hold for both walls at once. Seeding a circle would open the run with a
contact-angle relaxation transient far larger than the slow taper-driven
migration being measured. The annulus is already at the right contact angle —
verified at exactly 90.000000° at all four wall contacts.

### Two things the solver needed

- **Affine wall baselines** (`-wall_{bot,top}_{y0,slope}`). The solver rebuilds
  the wall from the bump lists only to map (u,v) → (x,y) for the IC. Bumps have
  compact support and **cannot express a linear ramp**, so without this the ice
  would be seeded against a flat channel while the mesh is a wedge. Defaults
  (0, 0, Ly, 0) are a no-op; all four geometries above regenerate byte-identical.
- **The band IC** (`-wedge_apex_*`, `-wedge_band_r{1,2}`), added as an additive
  feature in `FormInitialMultiGrains2D` alongside `ice_shell_*`/`ice_flat_*`
  rather than as a new `ic_type`.

A straight wall is represented **exactly** by the mesh (B-splines reproduce
linear functions at their Greville abscissae): the mesh boundary matches the
solver's wall formula to 4e-20 m, against ~9e-9 m for the bump walls.

### Running it

```bash
# config check first — confirm the 90 deg contact before spending an allocation
./scripts/Studio/run_lunar.sh  2D_wedge_band Tgrad_T-20_G0_10s config_check
# production
./scripts/HPC/submit_lunar.sh  2D_wedge_band Tgrad_T-20_G0_90d wedge_migration
# the answer
python3 postprocess/track_wedge_band.py --dir <run folder> \
    --save-csv wedge_band.csv --save-fig wedge_band.png --no-show
```

**Pair only with a zero-gradient experiment** — a thermal gradient would swamp
the effect entirely (the G50 runs drive ~1.2e-3 supersaturation; this drives
1.7e-6).

**Expect a modest result.** The driving force is `d₀(1/r₁ − 1/r₂) = 1.7e-6`,
about **0.13×** the 50/150 µm ripening pair that showed clear change in ~30
days. If it is too slow the knob is a **steeper wedge**, not a longer run — the
force scales as `4α²/w`:

| half-angle | throat | d₀Δκ | vs. ripening |
|---|---|---|---|
| 3.8° | 125 µm | 1.4e-7 | 0.01× |
| **14°** | **50 µm** | **1.7e-6** | **0.13×** |
| 17° | 50 µm | 7.2e-6 | 0.53× |
| 25° | 50 µm | 1.6e-5 | 1.1× |

The mesh is ~0.17M control points (`# DOF_GRID: 497 332`), about a third of the
divot geometry.

## Running

    # from the project root
    ./scripts/Studio/run_lunar.sh 2D_regolith_pore        30day_T-20_h1.00_arrh
    ./scripts/Studio/run_lunar.sh 2D_regolith_pore_throat 30day_T-20_h1.00_arrh
    ./scripts/Studio/run_lunar.sh 2D_regolith_pore_lining 30day_T-20_h1.00_arrh

(3 DOF; `run_lunar.sh` sizes the rank count from the `# DOF_GRID` comment.)

## Regenerating / varying the geometry

    python3 preprocess/build_geometry_regolith_pore.py --ice-placement flank_caps
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement throat_bridge --tag throat
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement pore_lining  --tag lining
    python3 preprocess/build_geometry_regolith_pore.py --ice-placement wall_divots  --tag divots

To vary the *pore* geometry (wall roughness, throat tightness, GSD), edit the
module constants and re-run; the `.dat` and all `.opts` regenerate.

## Parameters — eps is proper, but temperature-locked

`eps = 8.5840e-7` is **not** an arbitrary/loose value: it is exactly what
`preprocess/comp_eps.py` recommends for **T = −20 °C, α_c = 1.341e-2**, and it
reproduces the validated `2D_ripening_two_sided` reference run bit-for-bit
(same `beta_sub0=5.9216e5`, `d0_sub0=1.0166e-9`). The binding K&P constraint is
the temperature-dependent **kinetic bound**, so eps is grain-size independent
here — every strategy passes `eps/R_ave < 5%` (3.0%, 3.0%, 4.8%, 4.1%).

**Consequence:** these geometries are valid **only at −20 °C**, paired with a
−20 °C experiment (`30day_T-20_h1.00_arrh`, which sets the matching kinetics).
For any other temperature — notably lunar PSR temperatures (40–120 K,
`docs/material_parameters.md`) — **recompute eps** and regenerate:

    python3 preprocess/comp_eps.py --Lx 8.14e-4 --Ly 2.6933e-4 \
        --Rave <smallest_ice_R> --T0 <degC> --alpha <alpha_c>
    # then set T0_C / ALPHA_C / EPS in build_geometry_regolith_pore.py and re-run

## Status & first gate

- First gate: a short −20 °C run that passes mass conservation
  (`postprocess/plot_mass.py`), confirming the geometry + IC behave.
- Then: sweep ice placement × throat tightness, and quantify the resulting
  adhesion morphology (which configs coarsen toward the walls vs. the throat,
  SSA evolution, ice–wall contact length).
