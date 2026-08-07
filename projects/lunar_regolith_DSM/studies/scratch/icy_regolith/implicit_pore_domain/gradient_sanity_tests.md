---
marp: true
theme: default
paginate: true
title: Temperature-gradient sanity tests (icy-regolith, flat-domain pre-study)
---

<!--
This file is BOTH the study record and the slide source.
Render to slides:
  marp --pdf  gradient_sanity_tests.md        # or --pptx / --html
  pandoc -t pptx -o deck.pptx gradient_sanity_tests.md
One slide per "##"; "---" separates slides. Speaker notes are in HTML comments.
-->

# Temperature-gradient sanity tests
## Icy-regolith metamorphism — flat-domain pre-study (Effort 1)

Jackson Baglino · Sublimation phase-field model

<!--
Goal of the deck: show the advisor/postdoc the plan for a clean set of
gradient-driven metamorphism tests that precede the implicit regolith-pore runs.
-->

---

## Why these tests

- The icy-regolith study needs a **clean temperature-gradient baseline** before
  we add pore-wall geometry.
- Our prior two-grain runs placed grains **too close** (surface-to-surface gap
  ≈ 0.5 R), which **distorted the vapour transport pathway** between them.
- Fix: **shrink grains to 50 µm** and **space them well apart** (and off the
  vertical walls), so what we measure is physics, not proximity artifacts.

> Deliverable: 7 flat-domain simulations, one fixed gradient, **no solver or
> boundary-condition changes** — just new geometry inputs.

<!-- Figure to add: old close-packed pair vs. new well-separated pair. -->

---

## The physics we are probing

Two drivers of ice redistribution, made to **compete** or **reinforce**:

- **Temperature-gradient metamorphism** — warm ice has higher equilibrium
  vapour pressure → **warm grain sublimates, cold grain grows**.
- **Ostwald ripening (curvature)** — **large grains grow, small grains shrink**.

Fixed setup for the whole series:

- Gradient **dT/dx = −50 K/m** (LEFT warm → RIGHT cold). Mass moves left→right.
- −20 °C, saturated. 50 K/m is chosen to **maximize the rate** so we see the
  same morphological change in fewer timesteps.

<!-- Schematic: horizontal T ramp, left-warm arrow (sublimation) → right-cold
     arrow (deposition); a second arrow set for big↔small ripening. -->

---

## Method

- **Flat rectangular domain**, semicircular ice grains resting on the floor
  (analytic tanh initial condition — no meshed regolith yet).
- **eps + mesh** from `comp_eps.py` (Kaempfer & Plapp) at −20 °C, α_c = 1.341e-2:
  `eps = 8.584e-7`, `h = eps/√2` (~7.5 elements across the interface).
- Gradient imposed by pinning the **left/right walls** (Dirichlet T); top/bottom
  left insulating — no new physics.
- **Time-synchronized output** (100 snapshots at fixed physical times) → every
  run is directly comparable snapshot-for-snapshot.

Runs on HPC as one batch → `batch_<timestamp>_Tgrad_sanity/`.

---

## Tests 0–1 — baselines

**Test 0 · single grain** (Lx = 0.375 mm)
- One isolated 50 µm grain at the domain centre.
- *Expect:* grain **translates toward the cold (right)** side, total ice
  conserved (~1e-4%). Confirms the gradient + IC behave.

**Test 1 · well-separated equal pair** (Lx = 0.75 mm)
- Two 50 µm grains at Lx/4 & 3Lx/4 (5.5 R apart).
- *Expect:* warm (left) grain sublimates → redeposits on the cold (right) grain.
  Confirms the **separation fix** removed the vapour-pathway artifact.

---

## Tests 2 & 6 — competing vs. complementing  ★ headline

Big grain (150 µm) vs. small (50 µm), 1.0 mm domain. **Only the side differs.**

| | Big grain on… | Gradient wants it to… | Ripening wants it to… | Net |
|---|---|---|---|---|
| **Test 2** | WARM (left) | shrink | grow | **compete** |
| **Test 6** | COLD (right) | grow | grow | **reinforce** |

- **Test 2 (competing):** which driver wins, and does the big grain survive?
- **Test 6 (complementing):** both drivers agree → big grain grows fastest.

<!-- Put the two side by side; this is the scientifically interesting contrast. -->

---

## Tests 3–5 — three-grain series

Three grains, equally spaced (0.2 / 0.5 / 0.8 Lx), 1.125 mm domain.

- **Test 3 · equal** — three 50 µm grains; the gradient sets a warm→cold
  feeding chain. Baseline.
- **Test 4 · small centre** (25 µm) — ripening drives the middle grain to
  shrink; gradient adds direction. Drivers on the **most disadvantaged** grain.
- **Test 5 · large centre** (100 µm) — ripening makes the middle grain the
  **hub** that grows; gradient biases which flanker feeds it.

---

## What we expect (hypotheses to confirm/refute)

| # | Prediction |
|---|---|
| 0 | Single grain drifts cold-ward; mass conserved. |
| 1 | Left shrinks, right grows; clean transfer (no proximity artifact). |
| 2 | Competing — big warm grain shrinks **slower** than a small one would; watch for survival. |
| 6 | Complementing — big cold grain grows **fastest** of the series. |
| 3 | Monotone warm→cold mass shift across the chain. |
| 4 | Small centre disappears first (ripening + gradient both against it). |
| 5 | Large centre grows; net feed comes preferentially from the warm flanker. |

> Framed as **falsifiable** predictions so the discussion has hooks.

---

## Next steps

- Calibrate `t_final` from tests 0/1 (well-separated grains transfer slower than
  the old close-packed pair — 3 days may not reach the end-state).
- Carry the **validated gradient** into the **implicit regolith-pore geometry**
  (deformed walls = regolith), then add ice–regolith / vapour–regolith wall BCs.

---

## Reference — the 7 geometries

| # | File | Lx [mm] | Ly [mm] | Grains (R, µm) | Role |
|---|---|---|---|---|---|
| 0 | `tgrad_2D_L375um_eps0.86um_1grain`       | 0.375 | 0.20 | 50                | single baseline |
| 1 | `tgrad_2D_L750um_eps0.86um_2grain_equal` | 0.75  | 0.20 | 50, 50            | equal pair |
| 2 | `tgrad_2D_L1mm_eps0.86um_2grain_bigL`  | 1.00  | 0.30 | **150**, 50       | competing |
| 6 | `tgrad_2D_L1mm_eps0.86um_2grain_bigR`  | 1.00  | 0.30 | 50, **150**       | complementing |
| 3 | `tgrad_2D_L1.12mm_eps0.86um_3grain_equal` | 1.125 | 0.20 | 50, 50, 50        | 3-grain baseline |
| 4 | `tgrad_2D_L1.12mm_eps0.86um_3grain_smallC`| 1.125 | 0.20 | 50, **25**, 50    | small centre |
| 5 | `tgrad_2D_L1.12mm_eps0.86um_3grain_bigC`  | 1.125 | 0.25 | 50, **100**, 50   | large centre |

All pair with experiment `tgrad_T-20_h1.00_3d_G50` (dT/dx = −50 K/m, −20 °C, 3 days).

**Run (HPC, one batch):**

```bash
./scripts/HPC/submit_batch.sh --tag Tgrad_sanity --tests "\
tgrad_2D_L375um_eps0.86um_1grain:tgrad_T-20_h1.00_3d_G50,tgrad_2D_L750um_eps0.86um_2grain_equal:tgrad_T-20_h1.00_3d_G50,\
tgrad_2D_L1mm_eps0.86um_2grain_bigL:tgrad_T-20_h1.00_3d_G50,tgrad_2D_L1mm_eps0.86um_2grain_bigR:tgrad_T-20_h1.00_3d_G50,\
tgrad_2D_L1.12mm_eps0.86um_3grain_equal:tgrad_T-20_h1.00_3d_G50,tgrad_2D_L1.12mm_eps0.86um_3grain_smallC:tgrad_T-20_h1.00_3d_G50,\
tgrad_2D_L1.12mm_eps0.86um_3grain_bigC:tgrad_T-20_h1.00_3d_G50"
```
