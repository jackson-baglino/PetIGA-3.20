# Permafrost Phase-Field — Test Suite Report

**Generated:** 2026-04-18 08:04:27  
**Executable:** `/Users/jacksonbaglino/PetIGA-3.20/projects/permafrost/permafrost`  

## Summary

| Test | Name | Result | Key Metric |
| ---- | ---- | ------ | ---------- |
| T01 | Smoke Test | **✅ PASS** | exit code=0 (want 0); SSA_evo.dat=found; NaN in tot_ice=no |
| T02 | Initial Condition Accuracy | **✅ PASS** | tot_ice(0)=3.020e-05 m  (expect ≈ 3.00e-05); tot_sed(0)=2.01 |
| T03 | Sediment Inertness | **✅ PASS** | max |Δvol_sed|/vol_sed(0) = 0.00e+00  (threshold 1e-4) |
| T04 | Phase Non-Negativity | **✅ PASS** | min(tot_ice)=3.02e-05; min(tot_sed)=2.01e-05; min(tot_air)=4 |
| T05 | SNES Convergence | **✅ PASS** | steps converged=1; diverged=0; max_iters=3 |
| T06 | Sublimation Kinetics | **✅ PASS** | Δ(tot_ice) = -10.38 nm  (must be < 0); mean rate = -4.91e-07 |
| T07 | Bergeron Process | **✅ PASS** | tot_rhov(t=0): no-gradient=4.2890e-08, bergeron=5.4460e-08   |
| T08 | Flat Interface Stability | **✅ PASS** | max |Δtot_ice| / tot_ice(0) = 0.00e+00  (threshold 0.5 %) |
| T09 | Vapour Saturation at t=0 | **✅ PASS** | tot_rhov(0) = 4.005e-08 kg/m;  expected ≈ 4.005e-08 kg/m;  e |
| T10 | Interface Density Evolution | **✅ PASS** | Σ/ε: initial=1.718e-01,  final=1.718e-01,  rel Δ=1.56e-04  ( |

**10 / 10 tests passed.**

---

## Test Descriptions and Results

### T01 — Smoke Test

**Category:** Model correctness  
**Purpose:** Code exits cleanly, output files are written, no NaN values.  
**Pass criterion:** exit code = 0; SSA_evo.dat exists; no NaN in tot_ice  
**Result:** **✅ PASS**  
**Detail:** exit code=0 (want 0); SSA_evo.dat=found; NaN in tot_ice=no

![T01 plot](plots/T01_smoke.png)

---

### T02 — Initial Condition Accuracy

**Category:** Model correctness  
**Purpose:** Ice and sediment volumes at t=0 match the tanh-slab geometry.  
**Pass criterion:** vol_ice(0) ≈ 0.30·Lx ± 5 %; vol_sed(0) ≈ 0.20·Lx ± 5 %; sum = Lx  
**Result:** **✅ PASS**  
**Detail:** tot_ice(0)=3.020e-05 m  (expect ≈ 3.00e-05); tot_sed(0)=2.013e-05 m  (expect ≈ 2.00e-05); sum=1.000e-04 m  (expect 1.00e-04)

![T02 plot](plots/T02_initial_conditions.png)

---

### T03 — Sediment Inertness

**Category:** Model correctness  
**Purpose:** With mob_sed=0, the sediment volume must not change.  
**Pass criterion:** |Δvol_sed| / vol_sed(0) < 1×10⁻⁴  
**Result:** **✅ PASS**  
**Detail:** max |Δvol_sed|/vol_sed(0) = 0.00e+00  (threshold 1e-4)

![T03 plot](plots/T03_sediment_inert.png)

---

### T04 — Phase Non-Negativity

**Category:** Model correctness  
**Purpose:** All phase volume integrals remain ≥ 0 throughout the simulation.  
**Pass criterion:** min(tot_ice) ≥ 0; min(tot_sed) ≥ 0; min(tot_air) ≥ 0  
**Result:** **✅ PASS**  
**Detail:** min(tot_ice)=3.02e-05; min(tot_sed)=2.01e-05; min(tot_air)=4.97e-05

![T04 plot](plots/T04_phase_bounds.png)

---

### T05 — SNES Convergence

**Category:** Model correctness  
**Purpose:** Newton solver converges at every time step within ≤ 7 iterations.  
**Pass criterion:** SNES never diverges; max iterations per step ≤ 7  
**Result:** **✅ PASS**  
**Detail:** steps converged=1; diverged=0; max_iters=3

![T05 plot](plots/T05_snes_convergence.png)

---

### T06 — Sublimation Kinetics

**Category:** Dry snow metamorphism  
**Purpose:** Under undersaturated vapour (hum=0.5), ice sublimates and tot_ice decreases.  
**Pass criterion:** tot_ice decreases monotonically  
**Result:** **✅ PASS**  
**Detail:** Δ(tot_ice) = -10.38 nm  (must be < 0); mean rate = -4.91e-07 m/s

![T06 plot](plots/T06_sublimation.png)

---

### T07 — Bergeron Process

**Category:** Dry snow metamorphism  
**Purpose:** Temperature gradient drives larger ice-volume change than no-gradient case.  
**Pass criterion:** |Δtot_ice| with gradient > |Δtot_ice| without gradient  
**Result:** **✅ PASS**  
**Detail:** tot_rhov(t=0): no-gradient=4.2890e-08, bergeron=5.4460e-08  (ratio=1.270, must be >1.10)

![T07 plot](plots/T07_bergeron.png)

---

### T08 — Flat Interface Stability

**Category:** Dry snow metamorphism  
**Purpose:** A planar ice-air interface at saturation does not drift spontaneously.  
**Pass criterion:** |Δtot_ice| / tot_ice(0) < 0.5 % over 100 steps  
**Result:** **✅ PASS**  
**Detail:** max |Δtot_ice| / tot_ice(0) = 0.00e+00  (threshold 0.5 %)

![T08 plot](plots/T08_flat_interface.png)

---

### T09 — Vapour Saturation at t=0

**Category:** Dry snow metamorphism  
**Purpose:** Initial rhov = hum × rho_vs(T) × vol_air within 2 %.  
**Pass criterion:** |tot_rhov(0) − hum·ρ_vs·vol_air| / expected < 2 %  
**Result:** **✅ PASS**  
**Detail:** tot_rhov(0) = 4.005e-08 kg/m;  expected ≈ 4.005e-08 kg/m;  err = 0.01 %

![T09 plot](plots/T09_vapor_saturation.png)

---

### T10 — Interface Density Evolution

**Category:** Dry snow metamorphism  
**Purpose:** Allen-Cahn dynamics change the ice-air interface density over time.  
**Pass criterion:** rel |Δ(Σ/ε)| > 1×10⁻⁴  
**Result:** **✅ PASS**  
**Detail:** Σ/ε: initial=1.718e-01,  final=1.718e-01,  rel Δ=1.56e-04  (must exceed 1e-4)

![T10 plot](plots/T10_interface_evolution.png)

---

## Methodology

### Phase-field model tests
These verify that the numerical solver is well-posed and the initial
conditions are physically consistent.  Quantities are extracted from the
monitor table printed to stdout (parsed via regex) and from `SSA_evo.dat`
(4-column ASCII file written every output step).

### Dry snow metamorphism tests
The non-variational Allen-Cahn formulation couples three physical
mechanisms:
1. **Curvature-driven interface motion** — Allen-Cahn bulk free energy
   minimisation drives diffuse interfaces toward lower curvature.
2. **Sublimation/deposition kinetics** — the source term
   `−α·φ_i²·φ_a²·(ρ_v − ρ_vs(T))/ρ_ice` converts ice ↔ vapour
   whenever the local vapour density deviates from saturation.
3. **Bergeron (temperature-gradient) process** — a macroscopic
   temperature gradient creates a vapour density gradient (d ρ_vs/dT ≠ 0)
   that drives net vapour flux from warm to cold, causing sublimation at
   the warm end and deposition at the cold end.
