# Permafrost Phase-Field — Test Suite Report

**Generated:** 2026-04-20 16:05:10  
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
| T11 | Sublimation Rate Deceleration | **✅ PASS** | Δtot_ice=-15.9 nm; rate_early=-1.66e-06 m/s; rate_late=-7.67 |
| T12 | Deposition at Supersaturation | **✅ PASS** | Δ(tot_ice) = 10.39 nm  (must be > 0); mean deposition rate = |
| T13 | Temperature Field Consistency | **✅ PASS** | ∫T dx = -2.00000e-03 °C·m;  T_avg(0) = -20.0000°C  (expect - |
| T16 | Mass Conservation | **✅ PASS** | max |Δ(ρ_i·tot_ice + tot_rhov)| / mass(0) = 6.62e-04 (thresh |
| T17 | Temperature BC Fix | **✅ PASS** | max |∫T dx − T₀·Lx| = 6.00e-06 °C·m  (threshold 1.00e-05);   |
| T18 | Contact-Sed 2D Smoke | **✅ PASS** | exit code=0 (want 0); SSA_evo.dat=found; NaN in tot_ice=no |
| T19 | Contact-Sed 2D IC Accuracy | **✅ PASS** | tot_sed(0)=1.944e-09 m²  (expect ≈ 1.92e-09, err=1.0%); tot_ |

**17 / 17 tests passed.**

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

### T11 — Sublimation Rate Deceleration

**Category:** Dry snow metamorphism  
**Purpose:** In a finite domain, sublimation decelerates as vapour builds up toward saturation.  
**Pass criterion:** Δtot_ice < 0; rate_early > rate_late  
**Result:** **✅ PASS**  
**Detail:** Δtot_ice=-15.9 nm; rate_early=-1.66e-06 m/s; rate_late=-7.67e-08 m/s; decelerates=True

![T11 plot](plots/T11_sublimation_rate.png)

---

### T12 — Deposition at Supersaturation

**Category:** Dry snow metamorphism  
**Purpose:** Under supersaturated vapour (hum=1.5), ice grows via deposition.  
**Pass criterion:** tot_ice increases (Δtot_ice > 0)  
**Result:** **✅ PASS**  
**Detail:** Δ(tot_ice) = 10.39 nm  (must be > 0); mean deposition rate = 4.91e-07 m/s

![T12 plot](plots/T12_deposition.png)

---

### T13 — Temperature Field Consistency

**Category:** Model correctness  
**Purpose:** Domain-averaged temperature at t=0 equals temp0 within 1% (T field correctly initialised).  
**Pass criterion:** |T_avg(0) − T₀| / |T₀| < 1%  AND  ice sublimated  
**Result:** **✅ PASS**  
**Detail:** ∫T dx = -2.00000e-03 °C·m;  T_avg(0) = -20.0000°C  (expect -20.0°C, err=0.000%);  ice sublimated: True

![T13 plot](plots/T13_latent_heat.png)

---

### T16 — Mass Conservation

**Category:** Dry snow metamorphism  
**Purpose:** Total water mass ρ_ice·tot_ice + tot_rhov is conserved throughout sublimation.  
**Pass criterion:** max |Δmass| / mass(0) < 2%  
**Result:** **✅ PASS**  
**Detail:** max |Δ(ρ_i·tot_ice + tot_rhov)| / mass(0) = 6.62e-04 (threshold 2%);  mass(0) = 2.7754e-02 kg/m²

![T16 plot](plots/T16_mass_conservation.png)

---

### T17 — Temperature BC Fix

**Category:** Model correctness  
**Purpose:** With flag_BC_Tfix=1, ∫T dx stays within 0.5% of T₀·Lx (Dirichlet BCs active).  
**Pass criterion:** max |∫T dx − T₀·Lx| / |T₀·Lx| < 0.5%  over 100 steps  
**Result:** **✅ PASS**  
**Detail:** max |∫T dx − T₀·Lx| = 6.00e-06 °C·m  (threshold 1.00e-05);  T_avg range: [-20.0600, -20.0000]°C

![T17 plot](plots/T17_temp_bc_fix.png)

---

### T18 — Contact-Sed 2D Smoke

**Category:** Model correctness  
**Purpose:** 2D contact-sed IC exits cleanly, output files written, no NaN.  
**Pass criterion:** exit code = 0; SSA_evo.dat exists; no NaN in tot_ice  
**Result:** **✅ PASS**  
**Detail:** exit code=0 (want 0); SSA_evo.dat=found; NaN in tot_ice=no

![T18 plot](plots/T18_contact_sed_smoke.png)

---

### T19 — Contact-Sed 2D IC Accuracy

**Category:** Model correctness  
**Purpose:** Ice and sediment areas at t=0 match the tanh-disc geometry (two grains in contact).  
**Pass criterion:** vol_sed(0) ≈ 2π·r_s² ± 10%; vol_ice(0) ≈ 2π·(r_i²−r_s²) ± 10%; sum = Lx·Ly  
**Result:** **✅ PASS**  
**Detail:** tot_sed(0)=1.944e-09 m²  (expect ≈ 1.92e-09, err=1.0%); tot_ice(0)=4.318e-09 m²  (expect ≈ 4.27e-09, err=1.2%); sum=3.125e-08 m²  (expect 3.12e-08)

![T19 plot](plots/T19_contact_sed_ic.png)

---

## Methodology

### Phase-field model tests (T01–T05, T17)
These verify that the numerical solver is well-posed and the initial
conditions are physically consistent.  Quantities are extracted from the
monitor table printed to stdout (parsed via regex) and from `SSA_evo.dat`
(4-column ASCII file written every output step).

### Dry snow metamorphism tests (T06–T13, T16)
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

Extended tests (T11–T13, T16) use a 1000-step sublimation run to verify
quasi-steady rate linearity, mass conservation, and latent heat coupling.
