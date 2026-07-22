# Material Parameters: Water Ice, Water Vapor, and Lunar Regolith

> ⚠️ **PARTIALLY HISTORICAL (banner added 2026-07-21).** The ice/vapor
> parameters remain current. The **sediment / triple-line surface-energy**
> material (γ_iv/γ_is/γ_sv, contact-angle derivation, the η > 0 constraints)
> describes the removed three-phase model and is **reference for Effort 2**,
> not the current two-phase code. Note the diagnosed fix: parameterize by
> contact angle θ so only γ_is > (γ_iv/2)(1 − cos θ) is a free constraint —
> the prior γ_is = 0.033 violated it. See
> `studies/icy_regolith/explicit_sediment_phase/`.

**Date:** 2026-05-15 (banner 2026-07-21)  
**Model:** Sublimation phase-field (permafrost2.c)  
**Target environment:** Lunar permanently shadowed regions (PSRs), T ≈ 40–120 K

> **Status note (2026-05-27).** Material properties documented here are
> unchanged by the §27 vapor-penalty overhaul. For the *numerical* parameter
> defaults (`k_pen`, `α_pen`, `ξ_v`, `ξ_T`, `NRmin`, etc.) consult
> [`model_description.md`](model_description.md) §12 — that table tracks
> the current solver configuration.

---

## 1. Surface / Interfacial Energies

These enter the Kim-Steinbach three-phase Allen-Cahn model through the combination
energies η defined as:

```
ηᵢ = γ_iv + γ_is − γ_sv      (drives ice interface mobility)
ηₛ = γ_sv + γ_is − γ_iv      (drives sediment interface mobility)
ηₐ = γ_iv + γ_sv − γ_is      (drives air interface mobility)
```

and through the Young contact angle at the ice-sediment-vapor triple line:

```
cos θ = (γ_sv − γ_is) / γ_iv
```

### 1.1 Ice–Vapor Surface Energy (γ_iv)

| Value | T [K] | Method | Reference |
|-------|--------|--------|-----------|
| 0.109 J/m² | 273 K (0 °C) | contact angle / nucleation | Petrenko & Whitworth (1999) [1] |
| 0.105 J/m² | 253 K (−20 °C) | same source, slight T-dependence | [1] |
| 0.098 J/m² | 200 K | extrapolation of [1] fit | [1], [2] |
| 0.065 J/m² | 100 K | molecular dynamics | Sega et al. (2007) [3] |
| ~0.060 J/m² | 77 K | estimated / MD | [3], [4] |

**Temperature dependence** (Petrenko & Whitworth 1999):
```
γ_iv(T) ≈ 0.109 − 1.5×10⁻⁴ × (T − 273)   [J/m², T in K]
```

**Current model value:** `γ_iv = 0.109 J/m²` (0 °C reference).  
For lunar PSR temperatures (40–120 K) a reduced value ~0.060–0.070 J/m² is
more physically appropriate.

### 1.2 Ice–Regolith (Solid) Surface Energy (γ_is)

The ice-mineral interfacial energy depends on mineral composition. Lunar
regolith is dominated by plagioclase feldspar, pyroxene, ilmenite, and
agglutinates; the nearest terrestrial analogues are silicate glasses.

| Value | Material | Method | Reference |
|-------|----------|--------|-----------|
| 0.057 J/m² | ice on quartz (SiO₂) | contact angle | Israelachvili (2011) [5] |
| 0.030 J/m² | ice on basalt glass | contact angle | Jellinek et al. (1999) [6] |
| 0.020–0.060 J/m² | ice on silicates (general range) | review | Dash et al. (2006) [7] |
| 0.040 J/m² | ice on plagioclase | estimated / [5,6] | — |

**Current model value:** `γ_is = 0.020 J/m²`.  
A better estimate for lunar silicates would be **0.030–0.057 J/m²**; the
lower end of this range yields a contact angle consistent with ice being
highly non-wetting on dry lunar regolith.

### 1.3 Regolith–Vapor (Sediment–Air) Surface Energy (γ_sv)

Dry silicate surfaces have high surface energy relative to ice.

| Value | Material | Method | Reference |
|-------|----------|--------|-----------|
| 0.300–0.500 J/m² | SiO₂ / quartz | fracture / contact angle | Parks (1984) [8] |
| 0.400 J/m² | plagioclase feldspar | estimated from [8] | — |
| 0.600–1.000 J/m² | pyroxene, olivine | fracture energy | Momber (2015) [9] |
| 0.300 J/m² | lunar simulant (MLS-1) | not directly measured | — |

**Current model value:** `γ_sv = 0.040 J/m²`.  
This is **far below** reported values for silicate minerals. The physical
value is likely 0.3–0.5 J/m², but using it directly would make ηₛ = γ_sv + γ_is − γ_iv
very negative (ηₛ ≈ 0.35 + 0.04 − 0.109 = 0.28), which inverts the
sediment double-well and destabilizes the phase field. The current small
value is a numerical choice; see §4 for discussion.

### 1.4 Adopted Model Values and Derived Quantities

The values below are the **currently adopted** defaults in `src/permafrost2.c`
(`gamma_iv`, `gamma_im`, `gamma_mv` in the source). They satisfy all three
Kim-Steinbach positivity constraints simultaneously while preserving a
physically realistic contact angle.

```
γ_iv = 0.109 J/m²     ice–vapor      (Petrenko & Whitworth 1999)
γ_is = 0.057 J/m²     ice–sediment   (Israelachvili 2011, ice on quartz)
γ_sv = 0.080 J/m²     sediment–vapor (effective; keeps all η > 0 at θ = 77.8°)

ηᵢ  = γ_iv + γ_is − γ_sv = 0.109 + 0.057 − 0.080 = +0.086 J/m²  ✓
ηₛ  = γ_sv + γ_is − γ_iv = 0.080 + 0.057 − 0.109 = +0.028 J/m²  ✓
ηₐ  = γ_iv + γ_sv − γ_is = 0.109 + 0.080 − 0.057 = +0.132 J/m²  ✓

Contact angle θ:  cos θ = (γ_sv − γ_is)/γ_iv = (0.080 − 0.057)/0.109 = 0.211
                  θ ≈ 77.8°
```

All three combination energies are positive, ensuring the Kim-Steinbach
double-well is convex for every phase. The contact angle (77.8°) is unchanged
from the previously used values and is physically plausible for ice on dry
lunar silicate.

---

## 2. Thermal Properties

### 2.1 Ice

| Property | Value | Units | Notes | Reference |
|----------|-------|-------|-------|-----------|
| Density ρᵢ | 919 | kg/m³ | hexagonal ice Ih at 273 K | [1] |
| Density ρᵢ (cold) | 927 | kg/m³ | at 100 K (slight increase) | [1] |
| Thermal conductivity kᵢ | 2.22 | W/(m·K) | at 273 K | [1] |
| kᵢ at 200 K | 3.5 | W/(m·K) | T-dependent: kᵢ ∝ 1/T | [1] |
| kᵢ at 100 K | 7.1 | W/(m·K) | — | [1] |
| Specific heat cₚ,ᵢ | 2090 | J/(kg·K) | at 273 K | [1] |
| cₚ,ᵢ at 200 K | 1567 | J/(kg·K) | — | [1] |
| cₚ,ᵢ at 100 K | 800 | J/(kg·K) | — | [1] |
| Latent heat of sublimation Lₛ | 2.83×10⁶ | J/kg | at 273 K | [1] |
| Lₛ at 200 K | 2.87×10⁶ | J/kg | increases slightly at lower T | [10] |

### 2.2 Water Vapor

| Property | Value | Units | Notes | Reference |
|----------|-------|-------|-------|-----------|
| Vapor diffusivity Dᵥ | 2.2×10⁻⁵ | m²/s | in air at 273 K, 1 atm | [11] |
| Dᵥ at 200 K | ~9×10⁻⁶ | m²/s | scales roughly as T^1.75 | [11] |
| Saturation vapor density ρᵥₛ | 4.85×10⁻³ | kg/m³ | at 273 K (ice surface) | [1] |
| ρᵥₛ at 253 K (−20 °C) | 8.49×10⁻⁴ | kg/m³ | Antoine equation / Clausius-Clapeyron | [1] |
| ρᵥₛ at 200 K | ~4×10⁻⁸ | kg/m³ | — | [12] |
| ρᵥₛ at 100 K | ~10⁻²² | kg/m³ | effectively zero | [12] |

**Saturation vapor density** (Clausius-Clapeyron approximation, valid 150–273 K):
```
ρᵥₛ(T) = (M_w / R) × P₀ × exp[−Lₛ M_w / R × (1/T − 1/T₀)] / T
```
where M_w = 0.018 kg/mol, R = 8.314 J/(mol·K), T₀ = 273.15 K,
P₀ = 611 Pa (saturation pressure at 0 °C).

### 2.3 Lunar Regolith

Lunar regolith is a porous, fine-grained (median ~50 µm) mixture of
minerals and glass. Bulk properties vary with density and porosity.

| Property | Value | Units | Notes | Reference |
|----------|-------|-------|-------|-----------|
| Grain density ρₛ | 2900–3100 | kg/m³ | mineral grains | [13] |
| Bulk density (upper layer) | 1300–1800 | kg/m³ | ~40–50% porous | [13] |
| Thermal conductivity kₛ | 7.5×10⁻⁴ | W/(m·K) | at 250 K, bulk (porous) | [14] |
| kₛ (dense mineral) | 1.5–2.5 | W/(m·K) | solid grain, 250 K | [14] |
| kₛ at 100 K (bulk) | ~5×10⁻⁴ | W/(m·K) | colder → less radiation | [14] |
| Specific heat cₚ,ₛ | 840 | J/(kg·K) | at 300 K | [13] |
| cₚ,ₛ at 200 K | 560 | J/(kg·K) | — | [13] |
| cₚ,ₛ at 100 K | 280 | J/(kg·K) | — | [13] |

> **Current model:** `rho_sed = 7753 kg/m³` is a metal-placeholder density.
> For lunar regolith grain density use **2900–3100 kg/m³**; for bulk (porous)
> regolith use **1300–1800 kg/m³** depending on compaction.

---

## 3. Grain Geometry and Characteristic Scales

| Parameter | Symbol | Typical value | Notes |
|-----------|--------|---------------|-------|
| Ice grain radius | Rᵢ | 10–100 µm | PSR ice particle size estimates |
| Regolith grain radius | Rₛ | 20–100 µm | median lunar regolith grain |
| Interface width | ε | ~0.7 × grain radius | PetIGA phase-field parameter |
| Domain size | Lx × Ly | 5–50 grain radii | depends on test |

---

## 4. Discussion: Kim-Steinbach Positivity Constraints

The Kim-Steinbach three-phase model requires all three combination energies to
be positive for the double-wells to be convex:

```
ηᵢ = γ_iv + γ_is − γ_sv > 0    (ice double-well)
ηₛ = γ_sv + γ_is − γ_iv > 0    (sediment double-well)
ηₐ = γ_iv + γ_sv − γ_is > 0    (air double-well, trivially satisfied)
```

These are equivalent to the **triangle inequalities** on {γ_iv, γ_is, γ_sv}.
The admissible window for γ_sv given fixed γ_iv and γ_is is:

```
γ_iv − γ_is  <  γ_sv  <  γ_iv + γ_is
```

With γ_iv = 0.109 J/m² and γ_is = 0.057 J/m² (ice on quartz):
```
0.052  <  γ_sv  <  0.166   [J/m²]
```

The physical surface energy of dry silicate (~0.3–0.5 J/m²) lies above this
window, meaning the Kim-Steinbach model cannot simultaneously use physical
values for all three surfaces for an ice–silicate–vapor system. The adopted
value γ_sv = 0.080 J/m² is within the window and preserves the contact angle
that was used previously.

### Why γ_sv ≪ physical silicate values

If γ_sv were set to its physical value (~0.4 J/m²), ηᵢ = 0.109+0.057−0.4 = −0.234
would be strongly negative, inverting the ice double-well and causing φᵢ to
diverge. The Kim-Steinbach model is designed for symmetric systems (grain
growth in metals) and requires rescaling the surface energies to satisfy the
triangle inequalities. The physically meaningful output of the model is the
**contact angle**, which is preserved exactly at 77.8° with the adopted values.

---

## 5. Summary of Adopted Parameter Values

| Parameter | Adopted value | Physical range | Notes |
|-----------|--------------|----------------|-------|
| γ_iv [J/m²] | 0.109 | 0.060–0.109 | 0.109 at 0°C (Petrenko & Whitworth) |
| γ_is [J/m²] | 0.057 | 0.030–0.057 | Ice on quartz (Israelachvili 2011) |
| γ_sv [J/m²] | 0.080 | 0.052–0.166 (constraint window) | Effective; preserves θ = 77.8°, ηₛ > 0 |
| θ [deg] | 77.8 | — | cos θ = (γ_sv − γ_is)/γ_iv |
| ηᵢ [J/m²] | 0.086 | — | All positive ✓ |
| ηₛ [J/m²] | 0.028 | — | Positive ✓ (was −0.020 with old values) |
| ηₐ [J/m²] | 0.132 | — | Positive ✓ |
| ρᵢ [kg/m³] | 919 | 919–927 | OK; increases slightly at low T |
| ρₛ [kg/m³] | 7753 (placeholder) | 2900–3100 | Update for lunar regolith |
| kᵢ [W/(m·K)] | 2.22 | 2.22–7.1 | Strong T-dependence; use if flag_Tdep |
| kₛ [W/(m·K)] | — | 1.5–2.5 (grain) | Porous bulk ~7.5×10⁻⁴ |
| Dᵥ [m²/s] | 2.2×10⁻⁵ | 9×10⁻⁶–2.2×10⁻⁵ | Reduce for cold PSR conditions |
| Lₛ [J/kg] | 2.83×10⁶ | 2.83–2.87×10⁶ | OK |
| d0_sub0 [m] | 1×10⁻¹¹ | 1×10⁻¹¹ to 1×10⁻¹⁰ | Reduced to suppress AC coarsening |

---

## 6. References

[1] Petrenko, V. F. & Whitworth, R. W. *Physics of Ice.* Oxford University
Press, 1999. (thermal properties, surface energy, saturation vapor pressure)

[2] Lide, D. R. (ed.) *CRC Handbook of Chemistry and Physics*, 85th edn.
CRC Press, 2004. (vapor pressure of ice)

[3] Sega, M., Höfinger, S. & Jedlovszky, P. Pressure profile calculation
with mesh Ewald methods. *J. Chem. Phys.* **126**, 054501 (2007).
(ice surface tension by MD at low T)

[4] Conde, M. M., Vega, C. & Patrykiejew, A. The thickness of a liquid layer
on the free surface of ice as obtained from computer simulation. *J. Chem. Phys.*
**129**, 014702 (2008). (ice surface properties)

[5] Israelachvili, J. N. *Intermolecular and Surface Forces*, 3rd edn.
Academic Press, 2011. (ice–mineral interfacial energies, §17)

[6] Jellinek, H. H. G. Ice adhesion and abhesion: a review. *CRREL Special
Report 83-23*, 1983. (ice on mineral surfaces, contact angle data)

[7] Dash, J. G., Rempel, A. W. & Wettlaufer, J. S. The physics of premelted
ice and its geophysical consequences. *Rev. Mod. Phys.* **78**, 695 (2006).
(ice–mineral interfaces, premelting films)

[8] Parks, G. A. Surface and interfacial free energies of quartz. *J. Geophys.
Res.* **89**, 3997–4008 (1984). (γ_sv for silica surfaces)

[9] Momber, A. W. Quantifying the erosion of rock by high-velocity water jets.
*Wear* **338-339**, 122–130 (2015). (mineral surface energy estimates)

[10] Murphy, D. M. & Koop, T. Review of the vapour pressures of ice and
supercooled water for atmospheric applications. *Q. J. R. Meteorol. Soc.*
**131**, 1539–1565 (2005). (saturation vapor pressure of ice, wide T range)

[11] Massman, W. J. A review of the molecular diffusivities of H₂O, CO₂, CH₄,
CO, O₃, SO₂, NH₃, N₂O, NO and NO₂ in air, O₂ and N₂ near STP. *Atmos.
Environ.* **32**, 1111–1127 (1998). (vapor diffusivity, T-dependence)

[12] Fray, N. & Schmitt, B. Sublimation of ices of astrophysical interest:
a bibliographic review. *Planet. Space Sci.* **57**, 2053–2080 (2009).
(water ice sublimation rates at low T, saturation vapor pressure)

[13] Carrier, W. D., Olhoeft, G. R. & Mendell, W. Physical properties of the
lunar surface. In *Lunar Sourcebook*, Chapter 9. Cambridge University Press,
1991. (lunar regolith density, thermal properties)

[14] Hayne, P. O. et al. Global regolith thermophysical properties of the Moon
from the Diviner Lunar Radiometer Experiment. *J. Geophys. Res. Planets*
**122**, 2371–2400 (2017). (lunar surface thermal conductivity measurements)
