# calonneThermalConductivitySnow2019 — Calonne, Milliancourt, Burr, Philip, Martin, Flin & Geindreau, GRL 46, 13079–13089 (2019)
DOI (PDF p.1–2): 10.1029/2019GL085228. HAL version, read in full (12 pp: HAL cover + 11 journal pages).

## Summary
Effective conductivity tensors computed by homogenization (Auriault et al. 2009 cell problem, Geodict finite differences) on 64 3D µCT images: 37 seasonal snow, 27 Antarctic firn/porous ice; 102–888 kg/m³; at −3, −20, −60 °C with temperature-dependent k_ice, k_air. Quadratic k(ρ) for snow, linear for firn above ~450 kg/m³; anisotropy; comparison with formulas; new density–temperature formula (Eq. 5).

## Findings relevant to the manuscript
- Why conductivity matters: "Estimating thermal conductivity of snow, firn, and porous ice is key for modeling the thermal regime of alpine and polar glaciers" (abstract, p.2/13079); temperature propagation through snow/firn "is ... crucial to understand" for many applications (p.2).
- Measured vs effective: "experimental values of thermal conductivity often correspond to apparent thermal conductivities, which characterize all the heat transport processes involved and not purely conduction. Overcoming the experimental issues, numerical simulations on 3-D tomographic images ... allowed to retrieve the effective thermal conductivity tensor (Calonne et al., 2011; Kaempfer et al., 2005)" (p.3/13080). Directly supports the outline's "conductivity here is conductive only".
- Method: "A specific boundary value problem arising from a homogenization technique (Auriault et al., 2009; Calonne et al., 2014) is solved on representative elementary volumes" (p.4/13081).
- Temperature-dependent phase conductivities: k_ice = 2.107, 2.330, 2.900 and k_air = 0.024, 0.023, 0.019 W/m/K at −3, −20, −60 °C (p.4/13081). Effect on bulk k: for firn >600 kg/m³ ~10% from −3 to −20 °C; "Those variations decrease to about 20% and 5% at 300 kg/m³ ... and become negligible at lower densities" (p.5/13082).
- Bonds matter at low density: "At low densities, ice grains are insulated by the predominant surrounding pore space. The few bonds between ice grains play then a critical role in the overall heat conduction as they act as heat paths ... the portion of ice added or removed elsewhere than at ice bonds has little impact" (p.5/13082). Strongly relevant: physical reason why neck growth (sintering) raises k at fixed density — supports the manuscript's main result qualitatively.
- Anisotropy: "anisotropy of the thermal conductivity of snow is overall well correlated with the anisotropy of the correlation length: k ranges from 0.66 for a rounded grains sample with a lc of 0.67, to 1.57 for a depth hoar sample with a lc of 1.34" (p.6–7/13083–13084); firn/porous ice: anisotropy "tend[s] to vanish above 600 kg/m³", highest firn ratio 1.11 (p.6/13083). Supports outline L551–553.
- Snow formulas scatter: "A large scatter was however often reported when plotting conductivity versus density, which was attributed to additional influences of the snow microstructure other than density" (p.2–3).
- Off-diagonal terms negligible; kx ≃ ky (p.4/13081).

## Claims supported / contradicted / qualified
- Outline L294–300 (Intro para 2 addition: k_eff sets how fast heat crosses snow/firn and the temperature profile; Calonne 2011, 2019): supported by p.2/13079.
- Outline L551–553 (k_z/k_xy tracks lc ratio: 0.66 at 0.67; 1.57 at 1.34; vanishes in firn above ~600): supported verbatim (p.6–7).
- Outline L1355–1356 (Limitations: Calonne 2019 on temperature-dependent conductivities): supported; the manuscript holds k_ice = 2.29, k_air = 0.02 constant, whereas Calonne use k_ice from 2.107 (−3 °C) to 2.330 (−20 °C). Over −40 to −5 °C k_ice changes by roughly 15–20%; this is absent from the model and bears on the "temperature only rescales time" claim (it holds by construction for conduction because k_ice, k_air are fixed). The limitation should say so.
- Outline L1288 (Calonne 2011, 2019 for k vs density as the contrast): supported.
- Draft L256 cites only Calonne 2011 for heat flow importance; 2019 is closer for "thermal regime".

## Candidate additions
- Section 5.3 ("Density is not enough"): Calonne 2019 p.5 — at low density the bonds are the heat paths, so ice added away from bonds barely changes k. This is the published mechanism behind conductivity rising while density stays fixed in sintering. Closely matching.
- Section 2.3: Calonne 2019 p.3 for apparent vs effective conductivity.
- Limitations: temperature dependence of k_ice/k_air (p.4–5) is not in the model; with it, curves at different T would not collapse exactly in absolute k (normalized curves are less affected).
