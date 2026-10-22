# Units of symbols (semantic memory)

> **UNRELIABLE — needs verification.** Exported 2026-10-22 from the
> practice-round registry `Conventions/units.sqlite` (archived in
> `scratch/2026-10-21_practice_round/research_archive/`). Do not cite a unit
> from this file until the units-verification workflow has checked its row
> against the page image. Rows checked against our own source code say so.

What was fixed on 2026-10-22, without the papers:

- Units are written in one notation (`J m^-3 K^-1`). Every unit was checked
  for dimensional consistency with its meaning; none was wrong.
- Five pairs had one merged meaning ("density of ice, dry air" for both ρ_i
  and ρ_a). Each is split by its subscript and marked as inferred.
- `M0` and `alpha` were checked against the solver code.

Known problems still open:

- The `quote` column is the PDF text layer, which drops minus signs, powers
  of ten and micro signs: `1.8106, 1.4103 J m−3 K−1` is the page's
  1.8×10⁶ and 1.4×10³ (values, not units). Never take a number from a quote.
- `verified` says "image-checked (Opus)" for rows whose quotes are still
  mangled; that check covered the unit, not the quote.
- `Pa` (Kaempfer and Plapp) and `P` (Moure and Fu) are the same quantity;
  the symbol `Pa` also collides with the unit.
- `T_hot`/`T_cold`: which printed symbol (T0, T7) is which is not recorded.

| symbol | meaning | unit | source | page | quote (text layer) | verified |
|---|---|---|---|---|---|---|
| `M0` | mobility of the phase-field equation (mob_sub); same role and unit as M_j of Moure and Fu | m s^-1 | dimensional-analysis |  |  | checked against src/monitoring.c:58 (mob = eps/(3 tau_sub) → m s^-1) |
| `alpha` | sublimation attachment rate in the ice source term, alpha = lam/tau_sub (Moure and Fu alpha_j; SI Table S1) | s^-1 | user |  |  | checked against src/monitoring.c:59 (alph = lambda/tau_sub, lambda dimensionless → s^-1) |
| `Ca` | heat capacity per volume of dry air (printed as Ca) — meaning inferred from the subscript | J m^-3 K^-1 | kaempfer_plapp_2009 | 4 | 1.8106, 1.4103 J m−3 K−1 | image-checked (Opus) against the page image |
| `Ci` | heat capacity per volume of ice (printed as Ci) — meaning inferred from the subscript | J m^-3 K^-1 | kaempfer_plapp_2009 | 4 | 1.8106, 1.4103 J m−3 K−1 | image-checked (Opus) against the page image |
| `D_STP` | water vapor diffusion coefficient in air at standard pressure and temperature (printed as DSTP) | m^2 s^-1 | kaempfer_plapp_2009 | 4 | 2.17810−5 m2 s−1 | image-checked (Opus) against the page image |
| `D_v` | water vapor diffusion coefficient in air (printed as Dv) | m^2 s^-1 | kaempfer_plapp_2009 | 4 | 2.17810−5 m2 s−1 | image-checked (Opus) against the page image |
| `L_sg` | latent heat of sublimation of ice per volume (printed as Lsg) | J m^-3 | kaempfer_plapp_2009 | 4 | 2.60109 J m−3 26 | image-checked (Opus) against the page image |
| `Pa` | atmospheric pressure (printed as Pa) | Pa | kaempfer_plapp_2009 | 4 | 1.01325105 Pa | image-checked (Opus) against the page image |
| `R_da` | individual gas constant of dry air (printed as Rda) — meaning inferred from the subscript | J kg^-1 K^-1 | kaempfer_plapp_2009 | 4 | 286.9, 461.5 J Kg−1 K−1 | image-checked (Opus) against the page image |
| `R_v` | individual gas constant of water vapor (printed as Rv) — meaning inferred from the subscript | J kg^-1 K^-1 | kaempfer_plapp_2009 | 4 | 286.9, 461.5 J Kg−1 K−1 | image-checked (Opus) against the page image |
| `T0` | reference temperature (printed as T0) | K | kaempfer_plapp_2009 | 4 | reference temperature of T0=−10 °C, respectively, 263 K | image-checked (Opus) against the page image |
| `T_STP` | standard temperature (printed as TSTP) | K | kaempfer_plapp_2009 | 4 | 273.15 K | image-checked (Opus) against the page image |
| `a` | mean intermolecular spacing in ice (printed as a) | m | kaempfer_plapp_2009 | 4 | 3.1910−10 m | image-checked (Opus) against the page image |
| `beta` | interface kinetic coefficient (printed as β) | s m^-1 | kaempfer_plapp_2009 | 4 | in simulations =5.5105 s m−1 | image-checked (Opus) against the page image |
| `d0` | capillary length (printed as d0) | m | kaempfer_plapp_2009 | 4 | 1.310−9 m | image-checked (Opus) against the page image |
| `gamma` | interfacial free energy of ice (printed as γ) | J m^-2 | kaempfer_plapp_2009 | 4 | 1.09  10−1 J m−2 26 | image-checked (Opus) against the page image |
| `k` | Boltzmann's constant (printed as k) | J K^-1 | kaempfer_plapp_2009 | 4 | 1.3810−23 J K−1 | image-checked (Opus) against the page image |
| `kappa_a` | heat conductivity of dry air (printed as κa) — meaning inferred from the subscript | W m^-1 K^-1 | kaempfer_plapp_2009 | 4 | 2.29, 0.02 W m−1 K−1 47 | image-checked (Opus) against the page image |
| `kappa_i` | heat conductivity of ice (printed as κi) — meaning inferred from the subscript | W m^-1 K^-1 | kaempfer_plapp_2009 | 4 | 2.29, 0.02 W m−1 K−1 47 | image-checked (Opus) against the page image |
| `m` | mass of a water molecule (printed as m) | kg | kaempfer_plapp_2009 | 4 | 0.018/6.021023=2.9920−26 kg | image-checked (Opus) against the page image |
| `rho_a` | density of dry air (printed as ρa) — meaning inferred from the subscript | kg m^-3 | kaempfer_plapp_2009 | 4 | 918.9, 1.341 kg m−3 47 | image-checked (Opus) against the page image |
| `rho_i` | density of ice (printed as ρi) — meaning inferred from the subscript | kg m^-3 | kaempfer_plapp_2009 | 4 | 918.9, 1.341 kg m−3 47 | image-checked (Opus) against the page image |
| `rho_v` | water vapor concentration (density) (printed as ρv) | kg m^-3 | kaempfer_plapp_2009 | 4 | 10−3 to 10−2 kg m−3 24 | image-checked (Opus) against the page image |
| `v_n` | normal interface velocity (printed as vn) | m s^-1 | kaempfer_plapp_2009 | 4 | 10−9 m s−1 | image-checked (Opus) against the page image |
| `x_s` | specific humidity ratio of air at saturation (printed as xs) | 1 (dimensionless) | kaempfer_plapp_2009 | 4 | kg/kg | image-checked (Opus) against the page image |
| `L` | estimated pore size (printed as L) | m | kaempfer_plapp_2009 | 9 | estimated pore size of L=10−4 m | image-checked (Opus) against the page image |
| `W` | interface thickness (printed as W) | m | kaempfer_plapp_2009 | 9 | W510−7 m to W510−5 m | image-checked (Opus) against the page image |
| `t_D` | characteristic water vapor diffusion time (printed as tD) | s | kaempfer_plapp_2009 | 9 | tD=L2/Dv 510−4 s | image-checked (Opus) against the page image |
| `t_n` | characteristic interface migration time (printed as tn) | s | kaempfer_plapp_2009 | 9 | tn=L/vn3105 s | image-checked (Opus) against the page image |
| `Deltat` | time step (printed as Δt) | s | kaempfer_plapp_2009 | 11 | t = 5  10−5 s | image-checked (Opus) against the page image |
| `T_cold` | colder boundary temperature (the paper prints T0 and T7; which symbol is which is not established) (printed as T_cold) — meaning inferred from the subscript | K | kaempfer_plapp_2009 | 11 | T0=261 K, T7=260 K | image-checked (Opus) against the page image |
| `T_hot` | warmer boundary temperature (the paper prints T0 and T7; which symbol is which is not established) (printed as T_hot) — meaning inferred from the subscript | K | kaempfer_plapp_2009 | 11 | T0=261 K, T7=260 K | image-checked (Opus) against the page image |
| `j_T` | heat flux (printed as jT) | W m^-2 | kaempfer_plapp_2009 | 11 | =19.091Wm−2 | image-checked (Opus) against the page image |
| `j_v` | water vapor mass flux (printed as jv) | kg m^-2 s^-1 | kaempfer_plapp_2009 | 12 | j [10−6kgm−2s−1] v | image-checked (Opus) against the page image |
| `t` | time (printed as t) | s | kaempfer_plapp_2009 | 12 | time [s] | image-checked (Opus) against the page image |
| `h` | spatial grid resolution (printed as h) | m | kaempfer_plapp_2009 | 13 | resolution h = 5  10−7 m | image-checked (Opus) against the page image |
| `P` | atmospheric pressure (printed as P) | Pa | moure_fu_2024 | 3 | P is the atmospheric pressure (units: Pa) | image-checked (Opus) against the page image |
| `eps` | phase-field interface width parameter (printed as ε) | m | moure_fu_2024 | 6 | ε < 10−8 m, this severely restricts the domain size | image-checked (Opus) against the page image |
| `M_j` | phase-field mobility (printed as Mj) | m s^-1 | moure_fu_2024_SI | 7 | Mj (m/s) | image-checked (Opus) against the page image |
| `alpha_j` | phase-change coefficient (printed as αj) | s^-1 | moure_fu_2024_SI | 7 | αj (1/s) | image-checked (Opus) against the page image |
| `beta_sol` | kinetic attachment coefficient, solidification (printed as βsol) | s m^-1 | moure_fu_2024_SI | 7 | βj (s/m) | image-checked (Opus) against the page image |
| `d_j` | capillary length (Gibbs-Thomson), j = sol, sub, eva (printed as dj) | m | moure_fu_2024_SI | 7 | dj (m) | image-checked (Opus) against the page image |
| `u_j` | phase-change factor multiplying αjφl²φm² (printed as uj) | 1 (dimensionless) | moure_fu_2024_SI | 7 | uj (-) | image-checked (Opus) against the page image |
| `v_n_j` | normal interface velocity (printed as vn,j) | m s^-1 | moure_fu_2024_SI | 7 | vn,j (m/s) | image-checked (Opus) against the page image |
| `chi` | interface curvature (printed as χ) | m^-1 | moure_fu_2024_SI | 8 | typical ice curvatures χ ∼104 m−1 | image-checked (Opus) against the page image |
| `D_v0` | vapor diffusivity in air (printed as Dv0) | m^2 s^-1 | moure_fu_2024_SI | 9 | Dv0 2.178 × 10−5 m2/s | image-checked (Opus) against the page image |
| `K_a` | thermal conductivity of air (printed as Ka) | W m^-1 K^-1 | moure_fu_2024_SI | 9 | Ka 0.02 W/(m ◦C) | image-checked (Opus) against the page image |
| `K_i` | thermal conductivity of ice (printed as Ki) | W m^-1 K^-1 | moure_fu_2024_SI | 9 | Ki 2.29 W/(m ◦C) | image-checked (Opus) against the page image |
| `K_w` | thermal conductivity of water (printed as Kw) | W m^-1 K^-1 | moure_fu_2024_SI | 9 | Kw 0.554 W/(m ◦C) | image-checked (Opus) against the page image |
| `L_sol` | latent heat of solidification (printed as Lsol) | J kg^-1 | moure_fu_2024_SI | 9 | Lsol 3.34 × 105 J/Kg | image-checked (Opus) against the page image |
| `L_sub` | latent heat of sublimation (printed as Lsub) | J kg^-1 | moure_fu_2024_SI | 9 | Lsub 2.83 × 106 J/Kg | image-checked (Opus) against the page image |
| `Lambda` | model parameter Λ (printed as Λ) | N m^-1 | moure_fu_2024_SI | 9 | Λ 1 N/m | image-checked (Opus) against the page image |
| `alpha_nucl` | nucleation coefficient (printed as αnucl) | m^-1 | moure_fu_2024_SI | 9 | αnucl 5 × 106 1/m | image-checked (Opus) against the page image |
| `c_p_a` | specific heat of air (printed as cp,a) | J kg^-1 K^-1 | moure_fu_2024_SI | 9 | cp,a 1.044 × 103 J/(Kg ◦C) | image-checked (Opus) against the page image |
| `c_p_i` | specific heat of ice (printed as cp,i) | J kg^-1 K^-1 | moure_fu_2024_SI | 9 | cp,i 1.96 × 103 J/(Kg ◦C) | image-checked (Opus) against the page image |
| `c_p_w` | specific heat of water (printed as cp,w) | J kg^-1 K^-1 | moure_fu_2024_SI | 9 | cp,w 4.2 × 103 J/(Kg ◦C) | image-checked (Opus) against the page image |
| `gamma_ia` | ice-air surface tension (printed as σia) | N m^-1 | moure_fu_2024_SI | 9 | σia 0.109 N/m | image-checked (Opus) against the page image |
| `gamma_iw` | ice-water surface tension (printed as σiw) | N m^-1 | moure_fu_2024_SI | 9 | σiw 0.033 N/m | image-checked (Opus) against the page image |
| `gamma_wa` | water-air surface tension (printed as σwa) | N m^-1 | moure_fu_2024_SI | 9 | σwa 0.076 N/m | image-checked (Opus) against the page image |
| `rho_w` | water density (printed as ρw) | kg m^-3 | moure_fu_2024_SI | 9 | ρw 1000 Kg/m3 | image-checked (Opus) against the page image |
| `xi_T` | temperature time-scale factor (printed as ξT) | 1 (dimensionless) | moure_fu_2024_SI | 9 | ξT see Eq. (13) (-) | image-checked (Opus) against the page image |
| `xi_v` | vapor time-scale factor (printed as ξv) | 1 (dimensionless) | moure_fu_2024_SI | 9 | ξv 10−3 (-) | image-checked (Opus) against the page image |
