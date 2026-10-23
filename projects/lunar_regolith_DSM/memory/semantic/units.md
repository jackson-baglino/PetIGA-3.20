# Units of symbols (semantic memory)

Checked by the units-verification workflow, run `20261022-155335` (2026-10-22; record in
`runs/units-verification/20261022-155335/`). Each literature row's unit, printed symbol,
page and value were read blind from the page image (Sonnet), compared by script with
the table and with an OCR reading (Marker), disputed rows were judged against the
image (Opus), every meaning was reviewed against the page (Opus), and the user
inspected every disputed row and a random 10 % of the rest.

- 54 literature rows verified; 0 unresolved.
- Removed at the user's decision: t_n, T_cold, T_hot, beta_sol; earlier, `j_T` and `j_v`.
- `value (printed)` is the value as printed on the cited page, transcribed from the
  image; empty where the page prints none (model unknowns such as `rho_v`).
- `M0` and `alpha` are the solver's own parameters, checked against the code, not
  part of this run.

Notes (from `workflows/units-verification/table_notes.md`):

- `Pa` (Kaempfer and Plapp) and `P_a` (Moure and Fu) are the same quantity.
- `D_v_STP` (Kaempfer and Plapp, printed D_v^STP) and `D_v0` (Moure and Fu) are the
  same parameter.
- `T0` here is Moure and Fu's 273.15 offset (T_K = T_C + T_0). Kaempfer and Plapp's
  T_0 is a simulation's reference temperature, and the Moure and Fu SI also writes
  T_0(x) for an initial temperature field.

| symbol | meaning | unit | value (printed) | source | page | verified |
|---|---|---|---|---|---|---|
| `M0` | mobility of the phase-field equation (mob_sub); same role and unit as M_j of Moure and Fu | m s^-1 |  | dimensional-analysis |  | checked against src/monitoring.c:58 (mob = eps/(3 tau_sub) → m s^-1) |
| `alpha` | sublimation attachment rate in the ice source term, alpha = lam/tau_sub (Moure and Fu alpha_j; SI Table S1) | s^-1 |  | user |  | checked against src/monitoring.c:59 (alph = lambda/tau_sub, lambda dimensionless → s^-1) |
| `Ca` | heat capacity per volume of dry air (printed as Ca) | J m^-3 K^-1 | 1.4 × 10^3 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `Ci` | heat capacity per volume of ice (printed as Ci) | J m^-3 K^-1 | 1.8 × 10^6 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `D_v_STP` | water vapor diffusion coefficient in air at standard pressure and temperature; the same parameter as D_v0 of Moure and Fu (printed as D_v^STP) | m^2 s^-1 | 2.178 × 10^-5 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: corrected by the user (U1), then re-check agrees: K&P print D_v^STP; it is the same parameter as Moure and Fu's D_v0. Values and units found are correct. |
| `D_v` | water vapor diffusion coefficient in air (printed as Dv) | m^2 s^-1 | 2.178 × 10^-5 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `L_sg` | latent heat of sublimation of ice per volume (printed as Lsg) | J m^-3 | 2.60 × 10^9 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `Pa` | atmospheric pressure (printed as Pa) | Pa | 1.01325 × 10^5 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `R_da` | individual gas constant of dry air (printed as Rda) | J kg^-1 K^-1 | 286.9 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `R_v` | individual gas constant of water vapor (printed as Rv) | J kg^-1 K^-1 | 461.5 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `T0` | reference temperature converting °C to K, T_K = T_C + T_0 (printed as T0) | °C | 273.15 | moure_fu_2024 | 3 | verified 2026-10-22, run 20261022-155335: corrected by the user (U1), then re-check agrees: T_0 is overloaded. Use Moure and Fu's T_0, which converts °C to K (T_K = T_C + T_0); the K&P T_0 is a simulation's initial condition. Unit written as °C because T is in °C (°C and K differ only in where zero is). |
| `T_STP` | standard temperature (printed as TSTP) | K | 273.15 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `a` | mean intermolecular spacing in ice (printed as a) | m | 3.19 × 10^-10 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `beta` | interface kinetic coefficient (printed as β) | s m^-1 | in simulations β=5.5 × 10^5 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `d0` | capillary length (printed as d0) | m | 1.3 × 10^-9 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `gamma` | interfacial free energy of ice (printed as γ) | J m^-2 | 1.09 × 10^-1 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `k` | Boltzmann's constant (printed as k) | J K^-1 | 1.38 × 10^-23 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `kappa_a` | heat conductivity of dry air (printed as κa) | W m^-1 K^-1 | 0.02 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `kappa_i` | heat conductivity of ice (printed as κi) | W m^-1 K^-1 | 2.29 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `m` | mass of a water molecule (printed as m) | kg | 0.018/6.02 × 10^23 = 2.99 × 20^-26 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `rho_a` | density of dry air (printed as ρa) | kg m^-3 | 1.341 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `rho_i` | density of ice (printed as ρi) | kg m^-3 | 918.9 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `rho_v` | water vapor concentration (density) (printed as ρv) | kg m^-3 |  | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Vapor density, a model unknown, so no value is expected; the unit is correct. |
| `v_n` | normal interface velocity (printed as vn) | m s^-1 | ~10^-9 | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Correct unit; ~10^-9 is fine: it results from the time-dependent phase-field solution. |
| `x_s` | specific humidity ratio of air at saturation (printed as xs) | 1 (dimensionless) |  | kaempfer_plapp_2009 | 4 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Specific humidity ratio, dimensionless (kg/kg); correct. |
| `L` | estimated pore size (printed as L) | m | 10^-4 | kaempfer_plapp_2009 | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Correct. The value is chosen by the user; in practice of order mm (not stated in the papers). |
| `W` | interface thickness (printed as W) | m | ≪ 5 × 10^-7 m to ≪ 5 × 10^-5 m | kaempfer_plapp_2009 | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `t_D` | characteristic water vapor diffusion time (printed as tD) | s | L^2/D_v ≅ 5 × 10^-4 | kaempfer_plapp_2009 | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `Deltat` | time step (printed as Δt) | s | ξΔt=5 × 10^-5 (and ξΔt=1 × 10^-5 for finest discretization) | kaempfer_plapp_2009 | 11 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `t` | time (printed as t) | s | 300 | kaempfer_plapp_2009 | 12 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Time of the problem; t = 300 s only sets the time of K&P's figure. |
| `h` | spatial grid resolution (printed as h) | m | 5 × 10^-7 | kaempfer_plapp_2009 | 13 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `P_a` | atmospheric pressure (printed as Pa) | Pa | 101325 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: corrected by the user (U1), then re-check agrees: Moure and Fu write P_a for atmospheric pressure; its value, P_a = 101325, is in SI Table S3. |
| `eps` | phase-field interface width parameter (printed as ε) | m | < 10^-6 (simulation choice); most restrictive condition ε < 10^-8 | moure_fu_2024 | 6 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `M_j` | phase-field mobility (printed as Mj) | m s^-1 | 2.12 × 10^-5 (Solidification), 4.33 × 10^-7 (Sublimation), 1.34 × 10^-6 (Evaporation) | moure_fu_2024_SI | 7 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): M_j is the mobility, defined in Eq. (16) of Moure and Fu; the SI defines M_sol, M_sub and M_eva, listed as the M_j values in Table S1. |
| `alpha_j` | phase-change coefficient (printed as αj) | s^-1 | 7.96 × 10^4 (Solidification), 1.3 × 10^7 (Sublimation), 6.69 × 10^7 (Evaporation) | moure_fu_2024_SI | 7 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `d_j` | capillary length (Gibbs-Thomson), j = sol, sub, eva (printed as dj) | m | Literature: ~ 5 × 10^-10 (Solidification), ~ 10^-9 (Sublimation), ~ 5 × 10^-10 (Evaporation); this work: 4 × 10^-9, 10^-7, 6 × 10^-8 | moure_fu_2024_SI | 7 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `u_j` | phase-change factor multiplying αjφl²φm² (printed as uj) | 1 (dimensionless) | ~ 10^-3 (Solidification), ~ 10^-8 (Sublimation), ~ 10^-8 (Evaporation) | moure_fu_2024_SI | 7 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Correct; not used in this project. |
| `v_n_j` | normal interface velocity (printed as vn,j) | m s^-1 | ~ (10^-4, 10^-5) (Solidification), ~ (10^-7, 10^-9) (Sublimation), ~ (10^-6, 10^-8) (Evaporation) | moure_fu_2024_SI | 7 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `chi` | interface curvature (printed as χ) | m^-1 | ∼ 10^4 | moure_fu_2024_SI | 8 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `D_v0` | vapor diffusivity in air (printed as Dv0) | m^2 s^-1 | 2.178 × 10^-5 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct; the S4m difference (at the freezing point) is minimal and can be ignored here. |
| `K_a` | thermal conductivity of air (printed as Ka) | W m^-1 K^-1 | 0.02 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `K_i` | thermal conductivity of ice (printed as Ki) | W m^-1 K^-1 | 2.29 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `K_w` | thermal conductivity of water (printed as Kw) | W m^-1 K^-1 | 0.554 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `L_sol` | latent heat of solidification (printed as Lsol) | J kg^-1 | 3.34 × 10^5 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `L_sub` | latent heat of sublimation (printed as Lsub) | J kg^-1 | 2.83 × 10^6 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `Lambda` | model parameter Λ (printed as Λ) | N m^-1 | 1 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `alpha_nucl` | nucleation coefficient (printed as αnucl) | m^-1 | 5 × 10^6 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): Correct; not used in this project. |
| `c_p_a` | specific heat of air (printed as cp,a) | J kg^-1 K^-1 | 1.044 × 10^3 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `c_p_i` | specific heat of ice (printed as cp,i) | J kg^-1 K^-1 | 1.96 × 10^3 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: approved by the user (U1): The blind read is correct. |
| `c_p_w` | specific heat of water (printed as cp,w) | J kg^-1 K^-1 | 4.2 × 10^3 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `gamma_ia` | ice-air surface tension (printed as σia) | N m^-1 | 0.109 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `gamma_iw` | ice-water surface tension (printed as σiw) | N m^-1 | 0.033 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `gamma_wa` | water-air surface tension (printed as σwa) | N m^-1 | 0.076 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `rho_w` | water density (printed as ρw) | kg m^-3 | 1000 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `xi_T` | temperature time-scale factor (printed as ξT) | 1 (dimensionless) | see Eq. (13) | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
| `xi_v` | vapor time-scale factor (printed as ξv) | 1 (dimensionless) | 10^-3 | moure_fu_2024_SI | 9 | verified 2026-10-22, run 20261022-155335: unit, symbol and page agree in table, blind read and Marker; meaning checked against the page |
