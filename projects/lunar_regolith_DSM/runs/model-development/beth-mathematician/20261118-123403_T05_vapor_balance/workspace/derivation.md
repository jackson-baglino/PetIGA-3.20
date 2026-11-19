# T05 – Vapor mass balance

New symbols: rho_i (ice density, kg m^-3, constant, ~920), D(T) (diffusivity of water vapor in the pore gas, m^2 s^-1, ~1e-5 near 200 K at about 1 bar; much larger at low pressure).

**S1 (assumption, mixture mass).** Each control volume holds ice with volume fraction phi and gas with volume fraction 1-phi. The water mass per unit total volume is

$$ m = \rho_i\,\phi + (1-\phi)\,\rho_v \qquad [\mathrm{kg\,m^{-3}}] $$

**S2 (conservation principle).** For any fixed control volume V with no sources of water:
$$ \frac{d}{dt}\int_V m\,dV = -\oint_{\partial V} \mathbf J\cdot\mathbf n\,dA $$
Divergence theorem and the fact that V is arbitrary (fields smooth) give
$$ \partial_t m + \nabla\cdot\mathbf J = 0 . $$

**S3 (constitutive flux).** Ice does not move (no flux from ice). Vapor moves only by Fickian diffusion through the gas fraction:
$$ \mathbf J = -(1-\phi)\,D(T)\,\nabla\rho_v \qquad [\mathrm{kg\,m^{-2}\,s^{-1}}] $$

**S4 (substitute S1, S3 into S2; rho_i constant so rho_i passes through d/dt).**
$$ \partial_t\big[(1-\phi)\rho_v\big] + \rho_i\,\partial_t\phi - \nabla\cdot\big[(1-\phi)D(T)\nabla\rho_v\big] = 0 $$
The term rho_i d_t phi is the exchange with ice: d_t phi < 0 (sublimation) is a source of vapor.
Checked: `scratch/check.py` (balance − R = 0, SymPy); units of every term kg m^-3 s^-1 (pint).

**S5 (product rule on S4, for reading).**
$$ (1-\phi)\,\partial_t\rho_v + (\rho_i-\rho_v)\,\partial_t\phi - \nabla\cdot\big[(1-\phi)D\nabla\rho_v\big] = 0 $$
Checked: equal to S4 (SymPy). The factor (rho_i − rho_v) is the usual Stefan-type jump in density.

**S6 (Cartesian form).**
$$ R = \partial_t[(1-\phi)\rho_v] + \rho_i\partial_t\phi - \sum_{\xi\in\{x,y,z\}} \partial_\xi\big[(1-\phi)D(T)\,\partial_\xi\rho_v\big] = 0 $$

## Assumptions
- phi, T known and smooth (given).
- Ice is rigid and stationary; ice density rho_i constant (needed in S4).
- Vapor transport is Fickian diffusion in the gas, flux weighted by gas fraction (1-phi); D may depend on T.
- No bulk gas flow; no other sources/sinks of water.

## Neglected
- Thermal expansion of ice: volumetric coefficient ~1.5e-4 K^-1, so d ln rho_i over tens of K is <1%.
- Advection by bulk gas flow / Stefan flow: vapor mole fraction at 200 K is ~1e-6 to 1e-5 of ~1 bar air, so the induced velocity is negligible.
- Adsorption on regolith grains and grain motion: outside this pore-space model.
- Soret (thermal) diffusion: thermal diffusion factor for water in air is ~0.1 or less, a small correction compared with the concentration gradient driven flux.
- Knudsen effects: absorbed into D(T); fine for this form, though D would then also depend on pore size.
