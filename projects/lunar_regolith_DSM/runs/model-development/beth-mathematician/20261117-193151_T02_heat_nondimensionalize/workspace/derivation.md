# T02: Nondimensionalizing the 1D heat equation

Scales: x = L ξ, t = τ s with τ = ρ c L²/k0 (units: kg m⁻³ · J kg⁻¹ K⁻¹ · m² / (W m⁻¹ K⁻¹) = s), T = T0 + ΔT θ(ξ, s).

**S1 (given).** Every term has units of W m⁻³:

$$\rho c\,\partial_t T = k_0\,\partial_x^2 T + Q_0$$

**S2 (chain rule for t, with ∂s/∂t = 1/τ and T0 constant).**

$$\partial_t T = \frac{\Delta T}{\tau}\,\partial_s\theta = \frac{k_0\,\Delta T}{\rho c L^2}\,\partial_s\theta$$

**S3 (chain rule for x applied twice, with ∂ξ/∂x = 1/L).**

$$\partial_x^2 T = \frac{\Delta T}{L^2}\,\partial_\xi^2\theta$$

**S4 (substitute S2 and S3 into S1).** Both sides have units of W m⁻³:

$$\frac{k_0\,\Delta T}{L^2}\,\partial_s\theta = \frac{k_0\,\Delta T}{L^2}\,\partial_\xi^2\theta + Q_0$$

**S5 (divide by k0 ΔT/L², which has units of W m⁻³ and is positive).** Every term is dimensionless:

$$\partial_s\theta = \partial_\xi^2\theta + \frac{Q_0 L^2}{k_0\,\Delta T}$$

Group check: (W m⁻³)(m²)/((W m⁻¹ K⁻¹)(K)) = 1.

**S6 (move all terms to the left side).**

$$R \equiv \partial_s\theta - \partial_\xi^2\theta - \frac{Q_0 L^2}{k_0\,\Delta T} = 0$$

The result is exact. Nothing is approximated. I checked it with SymPy in `scratch/check.py`: the difference is 0.
