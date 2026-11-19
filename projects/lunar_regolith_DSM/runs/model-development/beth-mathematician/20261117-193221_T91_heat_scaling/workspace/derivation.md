# T91_heat_scaling: nondimensionalizing the 1D heat equation

**Status: blocked.** The units the work order gives for Q0 do not match the given equation.

## Units check of the given equation (step S1)

$$\rho c\,\partial_t T = k_0\,\partial_x^2 T + Q_0$$

- $\rho c\,\partial_t T$: (kg m^-3)(J kg^-1 K^-1)(K s^-1) = **W m^-3**
- $k_0\,\partial_x^2 T$: (W m^-1 K^-1)(K m^-2) = **W m^-3**
- $Q_0$: the work order gives **W m^-2**

The source term is off by one factor of length. The equation is consistent only if Q0 is a
volumetric source (W m^-3), or if the equation is missing a 1/length factor (for example, Q0/L or
Q0/h for some thickness h).

## Conditional derivation (holds only if Q0 is in W m^-3)

S2. Chain rule with $x = L\xi$, $t = \tau s$, $\tau = \rho c L^2/k_0$ (units s), $T = T_0 + \Delta T\,\theta$:

$$\partial_t T = \frac{\Delta T}{\tau}\,\partial_s\theta, \qquad \partial_x^2 T = \frac{\Delta T}{L^2}\,\partial_\xi^2\theta$$

S3. Substitute S2 into the given equation:

$$\rho c\,\frac{\Delta T\,k_0}{\rho c L^2}\,\partial_s\theta = k_0\,\frac{\Delta T}{L^2}\,\partial_\xi^2\theta + Q_0$$

S4. Simplify the left coefficient:

$$\frac{k_0\Delta T}{L^2}\,\partial_s\theta = \frac{k_0\Delta T}{L^2}\,\partial_\xi^2\theta + Q_0$$

S5. Divide by $k_0\Delta T/L^2$ (units W m^-3) and move all terms to one side:

$$R = \partial_s\theta - \partial_\xi^2\theta - \frac{Q_0 L^2}{k_0\,\Delta T} = 0$$

With Q0 in W m^-3, the group $Q_0L^2/(k_0\Delta T)$ has units 1. With Q0 in W m^-2 as the work
order states, it has units of m. In that case it is not dimensionless, and the task's demand that every
group be dimensionless cannot be met.

SymPy check: `scratch/check.py` confirms the algebra in S2 to S5 (the residual against a test
function is 0) and the units of the group in both cases.
