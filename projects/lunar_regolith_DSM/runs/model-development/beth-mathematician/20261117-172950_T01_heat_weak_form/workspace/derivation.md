# T01_heat_weak_form: weak form of the 1D transient heat equation

Domain: 0 < x < L. Assumptions used: rho and c are constant; T, w, k and Q are smooth on [0, L].
The test function w is **not** assumed to vanish at x = 0 or x = L, so the boundary terms are kept.

Units: every term of the strong form is in W m^-3. After multiplying by w (units 1) and integrating
over x (units m), every term of the weak form is in W m^-2 (heat flux per unit cross-section area).

## S1 (given): strong form

$$\rho c\,\frac{\partial T}{\partial t} = \frac{\partial}{\partial x}\Big(k\,\frac{\partial T}{\partial x}\Big) + Q$$

Units check: rho c dT/dt -> kg m^-3 · J kg^-1 K^-1 · K s^-1 = W m^-3;
d/dx(k dT/dx) -> m^-1 · W m^-1 K^-1 · K m^-1 = W m^-3; Q -> W m^-3.

## S2 (from S1): move every term to the left side

$$\rho c\,\frac{\partial T}{\partial t} - \frac{\partial}{\partial x}\Big(k\,\frac{\partial T}{\partial x}\Big) - Q = 0$$

## S3 (from S2): multiply by the test function w(x)

$$\rho c\,\frac{\partial T}{\partial t}\,w - \frac{\partial}{\partial x}\Big(k\,\frac{\partial T}{\partial x}\Big)w - Q\,w = 0$$

## S4 (from S3): integrate over (0, L)

$$\int_0^L \rho c\,\frac{\partial T}{\partial t}\,w\,dx - \int_0^L \frac{\partial}{\partial x}\Big(k\,\frac{\partial T}{\partial x}\Big)w\,dx - \int_0^L Q\,w\,dx = 0$$

Units: each term W m^-3 · m = W m^-2.

## S5 (product rule, smoothness assumption): integration-by-parts identity for the conduction term

$$\int_0^L \frac{\partial}{\partial x}\Big(k\,\frac{\partial T}{\partial x}\Big)w\,dx
= \Big[k\,\frac{\partial T}{\partial x}\,w\Big]_{x=0}^{x=L} - \int_0^L k\,\frac{\partial T}{\partial x}\,\frac{dw}{dx}\,dx$$

(From d/dx(k T_x w) = (k T_x)_x w + k T_x w_x, integrated over (0, L) with the fundamental theorem of calculus.)

## S6 (substitute S5 into S4)

$$\int_0^L \rho c\,\frac{\partial T}{\partial t}\,w\,dx
- \Big[k\,\frac{\partial T}{\partial x}\,w\Big]_{0}^{L}
+ \int_0^L k\,\frac{\partial T}{\partial x}\,\frac{dw}{dx}\,dx
- \int_0^L Q\,w\,dx = 0$$

## S7 (from S6): write out the bracket; result R = 0

$$R = \int_0^L \rho c\,\frac{\partial T}{\partial t}\,w\,dx
+ \int_0^L k\,\frac{\partial T}{\partial x}\,\frac{dw}{dx}\,dx
- \int_0^L Q\,w\,dx
- \Big(k\,\frac{\partial T}{\partial x}\,w\Big)\Big|_{x=L}
+ \Big(k\,\frac{\partial T}{\partial x}\,w\Big)\Big|_{x=0} = 0$$

Units: every term is W m^-2. No second derivative of T remains. No approximation was made.

## Check

`scratch/check.py` substitutes two sets of concrete smooth fields (T, w, k, Q) and confirms that R equals
the integral of (strong residual) · w over (0, L) exactly (difference simplifies to 0).
