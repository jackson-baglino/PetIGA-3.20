# T01: Weak form of the 1D transient heat equation

Domain 0 < x < L. Units: every term in the strong equation is in W m^-3. After multiplying by w (units 1)
and integrating over x (units m), every term is in W m^-2. Boundary terms w k dT/dx are also W m^-2.

**S1 (given).** Strong form:

$$\rho c\,\partial_t T = \partial_x\!\left(k\,\partial_x T\right) + Q \qquad [\mathrm{W\,m^{-3}}]$$

**S2 (move all terms to one side, from S1).**

$$\rho c\,\partial_t T - \partial_x\!\left(k\,\partial_x T\right) - Q = 0$$

**S3 (multiply S2 by w(x)).**

$$w\,\rho c\,\partial_t T - w\,\partial_x\!\left(k\,\partial_x T\right) - wQ = 0$$

**S4 (integrate S3 over (0, L)).**

$$\int_0^L \rho c\, w\,\partial_t T\,dx - \int_0^L w\,\partial_x\!\left(k\,\partial_x T\right)dx - \int_0^L wQ\,dx = 0 \qquad [\mathrm{W\,m^{-2}}]$$

**S5 (product rule).** Since T, w, k are smooth on [0, L]:

$$w\,\partial_x\!\left(k\,\partial_x T\right) = \partial_x\!\left(w\,k\,\partial_x T\right) - k\,\partial_x w\,\partial_x T$$

**S6 (integrate S5 over (0, L), fundamental theorem of calculus = integration by parts).**
w is not assumed to vanish at either end, so both boundary terms stay:

$$\int_0^L w\,\partial_x\!\left(k\,\partial_x T\right)dx = \Big[w\,k\,\partial_x T\Big]_{x=L} - \Big[w\,k\,\partial_x T\Big]_{x=0} - \int_0^L k\,\partial_x w\,\partial_x T\,dx$$

**S7 (substitute S6 into S4).** No second derivative of T remains:

$$R \equiv \int_0^L \rho c\, w\,\partial_t T\,dx + \int_0^L k\,\partial_x w\,\partial_x T\,dx - \int_0^L wQ\,dx - \Big[w\,k\,\partial_x T\Big]_{x=L} + \Big[w\,k\,\partial_x T\Big]_{x=0} = 0$$

Units check: each term is W m^-2.

Sign convention: R is defined as the integral of w times (rho c T_t - (k T_x)_x - Q), i.e. the
storage term carries a plus sign. The negative of R is an equally valid residual.

SymPy check (`scratch/check.py`): the product-rule identity in S5 simplifies to 0, and a concrete
test with T = x^2 t + sin(x) t^2, w = 1 + x + x^3, k = 2 + x^2, Q = cos(x) gives the same value
for S4 and S7.
