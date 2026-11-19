# T04_heat_latent: energy equation with latent heat

Checks: `scratch/check.py` (SymPy identity for S6 to S7, residual 0; pint units check for S9; size estimates for dropped terms).

**Assumptions.** (A1) phi is known and smooth. (A2) Ice is rigid and at rest, and its density rho_i is constant (needed in S7 to take rho_i out of the time derivative). (A3) Gas is at rest on average. Vapor moves only by a diffusive mass flux j (kg m^-2 s^-1), and air has no net flux. (A4) Local thermal equilibrium: one T for ice and gas. (A5) Low pressure, so enthalpy equals internal energy per unit mass up to p/rho, and pressure work is dropped.

**S1 (definition).** Enthalpy per unit volume (J m^-3):

    H = phi rho_i h_i(T) + (1-phi) rho_v h_v(T) + (1-phi) rho_a h_a(T)

**S2 (definition).** Energy flux (W m^-2): Fourier conduction plus the enthalpy carried by diffusing vapor:

    q = -k(phi,T) grad T + h_v(T) j

**S3 (energy conservation over any fixed control volume, then the divergence theorem):**

    dH/dt + div q = 0

**S4 (water mass conservation, ice plus vapor):**

    d/dt[phi rho_i + (1-phi) rho_v] + div j = 0

**S5 (air mass conservation, A3):**

    d/dt[(1-phi) rho_a] = 0

**S6 (subtract h_v times S4 and h_a times S5 from S3; both are zero):**

    dH/dt + div q - h_v(T)*S4 - h_a(T)*S5 = 0

**S7 (product and chain rule, using dh_k/dt = c_k dT/dt, with c_k = dh_k/dT and A2; SymPy residual 0):**

    C(phi,T) dT/dt - div(k grad T) - rho_i L(T) dphi/dt + c_v j . grad T = 0
    C = phi rho_i c_i + (1-phi)(rho_v c_v + rho_a c_a),   L = h_v - h_i  (latent heat of sublimation)

**S8 (drop the vapor enthalpy transport c_v j . grad T).** Relative to conduction it is about c_v D rho_v / k ~ 1850*1e-3*1.7e-6/2 ~ 2e-6.

    C(phi,T) dT/dt - div(k(phi,T) grad T) - rho_i L(T) dphi/dt = 0

**S9 (write it in 3D Cartesian coordinates as R = 0):**

    R = C dT/dt - d/dx(k dT/dx) - d/dy(k dT/dy) - d/dz(k dT/dz) - rho_i L dphi/dt

Units check (pint): every term is in W m^-3 (kg m^-1 s^-3). Sign check: during deposition dphi/dt > 0, the source term +rho_i L dphi/dt is positive, and T rises.
