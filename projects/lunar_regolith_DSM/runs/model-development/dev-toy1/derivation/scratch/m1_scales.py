"""M1 scratch: reference-case numbers for the toy 1D heat problem (proposed C1 values)."""
import math

L, k, rho, cp, q, T0 = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  # proposed unit reference case (C1)
kappa = k / (rho * cp)
t_diff = L**2 / kappa
t_slow = L**2 / (math.pi**2 * kappa)       # e-folding time of slowest mode sin(pi x/L)
dT_q = q * L**2 / k                        # source temperature scale
Q = q * L**2 / (k * T0)                    # the one dimensionless group
dT_max = q * L**2 / (8 * k)                # peak of steady parabola at x = L/2
t_99 = t_slow * math.log(100)              # time for slowest mode to decay to 1 %
print(f"kappa={kappa} t_diff={t_diff} t_slow={t_slow:.6f} t_99={t_99:.6f}")
print(f"dT_q={dT_q} Q={Q} dT_max={dT_max} dT_max/T0={dT_max/T0}")
