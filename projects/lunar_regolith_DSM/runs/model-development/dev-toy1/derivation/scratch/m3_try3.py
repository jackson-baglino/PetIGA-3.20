"""M3 scratch (try 3): every number quoted in ledger.json / derivation.md, plus
symbolic checks of the try-3 edits (S13 with T_ref, S21-S22 chain rule).

Not evidence for anyone else; the verifier writes its own checks.
"""
import math

import sympy as sp

x, t, xi, tau = sp.symbols('x t xi tau', real=True)
L, k, rho, cp, q, T0, DT = sp.symbols('L k rho c_p q T_0 DeltaT_q', positive=True)
kappa = k / (rho * cp)

# --- S13: storage term from e = rho c_p (T - T_ref) -------------------------
T = sp.Function('T')(x, t)
Tref_t = sp.Function('T_ref')(t)
e_var = rho * cp * (T - Tref_t)
e_const = rho * cp * (T - sp.Symbol('T_ref'))
print('S13 d_t e, T_ref(t)    :', sp.diff(e_var, t))
print('S13 d_t e, T_ref const :', sp.diff(e_const, t))

# --- S21, S22: chain rule on T = T0 + DT * theta(x/L, kappa t / L^2) ---------
th = sp.Function('theta')
Tsc = T0 + DT * th(x / L, kappa * t / L**2)
X, Y = sp.symbols('X Y')
d_tau = sp.Subs(sp.Derivative(th(X, Y), Y), (X, Y), (x / L, kappa * t / L**2))
d_xi = sp.Subs(sp.Derivative(th(X, Y), X), (X, Y), (x / L, kappa * t / L**2))
d_xixi = sp.Subs(sp.Derivative(th(X, Y), X, 2), (X, Y), (x / L, kappa * t / L**2))
print('S21 residual:', sp.simplify(sp.diff(Tsc, t) - kappa * DT / L**2 * d_tau))
print('S22 residual (1st):', sp.simplify(sp.diff(Tsc, x) - DT / L * d_xi))
print('S22 residual (2nd):', sp.simplify(sp.diff(Tsc, x, 2) - DT / L**2 * d_xixi))

# --- S31: reference values (C1 = (a)) ----------------------------------------
ref = {L: 1, k: 1, rho: 1, cp: 1, q: 1, T0: 1}
print('S31 kappa =', kappa.subs(ref), ' t_c =', (L**2 / kappa).subs(ref),
      ' DeltaT_q =', (q * L**2 / k).subs(ref), ' Q =', (q * L**2 / (k * T0)).subs(ref))

# --- S56 numbers ------------------------------------------------------------
pi2 = math.pi**2
g = lambda s: 4 / math.pi * math.exp(-pi2 * s) - 4 / (3 * math.pi) * math.exp(-9 * pi2 * s)
print('S56 ln3 =', math.log(3), ' 8 pi^2 =', 8 * pi2, ' tau_g =', math.log(3) / (8 * pi2))
print('S56 g(0+) = 8/(3 pi) =', 8 / (3 * math.pi))
print('S56 pi^2*0.1 =', pi2 * 0.1, ' exp(-pi^2*0.1) =', math.exp(-pi2 * 0.1),
      ' 4/pi =', 4 / math.pi, ' product =', 4 / math.pi * math.exp(-pi2 * 0.1))
print('S56 9pi^2*0.1 =', 9 * pi2 * 0.1, ' exp(-9pi^2*0.1) =', math.exp(-9 * pi2 * 0.1),
      ' 4/(3pi) =', 4 / (3 * math.pi), ' product =', 4 / (3 * math.pi) * math.exp(-9 * pi2 * 0.1))
print('S56 g(0.1) =', g(0.1))
# g' sign change only at tau_g (sample)
gp = lambda s: -4 * math.pi * math.exp(-pi2 * s) + 12 * math.pi * math.exp(-9 * pi2 * s)
tg = math.log(3) / (8 * pi2)
print('S56 g\' at tau_g/2, tau_g, 2tau_g, 0.1:', gp(tg / 2), gp(tg), gp(2 * tg), gp(0.1))

# --- S57 numbers ------------------------------------------------------------
print('S57 exp(-pi^2) =', math.exp(-pi2), ' 4/pi*exp =', 4 / math.pi * math.exp(-pi2),
      ' 4/pi^3 =', 4 / math.pi**3, ' 4/pi^3*exp =', 4 / math.pi**3 * math.exp(-pi2),
      ' exp(-8pi^2) =', math.exp(-8 * pi2))

# --- S32 / derivation.md section 8 numbers ----------------------------------
print('S32 example tau=0.3: pi^2*0.3 =', pi2 * 0.3, ' bound (4/pi)e^{-pi^2 tau}(1+...) =',
      4 / math.pi * math.exp(-pi2 * 0.3) * (1 + math.exp(-8 * pi2 * 0.3) / (3 * (1 - math.exp(-8 * pi2 * 0.3)))))
vb = 4 / math.pi**3 * math.exp(-pi2)
print('D4 note: transient bound / peak 1/8 =', vb / 0.125)
