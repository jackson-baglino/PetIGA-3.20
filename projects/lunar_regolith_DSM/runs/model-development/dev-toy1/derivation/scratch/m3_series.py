"""M3 scratch (try 2): checks for the new steps S48-S57 and the restated S32.

Not evidence for anyone else; the verifier writes its own checks.
"""
import math

import sympy as sp

xi, tau = sp.symbols('xi tau', real=True)
l = sp.symbols('ell', integer=True, positive=True)

# S52: the two integrals
I1 = sp.simplify(sp.integrate(xi * sp.sin(l * sp.pi * xi), (xi, 0, 1)))
I2 = sp.simplify(sp.integrate(xi**2 * sp.sin(l * sp.pi * xi), (xi, 0, 1)))
I1_claim = -(-1)**l / (l * sp.pi)
I2_claim = -(-1)**l / (l * sp.pi) + 2 * ((-1)**l - 1) / (l * sp.pi)**3
print('S52 I1 ok:', sp.simplify(I1 - I1_claim) == 0, ' I2 ok:', sp.simplify(I2 - I2_claim) == 0)

# S53: beta_l = 2 int theta_s sin = I1 - I2
thetas = xi * (1 - xi) / 2
beta = sp.simplify(2 * sp.integrate(thetas * sp.sin(l * sp.pi * xi), (xi, 0, 1)))
beta_claim = 2 * (1 - (-1)**l) / (l * sp.pi)**3
print('S53 beta ok:', sp.simplify(beta - beta_claim) == 0, ' beta(1)=', sp.nsimplify(beta.subs(l, 1)), ' beta(2)=', beta.subs(l, 2))

# S51: v = thetas - theta solves v_tau = v_xixi given theta_tau = theta_xixi + 1, thetas'' = -1
print('S51 thetas\'\' =', sp.diff(thetas, xi, 2))

# S54/S55: each term solves the heat equation; d/dtau of beta_l e^{-l^2 pi^2 tau} sin = -(4/(l pi)) ... for odd l
term = 4 / (l * sp.pi)**3 * sp.sin(l * sp.pi * xi) * sp.exp(-l**2 * sp.pi**2 * tau)
print('S54 term heat eq:', sp.simplify(sp.diff(term, tau) - sp.diff(term, xi, 2)) == 0)
print('S55 -d term/dtau:', sp.simplify(-sp.diff(term, tau)))

# numbers
def dtheta_dtau(x, t, n=4001):
    return sum(4 / (k * math.pi) * math.sin(k * math.pi * x) * math.exp(-(k * math.pi)**2 * t)
               for k in range(1, n, 2))

def v_val(x, t, n=4001):
    return sum(4 / (k * math.pi)**3 * math.sin(k * math.pi * x) * math.exp(-(k * math.pi)**2 * t)
               for k in range(1, n, 2))

pi2 = math.pi**2
g = lambda t: 4 / math.pi * math.exp(-pi2 * t) - 4 / (3 * math.pi) * math.exp(-9 * pi2 * t)
tc = math.log(3) / (8 * pi2)
print('S56: tau_c = ln3/(8pi^2) =', tc, ' g(0+) = 8/(3pi) =', 8 / (3 * math.pi), ' g(tau_c)=', g(tc), ' g(0.1) =', g(0.1))
# lower bound vs true value on (0, 0.1]
worst = min(dtheta_dtau(0.5, t) - g(t) for t in [i * 1e-3 for i in range(1, 101)])
mn = min(g(t) for t in [i * 1e-4 for i in range(1, 1001)])
print('S56: min_{(0,0.1]} g =', mn, ' min(true - bound) on grid =', worst)
print('S56: true dtheta/dtau(1/2, 0.1) =', dtheta_dtau(0.5, 0.1))

br = lambda t: 1 + math.exp(-8 * pi2 * t) / (3 * (1 - math.exp(-8 * pi2 * t)))
ub = 4 / math.pi * math.exp(-pi2) * br(1.0)
print('S57: storage bound at tau=1:', ub, ' series max over xi:', max(abs(dtheta_dtau(i / 200, 1.0)) for i in range(201)))
ubv = 4 / math.pi**3 * math.exp(-pi2) * (1 + math.exp(-8 * pi2) / (27 * (1 - math.exp(-8 * pi2))))
print('S57: |theta-thetas| bound at tau=1:', ubv, ' series max:', max(abs(v_val(i / 200, 1.0)) for i in range(201)))
print('S57: e^{-8 pi^2} =', math.exp(-8 * pi2))
# series check of theta >= 0 (S50) on a grid
print('S50: min theta on grid:', min(0.5 * x * (1 - x) - v_val(x, t, 801)
                                     for x in [i / 100 for i in range(101)] for t in (1e-3, 1e-2, 0.1, 1.0)))
