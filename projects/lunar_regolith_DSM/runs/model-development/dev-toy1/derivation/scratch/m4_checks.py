"""M4 scratch: numbers quoted in spec.md. My own checks, not evidence for anyone else.

Prints:
  1. a priori time scales for the reference case (C1, N = 64);
  2. the discrete eigenvalues lambda_1 and lambda_max;
  3. the dt-rule value for the reference tolerance, and the reverse (tolerance of dtau = 1e-3);
  4. measured max-over-time BE error (mode-1 prediction vs a fine-N, many-mode series) for a few dtau;
  5. a Jacobian check of the residual of spec.md (J v vs central difference), per row group;
  6. that the residual with the MMS forcing vanishes to O(h^2 + dtau) on the manufactured solution.
"""
import math
import sys

import numpy as np

# ---------------- 1. scales (C1: L = k = rho = c_p = q = T_0 = 1) ----------------
L, k, rho, cp, q, T0 = 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
kappa = k / (rho * cp)
N = 64
h = 1.0 / N                    # dimensionless mesh size
t_c = L**2 / kappa
t_slow = L**2 / (math.pi**2 * kappa)
t_mesh = (L / N) ** 2 / kappa
Q = q * L**2 / (k * T0)
print("1. kappa=%.6g t_c=%.6g t_slow=%.6g t_mesh=(L/N)^2/kappa=%.6g Q=%.6g" % (kappa, t_c, t_slow, t_mesh, Q))

# ---------------- 2. discrete eigenvalues ----------------
def lam(l, N):
    hh = 1.0 / N
    return 4.0 / hh**2 * math.sin(l * math.pi * hh / 2) ** 2

lam1 = lam(1, N)
lammax = lam(N - 1, N)
print("2. lambda_1=%.6f (pi^2=%.6f) lambda_max=%.6f  1/lambda_max=%.6e  h^2/4=%.6e" % (lam1, math.pi**2, lammax, 1 / lammax, h * h / 4))

# ---------------- 3. dt rule ----------------
b1 = 4.0 / math.pi**3
def dtau_rule(eps_theta, lam1):
    return 2.0 * math.e * eps_theta / (b1 * lam1)
def eps_of_dtau(dtau, lam1):
    return b1 * lam1 * dtau / (2.0 * math.e)
print("3. b1=%.6f  eps(dtau=1e-3)=%.6e  dtau(eps=1e-4)=%.6e  dtau(eps=2.34e-4)=%.6e" % (
    b1, eps_of_dtau(1e-3, lam1), dtau_rule(1e-4, lam1), dtau_rule(2.34e-4, lam1)))

# ---------------- helpers: residual / Jacobian of spec.md (theta' form) ----------------
def residual(U, Uold, dtau, h, Q, gL=1.0, gR=1.0, f=None):
    """R_j = (U_j - Uold_j)/dtau - (U_{j-1} - 2U_j + U_{j+1})/h^2 - Q - f_j, j = 1..N-1,
    with U_0 = gL, U_N = gR. dtau = inf means steady (no storage)."""
    n = U.size
    Up = np.concatenate(([gL], U, [gR]))
    lap = (Up[:-2] - 2 * Up[1:-1] + Up[2:]) / h**2
    store = (U - Uold) / dtau if np.isfinite(dtau) else 0.0 * U
    R = store - lap - Q
    if f is not None:
        R = R - f
    return R

def jac_dense(n, dtau, h):
    shift = 1.0 / dtau if np.isfinite(dtau) else 0.0
    J = np.zeros((n, n))
    for j in range(n):
        J[j, j] = shift + 2.0 / h**2
        if j > 0:
            J[j, j - 1] = -1.0 / h**2
        if j < n - 1:
            J[j, j + 1] = -1.0 / h**2
    return J

# ---------------- 4. measured BE time error vs prediction ----------------
def series_theta(xi, tau, K=4001):
    s = xi * (1 - xi) / 2
    for l in range(1, K + 1, 2):
        s = s - 4.0 / (l * math.pi) ** 3 * np.sin(l * math.pi * xi) * math.exp(-(l * math.pi) ** 2 * tau)
    return s

def discrete_series(N, tau, dtau=None):
    """Exact solution of the semi-discrete (dtau=None) or fully discrete BE problem in theta form,
    via the discrete sine eigenbasis. Used only to separate time error from space error."""
    hh = 1.0 / N
    j = np.arange(1, N)
    xi = j * hh
    # discrete steady state = theta_s at nodes (selection.md D2)
    ths = xi * (1 - xi) / 2
    # expand ths in discrete sine basis
    out = ths.copy()
    for l in range(1, N):
        phi = np.sin(l * math.pi * xi)
        c = 2 * hh * np.dot(ths, phi)
        lm = lam(l, N)
        if dtau is None:
            amp = math.exp(-lm * tau)
        else:
            nsteps = round(tau / dtau)
            amp = (1.0 + lm * dtau) ** (-nsteps)
        out = out - c * amp * phi
    return xi, out

print("4. time error of BE alone (fully discrete minus semi-discrete, same N=64), max over a tau grid:")
for dtau in (4e-3, 2e-3, 1e-3, 5e-4):
    worst, tworst = 0.0, None
    taus = [dtau * m for m in range(1, int(round(1.0 / dtau)) + 1)]
    for tau in taus:
        _, a = discrete_series(N, tau)
        _, b = discrete_series(N, tau, dtau)
        e = np.max(np.abs(a - b))
        if e > worst:
            worst, tworst = e, tau
    print("   dtau=%.1e  max_err=%.4e at tau=%.4e   mode-1 prediction=%.4e  ratio=%.3f" % (
        dtau, worst, tworst, eps_of_dtau(dtau, lam1), worst / eps_of_dtau(dtau, lam1)))

# ---------------- 5. Jacobian check ----------------
rng = np.random.default_rng(12345)
n = N - 1
dtau = 1e-3
U = 1.0 + 0.1 * rng.standard_normal(n)
Uold = 1.0 + 0.1 * rng.standard_normal(n)
J = jac_dense(n, dtau, h)
worst_int, worst_bnd = 0.0, 0.0
for trial in range(5):
    v = rng.uniform(-1, 1, n)
    for eps in (1e-6,):
        fd = (residual(U + eps * v, Uold, dtau, h, Q) - residual(U - eps * v, Uold, dtau, h, Q)) / (2 * eps)
        Jv = J @ v
        rows_b = [0, n - 1]
        rows_i = list(range(1, n - 1))
        eb = np.linalg.norm(fd[rows_b] - Jv[rows_b]) / np.linalg.norm(Jv[rows_b])
        ei = np.linalg.norm(fd[rows_i] - Jv[rows_i]) / np.linalg.norm(Jv[rows_i])
        worst_int, worst_bnd = max(worst_int, ei), max(worst_bnd, eb)
print("5. Jacobian check (dtau=1e-3, eps=1e-6): worst rel err interior rows=%.3e, boundary rows=%.3e" % (worst_int, worst_bnd))
# deliberately wrong Jacobian (shift term dropped): must be caught
Jbad = jac_dense(n, np.inf, h)
v = rng.uniform(-1, 1, n)
fd = (residual(U + 1e-6 * v, Uold, dtau, h, Q) - residual(U - 1e-6 * v, Uold, dtau, h, Q)) / 2e-6
print("   with shift term dropped: rel err=%.3e (must be large)" % (np.linalg.norm(fd - Jbad @ v) / np.linalg.norm(Jbad @ v)))

# ---------------- 6. MMS consistency ----------------
# theta'_m(xi,tau) = 1 + A sin(pi xi)(1 - exp(-tau)) + B xi(1-xi) sin(2 pi tau) ; f = d_tau u - d_xixi u - Q
A, B = 0.5, 0.25
def um(xi, tau):
    return 1 + A * np.sin(np.pi * xi) * (1 - np.exp(-tau)) + B * xi * (1 - xi) * np.sin(2 * np.pi * tau)
def fm(xi, tau, Q):
    ut = A * np.sin(np.pi * xi) * np.exp(-tau) + 2 * np.pi * B * xi * (1 - xi) * np.cos(2 * np.pi * tau)
    uxx = -A * np.pi**2 * np.sin(np.pi * xi) * (1 - np.exp(-tau)) - 2 * B * np.sin(2 * np.pi * tau)
    return ut - uxx - Q
print("6. MMS truncation (residual of the exact MMS field, max over rows) at tau=0.5:")
for NN, dt in ((32, 1e-3), (64, 1e-3), (128, 1e-3), (64, 4e-3), (64, 2e-3)):
    hh = 1.0 / NN
    xi = np.arange(1, NN) * hh
    tau1 = 0.5
    U1 = um(xi, tau1)
    U0 = um(xi, tau1 - dt)
    R = residual(U1, U0, dt, hh, Q, gL=um(0.0, tau1), gR=um(1.0, tau1), f=fm(xi, tau1, Q))
    print("   N=%4d dtau=%.0e  max|R|=%.4e" % (NN, dt, np.max(np.abs(R))))

# ---------------- 7. full MMS run: observed orders ----------------
def solve_mms(NN, dt, tau_end):
    hh = 1.0 / NN
    xi = np.arange(1, NN) * hh
    U = um(xi, 0.0)
    Jm = jac_dense(NN - 1, dt, hh)
    nst = int(round(tau_end / dt))
    for s in range(1, nst + 1):
        tau1 = s * dt
        R = residual(U, U, dt, hh, Q, gL=um(0.0, tau1), gR=um(1.0, tau1), f=fm(xi, tau1, Q))
        # one Newton step from U (problem is affine): U1 = U - J^{-1} R(U)
        U = U - np.linalg.solve(Jm, R)
    return np.max(np.abs(U - um(xi, tau_end)))

print("7. MMS run, error at tau_end=0.5:")
e1 = solve_mms(16, 1e-5, 0.5); e2 = solve_mms(32, 1e-5, 0.5)
print("   space (dtau=1e-5): N=16 %.4e N=32 %.4e order %.3f" % (e1, e2, math.log2(e1 / e2)))
e3 = solve_mms(256, 1e-2, 0.5); e4 = solve_mms(256, 5e-3, 0.5)
print("   time  (N=256):     dtau=1e-2 %.4e 5e-3 %.4e order %.3f" % (e3, e4, math.log2(e3 / e4)))
e5 = solve_mms(64, 1e-3, 0.5)
print("   reference (N=64, dtau=1e-3): %.4e" % e5)

# ---------------- 8. F(U,Udot) form Jacobian check, split into dF/dU and dF/dUdot ----------------
def F(U, Udot, h, Q, gL=1.0, gR=1.0, f=None):
    Up = np.concatenate(([gL], U, [gR]))
    lap = (Up[:-2] - 2 * Up[1:-1] + Up[2:]) / h**2
    R = Udot - lap - Q
    return R if f is None else R - f
Ud = rng.standard_normal(n)
print("8. F(U,Udot) Jacobian check, eps=1e-6:")
for shift in (0.0, 1e1, 1e3, 1e5):
    Jm = jac_dense(n, 1.0 / shift if shift > 0 else np.inf, h)
    v = rng.uniform(-1, 1, n)
    e = 1e-6
    full = (F(U + e * v, Ud + shift * e * v, h, Q) - F(U - e * v, Ud - shift * e * v, h, Q)) / (2 * e)
    dU = (F(U + e * v, Ud, h, Q) - F(U - e * v, Ud, h, Q)) / (2 * e)
    dUd = (F(U, Ud + e * v, h, Q) - F(U, Ud - e * v, h, Q)) / (2 * e)
    Jv = Jm @ v
    JU = jac_dense(n, np.inf, h) @ v
    print("   shift=%.0e  full rel=%.2e  dF/dU rel=%.2e  dF/dUdot (vs v) rel=%.2e" % (
        shift, np.linalg.norm(full - Jv) / np.linalg.norm(Jv), np.linalg.norm(dU - JU) / np.linalg.norm(JU),
        np.linalg.norm(dUd - v) / np.linalg.norm(v)))

# ---------------- 9. steady residual and tau_end comparison ----------------
xi = np.arange(1, N) * h
ths_p = 1 + Q * xi * (1 - xi) / 2
print("9. steady residual on theta'_s nodes: max|R|=%.3e" % np.max(np.abs(F(ths_p, 0 * ths_p, h, Q))))
_, a = discrete_series(N, 1.0, 1e-3)
print("   fully discrete (N=64,dtau=1e-3) max|theta - theta_s| at tau=1: %.4e" % np.max(np.abs(a - xi * (1 - xi) / 2)))
# end flux, half-cell, at steady state (theta' form): (U_1 - U_0)/h + h Q/2
print("   half-cell end flux at steady: %.15f (Q/2=%.15f)" % ((ths_p[0] - 1.0) / h + h * Q / 2, Q / 2))
sys.exit(0)
