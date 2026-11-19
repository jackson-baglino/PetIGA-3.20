"""Verifier's independent reference of spec.md (heat1d-fd-be), for calibrating the M5 tests only.

DEV TOY DRY RUN. Written by the verifier from spec.md alone; it is NOT the M6 implementation and must not
be copied into the feature worktree. It lives in the run folder so the modeler's M6 work stays independent.
"""
import math

import numpy as np

_EPS = np.finfo(float).eps


def _scales(params):
    L, k, rho, cp, q, T0 = (float(params[s]) for s in ("L", "k", "rho", "cp", "q", "T0"))
    kappa = k / (rho * cp)
    Q = q * L * L / (k * T0)
    return kappa, Q, L * L / kappa, k * T0 / (L * L)


def mms_forcing(xi, tau, A, B, Q):
    e = np.exp(-tau)
    return (A * np.sin(np.pi * xi) * (e + np.pi ** 2 * (1 - e))
            + 2 * np.pi * B * xi * (1 - xi) * np.cos(2 * np.pi * tau)
            + 2 * B * np.sin(2 * np.pi * tau) - Q)


def mms_exact(xi, tau, A, B):
    return 1 + A * np.sin(np.pi * xi) * (1 - np.exp(-tau)) + B * xi * (1 - xi) * np.sin(2 * np.pi * tau)


def residual(tau, U, Udot, n, params):
    _, Q, _, _ = _scales(params)
    h = 1.0 / n
    Ue = np.concatenate(([1.0], np.asarray(U, float), [1.0]))
    F = np.asarray(Udot, float) - (Ue[:-2] - 2 * Ue[1:-1] + Ue[2:]) / h ** 2 - Q
    m = params.get("mms")
    if m is not None:
        xi = np.arange(1, n) * h
        F = F - mms_forcing(xi, tau, m["A"], m["B"], Q)
    return F


def jacobian_tridiag(shift, n, params):
    h = 1.0 / n
    m = n - 1
    return np.full(m, -1 / h ** 2), np.full(m, shift + 2 / h ** 2), np.full(m, -1 / h ** 2)


def thomas_solve(sub, diag, sup, rhs):
    m = len(diag)
    c = np.empty(m)
    d = np.empty(m)
    c[0] = sup[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, m):
        den = diag[i] - sub[i] * c[i - 1]
        c[i] = sup[i] / den
        d[i] = (rhs[i] - sub[i] * d[i - 1]) / den
    x = np.empty(m)
    x[-1] = d[-1]
    for i in range(m - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


def _guard(R, shift, h, U):
    lim = 1e3 * _EPS / 2 * (shift + 4 / h ** 2) * max(1.0, np.max(np.abs(U)))
    if np.max(np.abs(R)) > lim:
        raise RuntimeError("residual guard: %g > %g" % (np.max(np.abs(R)), lim))


def step_be(U, tau_new, dtau, n, params):
    shift = 1.0 / dtau
    R = residual(tau_new, U, np.zeros_like(U), n, params)
    sub, diag, sup = jacobian_tridiag(shift, n, params)
    Un = U - thomas_solve(sub, diag, sup, R)
    _guard(residual(tau_new, Un, (Un - U) / dtau, n, params), shift, 1.0 / n, Un)
    return Un


def solve(n, dt, t_end, params):
    kappa, Q, tc, _ = _scales(params)
    L, T0 = float(params["L"]), float(params["T0"])
    r = t_end / dt
    ns = int(round(r))
    if ns < 0 or abs(r - ns) > 1e-10 * max(1.0, abs(r)):
        raise ValueError("t_end/dt = %r is not an integer" % r)
    dtau = kappa * dt / L ** 2
    U = np.ones(n - 1)
    for s in range(ns):
        U = step_be(U, (s + 1) * dtau, dtau, n, params)
    return L * np.arange(n + 1) / n, T0 * np.concatenate(([1.0], U, [1.0]))


def solve_steady(n, params):
    p = dict(params)
    p["mms"] = None
    L, T0 = float(params["L"]), float(params["T0"])
    U0 = np.ones(n - 1)
    R = residual(0.0, U0, np.zeros(n - 1), n, p)
    sub, diag, sup = jacobian_tridiag(0.0, n, p)
    U = U0 - thomas_solve(sub, diag, sup, R)
    _guard(residual(0.0, U, np.zeros(n - 1), n, p), 0.0, 1.0 / n, U)
    return L * np.arange(n + 1) / n, T0 * np.concatenate(([1.0], U, [1.0]))


def end_fluxes(U, U_prev, dtau, h, Q, f_end):
    f0, f1 = f_end
    F0 = (U[1] - U[0]) / h + h / 2 * (Q + f0) - h / 2 * (U[0] - U_prev[0]) / dtau
    F1 = (U[-2] - U[-1]) / h + h / 2 * (Q + f1) - h / 2 * (U[-1] - U_prev[-1]) / dtau
    return F0, F1


def dt_for_tolerance(eps_theta_prime, n, t_end, params):
    kappa, Q, _, _ = _scales(params)
    L = float(params["L"])
    tau_end = kappa * t_end / L ** 2
    lam1 = 4 * n * n * math.sin(math.pi / (2 * n)) ** 2
    b1 = 4 / math.pi ** 3
    dmax = tau_end if Q == 0 else 2 * math.e * eps_theta_prime / (abs(Q) * b1 * lam1)
    ns = math.ceil(tau_end / dmax)
    return tau_end / ns * L ** 2 / kappa


def check_jacobian(n, params, U=None, Udot=None, shifts=(0, 1e1, 1e3, 1e5), eps=1e-6, seed=0):
    rng = np.random.default_rng(seed)
    m = n - 1
    U = 1 + 0.1 * rng.standard_normal(m) if U is None else np.asarray(U, float)
    Udot = rng.standard_normal(m) if Udot is None else np.asarray(Udot, float)
    tau = 0.3
    out = {"groups": {}}
    worst = 0.0
    for shift in shifts:
        sub, diag, sup = (np.asarray(a, float) for a in jacobian_tridiag(shift, n, params))
        J = np.diag(diag) + np.diag(sub[1:], -1) + np.diag(sup[:-1], 1)
        err = 0.0
        for k in range(5):
            v = rng.standard_normal(m)
            fd = (residual(tau, U + eps * v, Udot + shift * eps * v, n, params)
                  - residual(tau, U - eps * v, Udot - shift * eps * v, n, params)) / (2 * eps)
            err = max(err, np.linalg.norm(J @ v - fd) / max(np.linalg.norm(fd), 1e-300))
        worst = max(worst, err)
        out["groups"]["shift=%g" % shift] = err
    out["max_rel_err"] = worst
    out["pass"] = bool(worst <= 1e-6)
    return out
