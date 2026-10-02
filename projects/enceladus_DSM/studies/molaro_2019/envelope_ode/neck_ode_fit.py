#!/usr/bin/env python3
"""neck_ode_fit.py -- back-of-the-envelope neck-growth ODE: what alpha_c and
local undersaturation would reproduce Molaro et al. (2019) at -20 C?

Run from enceladus_DSM/:
    python studies/molaro_2019/envelope_ode/neck_ode_fit.py

THE MODEL (vapour route only, quasi-steady, two-sphere geometry)
----------------------------------------------------------------
Neck radius x, grain radius R, Kuczynski fillet rho = x^2 / (2(R - x)).
The neck surface (mean curvature 1/x - 1/rho, concave) is fed by vapour from
the convex grain surface (2/R) through two resistances in series, attachment
and gas diffusion over a length ~ c*rho:

    v_n   = [ d0 (1/rho - 1/x + 2/R) - s ] / ( beta_sub + c * rho * k_D )
    dx/dt = (pi/2) v_n           (small-neck: dV/dx = 2 pi x^3/R, A = pi^2 x^3/R)

    beta_sub = rho_i / (alpha_c rho_vs v_HK)   K&P beta_0 [s/m], prop. to 1/alpha_c
    k_D      = rho_i / (D_v rho_vs)            gas-diffusion resistance per metre
    s        = local undersaturation at the neck (0 = saturated)

Limits: attachment-limited (beta_sub >> c rho k_D) gives x^3 ~ t (a = 1/3);
diffusion-limited gives x^5 ~ t (a = 1/5); s > 0 caps the neck where
d0/rho falls to s, which bends the late-time slope down further.

Grain shrinkage uses the same physics with the far-field undersaturation s_inf
and diffusion length R:  dR/dt = -(s_inf + 2 d0/R) / (beta_sub + R k_D).

PROCEDURE
---------
1. Calibrate the one geometric unknown, c, on OUR -20 C Fig. 2 run (alpha_c = 0.1
   known), fitting c and s to its neck_width.csv.
2. Validate on the mesh_pair fine run (alpha_c = 1e-3, sealed h = 1, tangent
   start): predict time to the 32.81 um anchor and the post-anchor exponent.
3. Fit the Molaro -20 C DATA for alpha_c and s with c fixed. Profile alpha_c
   over [1e-3, 1] and report the best achievable fit at each value.
4. Read s_inf from the large grain's shrinkage, as a consistency bound on s.
Outputs: results/ (CSV + figure).
"""

import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]                                    # enceladus_DSM/
sys.path.insert(0, str(ROOT / "preprocess"))
import comp_eps as ce                                     # noqa: E402

RES = HERE / "results"
RHO_I = 917.0
T_C = -20.0
R_EFF = 2 * 72.5e-6 * 101e-6 / (72.5e-6 + 101e-6)      # 84.4 um, Molaro pair
FIG2_RUN = Path("~/SimulationResults/HPC_results/enceladus_DSM/GrainPairSintering/"
                "batch_2026-09-08__17.20.46_molaro_T-20_round2/"
                "molaro_2D_L450x225um_eps0.12um_axisym_T-20pair_r14um_dom2__"
                "molaro_T-20_h0.99715_2h_a1e-1_dirichlet/neck_width.csv").expanduser()
DATA = ROOT / "inputs/validation/molaro2019_fig11_T-20.csv"


def thermo(T):
    rvs = ce.rho_vs_sat(T)
    return dict(d0=ce.capillary_length(T), k_D=RHO_I / (ce.Dv_T(T) * rvs),
                bsub_per_inv_alpha=ce.beta_HK(T, 1.0) * RHO_I / rvs)


TH = thermo(T_C)


def beta_sub(alpha):
    return TH["bsub_per_inv_alpha"] / alpha


def rhs(t, y, alpha, s, c, R):
    x = y[0]
    rho = x * x / (2.0 * (R - x))
    drive = TH["d0"] * (1.0 / rho - 1.0 / x + 2.0 / R) - s
    return [0.5 * math.pi * drive / (beta_sub(alpha) + c * rho * TH["k_D"])]


def integrate(x0, t_eval, alpha, s, c, R=R_EFF):
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), [x0], t_eval=t_eval,
                    args=(alpha, s, c, R), method="LSODA", rtol=1e-8, atol=1e-12)
    if not sol.success or sol.y.shape[1] != len(t_eval):
        return np.full(len(t_eval), np.nan)
    return sol.y[0]


def d_free(t, w):
    """C (t + t0)^a with t0 free, the protocol used throughout the repo."""
    from scipy.optimize import curve_fit
    f = lambda t, C, t0, a: C * (t + t0) ** a
    p, cov = curve_fit(f, t, w, p0=[w[0] / 10, 60.0, 0.2],
                       bounds=([0, 1e-6, 0.01], [1, 1e7, 2]), maxfev=50000)
    return p[2], 1.96 * math.sqrt(cov[2, 2])


def load_data():
    d = np.loadtxt(DATA, delimiter=",", comments="#")
    return d[:, 0] * 60.0, d[:, 1] * 1e-6 / 2, d[:, 4] * 1e-6 / 2   # t[s], x, R_large


def main():
    RES.mkdir(exist_ok=True)
    rows = []
    print(f"-20 C: d0={TH['d0']:.4e} m  k_D={TH['k_D']:.3e} s/m^2  "
          f"beta_sub(alpha=0.1)={beta_sub(0.1):.3e} s/m  R_eff={R_EFF*1e6:.1f} um")

    # ---- 1. calibrate c on our Fig. 2 run (alpha = 0.1 known) -------------
    sim = np.loadtxt(FIG2_RUN, delimiter=",", skiprows=1)
    ts, xs = sim[:, 0], sim[:, 1] / 2
    m = ts >= 60.0                               # skip the IC fillet transient
    ts, xs = ts[m] - ts[m][0], xs[m]
    sel = np.unique(np.linspace(0, len(ts) - 1, 80).astype(int))
    ts, xs = ts[sel], xs[sel]

    def res_sim(p):
        c, s = math.exp(p[0]), p[1]
        return (integrate(xs[0], ts, 0.1, s, c) - xs) / xs
    fit = least_squares(res_sim, [0.0, 1e-3], bounds=([-5, 0], [5, 5e-2]))
    c_cal, s_sim = math.exp(fit.x[0]), fit.x[1]
    rms = np.sqrt(np.mean(fit.fun ** 2))
    a_sim, _ = d_free(ts + 1.0, 2 * xs)
    a_ode, _ = d_free(ts + 1.0, 2 * integrate(xs[0], ts, 0.1, s_sim, c_cal))
    print(f"[1] Fig.2 run: c = {c_cal:.3f}, s = {s_sim:.3e} "
          f"(wall 1-h = 2.85e-3), rms = {rms*100:.2f} %, a sim {a_sim:.3f} / ODE {a_ode:.3f}")
    rows.append(dict(step="1 calibrate on Fig.2 run", alpha_c=0.1, c=c_cal, s=s_sim,
                     rms_pct=100 * rms, a_target=a_sim, a_ode=a_ode, note="alpha fixed"))

    # ---- 2. validate on mesh_pair fine (alpha 1e-3, sealed, tangent) ------
    t_grid = np.geomspace(1.0, 3.0e5, 400)
    xg = integrate(0.3e-6, t_grid, 1e-3, 0.0, c_cal)
    t_anchor = np.interp(16.405e-6, xg, t_grid)
    t_end = np.interp(25.55e-6, xg, t_grid)
    mpost = (xg >= 16.405e-6) & (xg <= 25.55e-6)
    a_mp, _ = d_free(t_grid[mpost] - t_anchor + 1.0, 2 * xg[mpost])
    print(f"[2] mesh_pair fine: anchor at {t_anchor:.0f} s (run 54174 s), "
          f"51.1 um at {t_end/3600:.1f} h (run 78.6 h), a = {a_mp:.3f} (run 0.283)")
    rows.append(dict(step="2 validate on mesh_pair fine", alpha_c=1e-3, c=c_cal, s=0.0,
                     rms_pct=float("nan"), a_target=0.283, a_ode=a_mp,
                     note=f"t_anchor {t_anchor:.0f}s vs 54174s; t_end {t_end/3600:.1f}h vs 78.6h"))

    # ---- 3. fit the Molaro data -------------------------------------------
    td, xd, Rl = load_data()
    a_data, a_data_ci = d_free(td + 1.0, 2 * xd)

    def fit_data(alpha=None):
        if alpha is None:
            f = lambda p: (integrate(xd[0], td, 10 ** p[0], p[1], c_cal) - xd) / xd
            r = least_squares(f, [-1.0, 1e-4], bounds=([-4, 0], [2, 5e-2]))
            return 10 ** r.x[0], r.x[1], r
        f = lambda p: (integrate(xd[0], td, alpha, p[0], c_cal) - xd) / xd
        r = least_squares(f, [1e-5], bounds=([0], [5e-2]))
        return alpha, r.x[0], r

    alpha_b, s_b, r_b = fit_data()
    xb = integrate(xd[0], td, alpha_b, s_b, c_cal)
    a_b, _ = d_free(td + 1.0, 2 * xb)
    print(f"[3] Molaro data free fit: alpha_c = {alpha_b:.3g}, s = {s_b:.3e}, "
          f"rms {100*np.sqrt(np.mean(r_b.fun**2)):.2f} %, a ODE {a_b:.3f} vs data "
          f"{a_data:.3f} +- {a_data_ci:.3f}")
    rows.append(dict(step="3 free fit to Molaro data", alpha_c=alpha_b, c=c_cal, s=s_b,
                     rms_pct=100 * np.sqrt(np.mean(r_b.fun ** 2)), a_target=a_data,
                     a_ode=a_b, note="alpha and s free"))

    prof = []
    for alpha in [1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0]:
        _, s_a, r_a = fit_data(alpha)
        x_a = integrate(xd[0], td, alpha, s_a, c_cal)
        a_a = d_free(td + 1.0, 2 * x_a)[0] if np.all(np.isfinite(x_a)) else float("nan")
        rms_a = 100 * np.sqrt(np.mean(r_a.fun ** 2))
        prof.append((alpha, s_a, rms_a, a_a, x_a))
        print(f"    alpha {alpha:7.0e}: best s {s_a:.2e}, rms {rms_a:6.2f} %, a {a_a:.3f}")
        rows.append(dict(step="3 profile", alpha_c=alpha, c=c_cal, s=s_a, rms_pct=rms_a,
                         a_target=a_data, a_ode=a_a, note="s fitted, alpha fixed"))

    # ---- 4. shrinkage -> far-field undersaturation -------------------------
    sl = stats.linregress(td, Rl)
    for alpha in (1e-3, 0.1):
        s_inf = -sl.slope * (beta_sub(alpha) + Rl[0] * TH["k_D"]) - 2 * TH["d0"] / Rl[0]
        print(f"[4] large grain dR/dt = {sl.slope:.3e} m/s -> s_inf = {s_inf:.2e} "
              f"at alpha {alpha:g}")
        rows.append(dict(step="4 shrinkage -> s_inf", alpha_c=alpha, c=float("nan"),
                         s=s_inf, rms_pct=float("nan"), a_target=float("nan"),
                         a_ode=float("nan"), note=f"dR/dt {sl.slope:.3e} m/s"))

    with open(RES / "envelope_fits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # ---- figure ------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
    ax[0].plot(ts / 60, 2e6 * xs, color="#5c5c5c", lw=2.5, label="Fig. 2 run (α_c = 0.1)")
    ax[0].plot(ts / 60, 2e6 * integrate(xs[0], ts, 0.1, s_sim, c_cal), "--",
               color="#0072B2", lw=1.5, label=f"ODE, c = {c_cal:.2f}, s = {s_sim:.1e}")
    ax[0].set_title("(a) calibration on our −20 °C run", loc="left", fontsize=10)
    tt = np.linspace(0, td[-1], 200)
    ax[1].errorbar(td / 60, 2e6 * xd, fmt="o", color="#1a1a1a", ms=5,
                   label=f"Molaro data, a = {a_data:.2f}")
    cols = {1e-3: "#CC79A7", 1e-2: "#E69F00", 0.1: "#0072B2", 1.0: "#009E73"}
    for alpha, s_a, rms_a, a_a, _ in prof:
        if alpha in cols:
            ax[1].plot(tt / 60, 2e6 * integrate(xd[0], tt, alpha, s_a, c_cal),
                       color=cols[alpha], lw=1.5,
                       label=(f"α_c = {alpha:g}: rms {rms_a:.0f} %, a = {a_a:.2f}"
                              if alpha >= 0.1 else  # too little growth to fit an a
                              f"α_c = {alpha:g}: rms {rms_a:.0f} %"))
    ax[1].set_title("(b) best fit at each α_c (s fitted)", loc="left", fontsize=10)
    for a_ in ax:
        a_.set_xlabel("time [min]")
        a_.set_ylabel("neck width, 2x [µm]")
        a_.legend(fontsize=8, frameon=False)
        a_.grid(color="#d8d8d8", lw=0.5)
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    for ext in ("png", "pdf"):
        fig.savefig(RES / f"envelope_fit.{ext}", dpi=200)
    print(f"wrote {RES}/envelope_fits.csv, envelope_fit.png/.pdf")


if __name__ == "__main__":
    main()
