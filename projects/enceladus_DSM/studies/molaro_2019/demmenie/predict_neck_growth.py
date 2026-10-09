#!/usr/bin/env python3
"""predict_neck_growth.py — neck growth predicted from the Gibbs–Thomson condition.

    venv_enceladus/bin/python studies/molaro_2019/demmenie/predict_neck_growth.py <run dir>

The interface condition the solver obeys (lunar docs/curvature_driven_growth.md;
projection of the phase-field residual on the translation mode) is

    sigma = d0 kappa + beta v_n ,

so a surface in vapour of supersaturation sigma moves at
v_n = (sigma - d0 kappa) / beta. For a neck of radius x between grains of
radius R, with the Kuczynski fillet rho = x^2 / (2 (R - x)), fed from the grain
surface through attachment in series with gas diffusion over a length c rho,

    v_n   = d0 (1/rho - 1/x + 2/R) / (beta_sub + c rho k_D) ,   dx/dt = (pi/2) v_n

(studies/molaro_2019/envelope_ode/: the same one-equation model, c = 0.275
calibrated there on the alpha_c = 0.1 Molaro run and NOT refitted here).
beta_sub ~ 1/alpha_c is attachment, k_D = rho_ice / (D_v rho_vs) gas diffusion.
Alone, attachment gives x^3 ~ t (a = 1/3) and diffusion x^5 ~ t (a = 1/5).

What this script does, for the saturated equal-grain run (-20 C, R = 84.4 um,
alpha_c = 1e-3, sigma_far = 2 d0 / R):
  (a) the predicted neck velocity v_n against neck width -- attachment only,
      diffusion only, both -- with the run's own (2/pi) d(w/2)/dt
  (b) neck width against time: run and predictions, started from the run's
      first sample after the relaxation period
  (c) the exponent of C (t + t0)^a against the start of the fit window: run,
      prediction, and the range of Demmenie et al. (2025)
  (d) the prediction at Demmenie's OWN conditions (-3 C, tangent start, 2.5 h)
      against alpha_c, for two grain radii: what the vapor route alone gives
      there

THE THIN-INTERFACE TERM (lunar docs/enceladus_carryover.md, docs/gt_deficit/).
Until commit 5ecdd236 (2026-09-09) tau_sub carried Karma counter-terms that the
one-sided diffusivity never subtracts back, so the realised beta was LARGER
than requested by an additive amount ~ eps: the interface was slower, and the
relative error grows with alpha_c. The script reads each run's tau_sub and
reports beta_realised / beta_requested = tau_sub d0 / (eps^2 beta_sub0). The
saturated run postdates the fix (ratio 1.000).

Writes prediction/neck_prediction.png and prediction/neck_prediction.txt.
"""
from __future__ import annotations

import math, re, sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "preprocess")); sys.path.insert(0, str(HERE))
import comp_eps as ce  # noqa: E402
from analyze_demmenie import fit_free, window_scan, load, read_tau_sub, DEMMENIE, HOUR  # noqa: E402

RHO_I, C_GEO, UM = 919.0, 0.275, 1e-6


def thermo(T_C):
    rvs = ce.rho_vs_sat(T_C)
    return dict(d0=ce.capillary_length(T_C), k_D=RHO_I / (ce.Dv_T(T_C) * rvs),
                b1=ce.beta_HK(T_C, 1.0) * RHO_I / rvs)        # beta_sub at alpha_c = 1


def v_neck(x, R, th, alpha, c=C_GEO, att=True, dif=True):
    """Normal velocity of the neck surface [m/s] from Gibbs–Thomson."""
    rho = x * x / (2.0 * (R - x))
    drive = th["d0"] * (1.0 / rho - 1.0 / x + 2.0 / R)
    return drive / ((th["b1"] / alpha if att else 0.0) + (c * rho * th["k_D"] if dif else 0.0))


def grow(x0, t, R, th, alpha, **kw):
    f = lambda tt, y: [0.5 * math.pi * v_neck(y[0], R, th, alpha, **kw)]
    s = solve_ivp(f, (t[0], t[-1]), [x0], t_eval=t, method="LSODA", rtol=1e-9, atol=1e-14)
    return s.y[0]


def beta_ratio(run):
    txt = (run / "outp.txt").read_text(errors="replace")
    g = lambda k: float(re.search(rf"^\s*{k}\s+([0-9.eE+-]+)", txt, re.M).group(1))
    return g("tau_sub") * g("d0_sub0") / (g("eps") ** 2 * g("beta_sub"))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    run = Path(sys.argv[1])
    out = HERE / "prediction"; out.mkdir(exist_ok=True)
    L = []
    P = lambda s="": (print(s), L.append(s))

    T_C, R, alpha = -20.0, 84.4e-6, 1e-3
    th = thermo(T_C)
    t, w, _ = load(run); tau = read_tau_sub(run)
    keep = t >= 11 * tau
    tk, xk = t[keep], 0.5 * w[keep] * UM
    P(f"run: {run.name[:70]}")
    P(f"  beta_realised / beta_requested = {beta_ratio(run):.4f}  (1.0000 = no thin-interface inflation)")
    P(f"  d0 = {th['d0']:.4e} m, beta_sub = {th['b1'] / alpha:.4e} s/m, k_D = {th['k_D']:.4e} s/m2, c = {C_GEO}")
    xr = np.array([xk[0], xk[-1]]); rr = xr ** 2 / (2 * (R - xr))
    P(f"  diffusion / attachment resistance, c rho k_D / beta_sub: {C_GEO * rr[0] * th['k_D'] / (th['b1'] / alpha):.3f} at w = {2 * xr[0] / UM:.1f} um,"
      f" {C_GEO * rr[1] * th['k_D'] / (th['b1'] / alpha):.3f} at w = {2 * xr[1] / UM:.1f} um")

    preds = {"attachment + diffusion": dict(), "attachment only": dict(dif=False), "diffusion only": dict(att=False)}
    X = {k: grow(xk[0], tk, R, th, alpha, **kw) for k, kw in preds.items()}
    col = {"attachment + diffusion": "#2a78d6", "attachment only": "#d1495b", "diffusion only": "#3a9a5c"}
    P()
    P("  neck width at the end of the run (100 h), from the first sample after relaxation:")
    P(f"    run                       {2 * xk[-1] / UM:6.2f} um")
    for k in preds:
        P(f"    {k:25s} {2 * X[k][-1] / UM:6.2f} um")
    P()
    P("  exponent a of C (t + t0)^a over the samples after relaxation:")
    fr = fit_free(tk, 2 * xk / UM); P(f"    run                       {fr['a']:.3f}")
    A = {}
    for k in preds:
        A[k] = fit_free(tk, 2 * X[k] / UM); P(f"    {k:25s} {A[k]['a']:.3f}")

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
    # (a) velocity at the neck
    xs = np.linspace(0.8 * xk[0], 1.05 * xk[-1], 300)
    for k, kw in preds.items():
        ax[0, 0].semilogy(2 * xs / UM, v_neck(xs, R, th, alpha, **kw), color=col[k], lw=1.8, label=k)
    vr = (2 / math.pi) * np.gradient(xk, tk)
    ax[0, 0].semilogy(2 * xk[1:-1] / UM, vr[1:-1], "o", ms=3.5, mfc="white", mec="k", mew=0.8,
                      label=r"run, $(2/\pi)\,\mathrm{d}(w/2)/\mathrm{d}t$")
    ax[0, 0].set(xlabel=r"neck width $w$ [$\mu$m]", ylabel=r"neck velocity $v_n$ [m/s]")
    ax[0, 0].legend(frameon=False, fontsize=9)
    # (b) w(t)
    ax[0, 1].plot(tk / HOUR, 2 * xk / UM, "k-", lw=2.2, label="run")
    for k in preds:
        ax[0, 1].plot(tk / HOUR, 2 * X[k] / UM, color=col[k], lw=1.5, ls="--", label=f"{k} (a = {A[k]['a']:.2f})")
    ax[0, 1].set(xlabel="time [h]", ylabel=r"neck width $w$ [$\mu$m]", ylim=(30, 110))
    ax[0, 1].legend(frameon=False, fontsize=9, title=f"run a = {fr['a']:.2f}", title_fontsize=9)
    # (c) exponent against window start
    s = window_scan(tk, 2 * xk / UM)
    ax[1, 0].axhspan(*DEMMENIE, color="#d1495b", alpha=0.13, lw=0)
    ax[1, 0].axhline(1 / 3, color="0.4", lw=0.8, ls=":"); ax[1, 0].axhline(1 / 5, color="0.4", lw=0.8, ls=":")
    ax[1, 0].plot(s[:, 0] / HOUR, s[:, 2], "k-", lw=2.2, label="run")
    for k in preds:
        sp = window_scan(tk, 2 * X[k] / UM)
        ax[1, 0].plot(sp[:, 0] / HOUR, sp[:, 2], color=col[k], lw=1.5, ls="--", label=k)
    ax[1, 0].text(1, 1 / 3 + 0.004, "1/3", fontsize=8); ax[1, 0].text(1, 1 / 5 + 0.004, "1/5", fontsize=8)
    ax[1, 0].text(40, np.mean(DEMMENIE), "Demmenie et al. (2025)", fontsize=8, va="center")
    ax[1, 0].set(xlabel="start of fit window [h]", ylabel="exponent $a$", ylim=(0.15, 0.37))
    ax[1, 0].legend(frameon=False, fontsize=9, loc="lower right")
    # (d) Demmenie's own conditions
    P()
    P("  vapor route at Demmenie's conditions (-3 C, tangent start, 2.5 h; c as above):")
    thd = thermo(-3.0)
    td = np.geomspace(60.0, 2.5 * 3600, 60)
    al = np.geomspace(1e-4, 1.0, 17)
    for Rd, ls in ((0.5e-3, "-"), (1.0e-3, "--")):
        aa, ratio = [], []
        for a_c in al:
            xd = grow(1e-6, np.concatenate([[0.0], td]), Rd, thd, a_c)[1:]
            aa.append(fit_free(td, 2 * xd / UM)["a"])
            rho_end = xd[-1] ** 2 / (2 * (Rd - xd[-1]))
            ratio.append(C_GEO * rho_end * thd["k_D"] / (thd["b1"] / a_c))
        ax[1, 1].semilogx(al, aa, "k" + ls, lw=1.8, label=f"R = {Rd * 1e3:g} mm")
        for a_c in (1e-3, 1e-2, 1e-1, 1.0):
            i = int(np.argmin(np.abs(al - a_c)))
            P(f"    R = {Rd * 1e3:g} mm, alpha_c = {a_c:g}: a = {aa[i]:.3f}; diffusion/attachment resistance at 2.5 h = {ratio[i]:.3g}")
    ax[1, 1].axhspan(*DEMMENIE, color="#d1495b", alpha=0.13, lw=0)
    ax[1, 1].axhline(1 / 3, color="0.4", lw=0.8, ls=":"); ax[1, 1].axhline(1 / 5, color="0.4", lw=0.8, ls=":")
    ax[1, 1].text(1.2e-4, np.mean(DEMMENIE), "Demmenie et al. (2025), measured", fontsize=8, va="center")
    ax[1, 1].set(xlabel=r"condensation coefficient $\alpha_c$", ylabel="predicted exponent $a$", ylim=(0.15, 0.37))
    ax[1, 1].legend(frameon=False, fontsize=9, loc="lower left")
    for a_, lab, ttl in zip(ax.ravel(), "abcd", ("neck velocity from Gibbs–Thomson, our run", "neck width, our run",
                                                 "fitted exponent, our run", "vapor route at Demmenie's conditions")):
        a_.set_title(f"({lab}) {ttl}", fontsize=10, loc="left")
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(out / "neck_prediction.png", dpi=180); plt.close(fig)
    (out / "neck_prediction.txt").write_text("\n".join(L) + "\n")
    P(f"\nwrote {out}/neck_prediction.png and .txt")


if __name__ == "__main__":
    main()
