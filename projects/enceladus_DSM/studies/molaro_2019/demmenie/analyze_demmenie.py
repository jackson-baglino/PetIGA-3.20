#!/usr/bin/env python3
"""analyze_demmenie.py — neck growth of the saturated equal-grain run against t^(1/3).

    venv_enceladus/bin/python studies/molaro_2019/demmenie/analyze_demmenie.py <run dir>
        [--relax-tau 11] [--out <dir>]

The run is the Demmenie-conditions test (studies/molaro_2019/TODO.md): two
equal grains in vapour saturated over their own curvature, so nothing drives
mass transfer but the neck. Demmenie, Woutersen & Bonn (2025) measure
r ~ t^a with a = 0.26-0.33 there, the evaporation-condensation law a = 1/3.
This run is NOT compared with the Molaro et al. (2019) cryostage data: that
experiment was unsaturated and its grains shrank, a different problem.

What is fitted. The run starts from a 14 um neck, not from tangent contact, so
it has no physical t = 0 and a bare C t^a fit is meaningless. Demmenie's own
form w = C (t + t0)^a is used, t0 free, over the samples AFTER the relaxation
period: the first --relax-tau phase-change times tau_sub (default 11, the
k_eff campaign's convention), during which the sharp initial cusp relaxes to
the diffuse profile. Three fits are reported on that window:

    free      w = C (t + t0)^a          C, t0, a free
    one-third w = C (t + t0)^(1/3)      C, t0 free -- the law being tested
    Kuczynski w^m - w_r^m = K (t - t_r) m, K free; w_r, t_r = first sample kept

and the free exponent again for fit windows that start later and later, which
is the honest way to show a fitted exponent: it is a property of the curve AND
the window.

Writes into <run>/plots/demmenie/ (or --out):
    neck_onethird.png    w(t) with the free and the one-third fit
    neck_loglog.png      w against t + t0 (free fit's t0), log-log, slope guides
    exponent_window.png  fitted a against the start of the fit window
    grain_radius.png     grain radius against time: the saturation check
    fits.txt             the numbers
Reads neck_width.csv and grain_shrinkage.csv (postprocess/run_batch_measure.sh).
"""
from __future__ import annotations

import argparse, re, sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

UM, HOUR = 1e-6, 3600.0
DEMMENIE = (0.26, 0.33)          # the four measured exponents span this


def read_tau_sub(run: Path) -> float:
    m = re.search(r"^\s*tau_sub\s+([0-9.eE+-]+)\s*s", (run / "outp.txt").read_text(errors="replace"), re.M)
    if not m:
        sys.exit("tau_sub not found in outp.txt")
    return float(m.group(1))


def load(run: Path):
    n = np.genfromtxt(run / "neck_width.csv", delimiter=",", names=True)
    t, w = n["t_s"], n["neck_width_m"] / UM
    g = None
    if (run / "grain_shrinkage.csv").is_file():
        gg = np.genfromtxt(run / "grain_shrinkage.csv", delimiter=",", names=True)
        g = (gg["t_s"], gg["R_large_m"] / UM)
    return t, w, g


def _pl(t, C, t0, a):
    return C * (t + t0) ** a


def _p3(t, C, t0):
    return C * (t + t0) ** (1.0 / 3.0)


def _rms(model, w):
    r = model / w - 1.0
    return 100 * float(np.sqrt(np.mean(r ** 2))), 100 * float(np.max(np.abs(r)))


def fit_free(t, w):
    p, c = curve_fit(_pl, t, w, p0=(5.0, 1e4, 0.2), maxfev=100000,
                     bounds=([0, -t.min() + 1e-6, 0.01], [np.inf, np.inf, 1.0]))
    rms, mx = _rms(_pl(t, *p), w)
    return dict(C=p[0], t0=p[1], a=p[2], a_se=float(np.sqrt(c[2, 2])), rms=rms, max=mx)


def fit_third(t, w):
    p, c = curve_fit(_p3, t, w, p0=(1.0, 1e4), maxfev=100000,
                     bounds=([0, -t.min() + 1e-6], [np.inf, np.inf]))
    rms, mx = _rms(_p3(t, *p), w)
    return dict(C=p[0], t0=p[1], a=1 / 3, rms=rms, max=mx)


def fit_kuczynski(t, w):
    tr, wr = t[0], w[0]
    f = lambda tt, K, m: (wr ** m + K * (tt - tr)) ** (1.0 / m)
    p, c = curve_fit(f, t[1:], w[1:], p0=(1e3, 4.0), maxfev=100000)
    rms, mx = _rms(f(t[1:], *p), w[1:])
    return dict(K=p[0], m=p[1], m_se=float(np.sqrt(c[1, 1])), rms=rms, max=mx)


def window_scan(t, w, min_pts=25):
    """Free-fit exponent for windows [t_i, end], t_i walking forward."""
    out = []
    for i in range(0, len(t) - min_pts):
        try:
            f = fit_free(t[i:], w[i:])
        except Exception:
            continue
        out.append((t[i], w[i], f["a"], f["a_se"]))
    return np.array(out)


def analyse(run: Path, relax_tau: float):
    t, w, g = load(run)
    tau = read_tau_sub(run)
    keep = t >= relax_tau * tau
    tk, wk = t[keep], w[keep]
    return dict(t=t, w=w, g=g, tau=tau, t_relax=relax_tau * tau, tk=tk, wk=wk,
                free=fit_free(tk, wk), third=fit_third(tk, wk), kucz=fit_kuczynski(tk, wk),
                scan=window_scan(tk, wk))


def report(A, relax_tau):
    f, h, k, s = A["free"], A["third"], A["kucz"], A["scan"]
    L = [f"tau_sub = {A['tau']:.1f} s; relaxation period = {relax_tau:g} tau_sub = {A['t_relax'] / 60:.0f} min",
         f"samples: {len(A['t'])} in all, {len(A['tk'])} after relaxation "
         f"(t = {A['tk'][0] / HOUR:.2f}-{A['tk'][-1] / HOUR:.1f} h, w = {A['wk'][0]:.1f}-{A['wk'][-1]:.1f} um)",
         "",
         "fits over the samples after relaxation:",
         f"  free       w = C (t + t0)^a      a = {f['a']:.3f} +- {f['a_se']:.3f}   t0 = {f['t0'] / HOUR:.2f} h   "
         f"rms {f['rms']:.2f} %  max {f['max']:.2f} %",
         f"  one-third  w = C (t + t0)^(1/3)  a = 0.333 (fixed)     t0 = {h['t0'] / HOUR:.2f} h   "
         f"rms {h['rms']:.2f} %  max {h['max']:.2f} %",
         f"  Kuczynski  w^m - w_r^m = K dt    m = {k['m']:.2f} +- {k['m_se']:.2f}  (1/m = {1 / k['m']:.3f})   "
         f"rms {k['rms']:.2f} %  max {k['max']:.2f} %",
         "",
         "free exponent against the start of the fit window:"]
    for frac in (0.0, 0.1, 0.25, 0.5, 0.7):
        i = int(np.argmin(np.abs(s[:, 0] - (A["tk"][0] + frac * (s[-1, 0] - A["tk"][0])))))
        L.append(f"  from t = {s[i, 0] / HOUR:5.1f} h (w = {s[i, 1]:.1f} um): a = {s[i, 2]:.3f} +- {s[i, 3]:.3f}")
    L += ["", f"Demmenie et al. (2025): a = {DEMMENIE[0]}-{DEMMENIE[1]} (four runs); evaporation-condensation law a = 1/3."]
    if A["g"] is not None:
        gt, gR = A["g"]
        L.append(f"grain radius: {gR[0]:.3f} -> {gR[-1]:.3f} um ({100 * (gR[-1] / gR[0] - 1):+.2f} %) over the run")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--relax-tau", type=float, default=11.0)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    A = analyse(a.run, a.relax_tau)
    out = a.out or a.run / "plots" / "demmenie"
    out.mkdir(parents=True, exist_ok=True)
    txt = report(A, a.relax_tau)
    (out / "fits.txt").write_text(txt + "\n")
    print(txt)

    t, w, tk, wk, f, h = A["t"], A["w"], A["tk"], A["wk"], A["free"], A["third"]
    tt = np.linspace(tk[0], tk[-1], 400)
    c_m, c_f, c_3 = "#1a1a1a", "#2a78d6", "#d1495b"

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axvspan(0, A["t_relax"] / HOUR, color="0.9", lw=0)
    ax.plot(t / HOUR, w, "o", ms=3.5, mfc="white", mec=c_m, mew=0.8, label="simulation")
    ax.plot(tt / HOUR, _pl(tt, f["C"], f["t0"], f["a"]), color=c_f, lw=1.6,
            label=rf"$C(t+t_0)^a$, $a={f['a']:.3f}$  (rms {f['rms']:.2f} %)")
    ax.plot(tt / HOUR, _p3(tt, h["C"], h["t0"]), color=c_3, lw=1.6, ls="--",
            label=rf"$C(t+t_0)^{{1/3}}$  (rms {h['rms']:.2f} %)")
    ax.set(xlabel="time [h]", ylabel=r"neck width $w$ [$\mu$m]")
    ax.legend(frameon=False, fontsize=9); fig.tight_layout(); fig.savefig(out / "neck_onethird.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    x = tk + f["t0"]
    ax.loglog(x / HOUR, wk, "o", ms=3.5, mfc="white", mec=c_m, mew=0.8, label="simulation (after relaxation)")
    ax.loglog(x / HOUR, _pl(tk, f["C"], f["t0"], f["a"]), color=c_f, lw=1.4, label=rf"slope $a={f['a']:.3f}$")
    ax.loglog(x / HOUR, wk[-1] * (x / x[-1]) ** (1 / 3), color=c_3, lw=1.4, ls="--", label="slope 1/3 through the last sample")
    ax.set(xlabel=rf"$t+t_0$ [h]   ($t_0$ = {f['t0'] / HOUR:.2f} h, free fit)", ylabel=r"neck width $w$ [$\mu$m]")
    ax.legend(frameon=False, fontsize=9); fig.tight_layout(); fig.savefig(out / "neck_loglog.png", dpi=200); plt.close(fig)

    s = A["scan"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.axhspan(*DEMMENIE, color="#d1495b", alpha=0.15, lw=0, label="Demmenie et al. (2025), four runs")
    ax.axhline(1 / 3, color=c_3, lw=1.0, ls="--", label="1/3")
    ax.plot(s[:, 0] / HOUR, s[:, 2], color=c_f, lw=1.6, label="free-fit exponent, window from here to 100 h")
    ax.fill_between(s[:, 0] / HOUR, s[:, 2] - 2 * s[:, 3], s[:, 2] + 2 * s[:, 3], color=c_f, alpha=0.2, lw=0)
    ax.set(xlabel="start of the fit window [h]", ylabel=r"exponent $a$", ylim=(0.15, 0.37))
    ax.legend(frameon=False, fontsize=9); fig.tight_layout(); fig.savefig(out / "exponent_window.png", dpi=200); plt.close(fig)

    if A["g"] is not None:
        gt, gR = A["g"]
        fig, ax = plt.subplots(figsize=(6.4, 3.6))
        ax.plot(gt / HOUR, 100 * (gR / gR[0] - 1), color=c_m, lw=1.4)
        ax.set(xlabel="time [h]", ylabel=r"change of grain radius [%]")
        fig.tight_layout(); fig.savefig(out / "grain_radius.png", dpi=200); plt.close(fig)
    print(f"\nwrote {out}/neck_onethird.png, neck_loglog.png, exponent_window.png, grain_radius.png, fits.txt")


if __name__ == "__main__":
    main()
