#!/usr/bin/env python3
"""fit_molaro_exponents.py — neck-growth exponents of the two Molaro runs.

    venv_enceladus/bin/python studies/molaro_2019/exponents/fit_molaro_exponents.py

The same fit as the saturated (Demmenie-conditions) run, so the three
exponents can be set side by side: w = C (t + t0)^a with t0 free, over the
samples AFTER the relaxation period (the first 11 tau_sub, the convention of
studies/molaro_2019/demmenie/analyze_demmenie.py, whose fit routines are used
here unchanged).

For each of the two final manuscript runs (alpha_c = 0.1, nominal D_v,
thin-interface term on) it reports the exponent over
  (i)   every sample after relaxation (the whole run),
  (ii)  the part of the run that overlaps the experiment: from the time the
        model neck equals Molaro's first measurement (their time anchor) for
        the length of their series (78 min at -20 C, 48 min at -5 C),
and the same fit to Molaro's own measurements, plus how the model exponent
moves with the start of the fit window.

Writes exponents.txt and exponents.png next to this file.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "demmenie"))
import analyze_demmenie as ad                                   # noqa: E402

G = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/GrainPairSintering"
RUNS = {
    "-20 C": (G / "batch_2026-09-08__17.20.46_molaro_T-20_round2"
                / "molaro_2D_L450x225um_eps0.12um_axisym_T-20pair_r14um_dom2__molaro_T-20_h0.99715_2h_a1e-1_dirichlet",
              REPO / "inputs/validation/molaro2019_fig11_T-20.csv"),
    "-5 C": (G / "batch_2026-10-09__08.51.38_molaro_T-5_corr_on"
               / "molaro_2D_L450x225um_eps0.11um_axisym_T-5pair_r14um_dom2__molaro_T-5_h0.99674_2h_a1e-1_dirichlet",
             REPO / "inputs/validation/molaro2019_fig11_T-5.csv"),
}
RELAX_TAU = 11.0


def main():
    L = []
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    for k, (name, (run, data)) in enumerate(RUNS.items()):
        t, w, _ = ad.load(run)
        tau = ad.read_tau_sub(run)
        keep = t >= RELAX_TAU * tau
        tk, wk = t[keep], w[keep]
        d = np.genfromtxt(data, delimiter=",", comments="#")
        td, wd = d[:, 0] * 60.0, d[:, 1]
        t_star = float(np.interp(wd[0], wk, tk))                 # model time at Molaro's first width
        ov = (tk >= t_star) & (tk <= t_star + td[-1])
        f_all, f_ov, f_dat = ad.fit_free(tk, wk), ad.fit_free(tk[ov], wk[ov]), ad.fit_free(td + 1.0, wd)
        scan = ad.window_scan(tk, wk)
        L += [f"{name}:  tau_sub = {tau:.3f} s; relaxation period = {RELAX_TAU:g} tau_sub = {RELAX_TAU * tau:.1f} s",
              f"  samples: {len(t)} in all, {len(tk)} after relaxation, {int(ov.sum())} in the experiment's window",
              f"  model reaches Molaro's first width ({wd[0]:.2f} um) at t* = {t_star / 60:.1f} min; run ends at {t[-1] / 60:.0f} min",
              f"  model, whole run after relaxation ({tk[0] / 60:.2f}-{tk[-1] / 60:.0f} min):  a = {f_all['a']:.3f} +- {f_all['a_se']:.3f}  (t0 = {f_all['t0'] / 60:.1f} min, rms {f_all['rms']:.2f} %)",
              f"  model, experiment's window (t* to t* + {td[-1] / 60:.0f} min):      a = {f_ov['a']:.3f} +- {f_ov['a_se']:.3f}  (t0 = {f_ov['t0'] / 60:.1f} min, rms {f_ov['rms']:.2f} %)",
              f"  Molaro's measurements ({len(td)} points, {td[-1] / 60:.0f} min):          a = {f_dat['a']:.3f} +- {f_dat['a_se']:.3f}  (t0 = {f_dat['t0'] / 60:.1f} min, rms {f_dat['rms']:.2f} %)",
              "  model exponent against the start of the fit window (window runs to the end of the run):"]
        for frac in (0.0, 0.1, 0.25, 0.5):
            i = min(int(frac * (len(scan) - 1)), len(scan) - 1)
            L.append(f"    from {scan[i, 0] / 60:6.1f} min (w = {scan[i, 1]:.1f} um):  a = {scan[i, 2]:.3f}")
        L.append("")
        ax[k].plot(scan[:, 0] / 60, scan[:, 2], color="k", lw=1.6, label="model, window from here to end of run")
        ax[k].axhline(f_dat["a"], color="C3", lw=1.0, ls="--", label=f"Molaro data, a = {f_dat['a']:.2f} $\\pm$ {f_dat['a_se']:.2f}")
        ax[k].axvline(t_star / 60, color="0.6", lw=0.8, ls=":")
        for y, s_ in ((1 / 3, "1/3"), (1 / 5, "1/5"), (1 / 7, "1/7")):
            ax[k].axhline(y, color="0.75", lw=0.6); ax[k].text(ax[k].get_xlim()[1], y, " " + s_, va="center", fontsize=8, color="0.4")
        ax[k].set(xlabel="start of the fit window [min]", ylabel="exponent $a$", title=name, ylim=(0.05, 0.40))
        ax[k].legend(frameon=False, fontsize=8)
    txt = "\n".join(L)
    print(txt)
    (HERE / "exponents.txt").write_text(txt + "\n")
    fig.tight_layout(); fig.savefig(HERE / "exponents.png", dpi=200)
    print(f"wrote {HERE}/exponents.txt, exponents.png")


if __name__ == "__main__":
    main()
