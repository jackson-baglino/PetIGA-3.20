"""Plan issue 1: is the 0.86 um mesh accurate enough for the fast narrow-plug
phase? Compares the eps = 0.328 um rerun (I1-fine, 1 day) with 0c-C1
(eps = 0.858 um, same walls, plugs and boundary conditions) over day 1.

Measured, all from the solver's own files and ice_bodies.csv:

  A(t)        total ice area, SSA_evo.dat (every step); relative difference
              (fine - coarse) / coarse on a common time grid
  lost(1 d)   ice lost by day 1 as a fraction of the start, and the relative
              difference of the two
  plug(t)     area of the narrow plug: every body whose com_x is below
              --split (default 450 um, between the plugs) in ice_bodies.csv
  t_gone      narrow-plug vanishing time: the first snapshot with no body
              left below --split (bracketed by the previous snapshot), and the
              time of the last step with Pi > 0.5 after it began (resolution.dat)
  v_front     median v_n over 2e3-2e4 s (resolution.dat), the steady
              narrow-plug front, and the Pi median there

The acceptance rule (RUNS.md I1-fine): under 5 % difference in the vanishing
time and the area curve keeps the 0.86 um mesh for production.

    venv_cryosuction/bin/python studies/cryosuction/scripts/compare_issue1.py \
        --fine <I1 run> --coarse <0c-C1 run> --out <folder>
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "postprocess"))
import pplib  # noqa: E402

DAY = 86400.0


def area(run):
    s = pplib.load_ssa(run)
    o = np.argsort(s[:, pplib.TIME])
    return s[o, pplib.TIME], s[o, pplib.ICE]


def plug(run, split):
    import pandas as pd
    b = pd.read_csv(os.path.join(run, "ice_bodies.csv"))
    g = b[b.com_x < split].groupby("t_s").area_m2.sum()
    t = np.sort(b.t_s.unique())
    return t, g.reindex(t, fill_value=0.0).to_numpy()


def resolution(run):
    a = np.genfromtxt(os.path.join(run, "resolution.dat"), comments="#", invalid_raise=False)
    return a[:, 0], a[:, 3], a[:, 4]


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fine", required=True)
    p.add_argument("--coarse", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--split", type=float, default=4.5e-4, help="x between the plugs [m]")
    a = p.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    runs = {"fine (eps 0.328 µm)": a.fine, "coarse 0c-C1 (eps 0.858 µm)": a.coarse}
    print("fine:", a.fine, "\ncoarse:", a.coarse, f"\nsplit x = {a.split*1e6:.0f} um")

    tf, Af = area(a.fine)
    tc, Ac = area(a.coarse)
    tg = np.linspace(0, DAY, 1441)
    af, ac = np.interp(tg, tf, Af), np.interp(tg, tc, Ac)
    rel = (af - ac) / ac
    lost_f, lost_c = (af[0] - af[-1]) / af[0], (ac[0] - ac[-1]) / ac[0]
    print(f"A0 fine {af[0]:.6e} coarse {ac[0]:.6e} m^2 (rel {(af[0]-ac[0])/ac[0]:+.2e})")
    print(f"max |A_fine - A_coarse| / A_coarse over day 1: {np.abs(rel).max():.4f} at t = {tg[np.abs(rel).argmax()]:.0f} s")
    print(f"ice lost by day 1: fine {lost_f:.4f}, coarse {lost_c:.4f}, relative difference {(lost_f-lost_c)/lost_c:+.4f}")

    res = {}
    for k, r in runs.items():
        tp, ap = plug(r, a.split)
        gone = np.flatnonzero((ap == 0) & (tp > 0))
        i = gone[0] if gone.size else None
        tr, vn, Pi = resolution(r)
        w = (tr > 2e3) & (tr < 2e4)
        late = tr[(Pi > 0.5) & (tr > 2e4) & (tr < DAY)]
        res[k] = dict(tp=tp, ap=ap, tr=tr, Pi=Pi,
                      gone=(tp[i - 1], tp[i]) if i else (np.nan, np.nan),
                      last_pi=late.max() if late.size else np.nan,
                      vn=np.median(vn[w]), Pi_med=np.median(Pi[w]), Pi_max=Pi[w].max())
        print(f"{k}: narrow plug gone between {res[k]['gone'][0]:.0f} and {res[k]['gone'][1]:.0f} s; "
              f"last Pi > 0.5 at {res[k]['last_pi']:.0f} s; front v_n (2e3-2e4 s) median {res[k]['vn']:.4e} m/s, "
              f"Pi median {res[k]['Pi_med']:.3f} max {res[k]['Pi_max']:.3f}")
    kf, kc = list(runs)
    dt = (res[kf]["last_pi"] - res[kc]["last_pi"]) / res[kc]["last_pi"]
    print(f"vanishing time (last Pi > 0.5): fine - coarse = {res[kf]['last_pi'] - res[kc]['last_pi']:+.0f} s ({dt:+.4f})")
    print(f"front speed fine / coarse: {res[kf]['vn'] / res[kc]['vn']:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    fig, ax = plt.subplots(2, 2, figsize=(7.5, 5.2), constrained_layout=True)
    for k, r in runs.items():
        t, A = area(r)
        m = t <= DAY
        ax[0, 0].plot(t[m] / 3600, A[m] * 1e12, lw=1.4, label=k)
        ax[1, 0].plot(res[k]["tp"] / 3600, res[k]["ap"] * 1e12, "o-", ms=2.5, lw=1.2)
        mm = (res[k]["tr"] > 179) & (res[k]["tr"] <= DAY)
        ax[1, 1].semilogy(res[k]["tr"][mm] / 3600, res[k]["Pi"][mm], lw=1.0)
    ax[0, 1].plot(tg / 3600, 100 * rel, color="k", lw=1.2)
    ax[0, 0].set_ylabel("total ice area [µm²]")
    ax[0, 1].set_ylabel("(fine − coarse) / coarse [%]")
    ax[1, 0].set_ylabel("narrow-plug area [µm²]")
    ax[1, 1].set_ylabel("Π (K&P Eq. 45)")
    ax[1, 1].axhline(0.5, color="0.5", lw=0.8, ls="--")
    ax[1, 0].set_xlim(10, 18)
    ax[1, 0].set_ylim(-30, 900)
    for x in ax.flat:
        x.set_xlabel("time [h]")
    fig.legend(loc="outside upper center", ncol=2, fontsize=7, frameon=False)
    bad = check_figure(fig)
    fig.savefig(os.path.join(a.out, "issue1_compare.png"), dpi=200)
    print("figcheck:", bad or "OK", "->", os.path.join(a.out, "issue1_compare.png"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
