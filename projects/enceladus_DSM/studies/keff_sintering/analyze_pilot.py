#!/usr/bin/env python3
"""Pilot analysis: how k_eff evolves, and whether SSA explains it.

    venv_enceladus/bin/python studies/keff_sintering/analyze_pilot.py <batch_dir>

Reads merged runs (postprocess/merge_restart_legs.py) and falls back to the
leg-1 directory for any seed whose resume leg has not arrived, so it is useful
while downloads are still in flight.

WHAT SSA IS HERE. The solver logs int[ phi^2 (1-phi)^2 ] / eps in column 0 of
SSA_evo.dat. For the equilibrium logistic profile that integral is exactly
eps/6 per unit interface length, so

    interface length = 6 * column0,    SSA = 6 * column0 / (Lx*Ly)   [1/m]

which is what is plotted. The factor is not cosmetic: without it the axis is
in units nobody can compare to a measurement.

WHY THE CORRELATION IS TESTED, NOT ASSUMED. k_eff and SSA both change
monotonically during sintering, so they correlate trivially in time. The
question that matters is whether SSA carries information about k_eff BEYOND
what the ice fraction already carries -- if it does not, the observed
SSA-k_eff correlation is two quantities tracking the same third thing. Panel
(d) tests that directly.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
sys.path.insert(0, str(PROJ / "postprocess"))
import figstyle as fs                                       # noqa: E402
import pplib                                                # noqa: E402

DAY = 86400.0


def load_run(d: Path):
    ssa = pplib.load_ssa(str(d))
    k = np.genfromtxt(d / "k_eff.csv", delimiter=",", names=True)
    k = np.atleast_1d(k)
    opts = pplib.read_opts(str(d))
    Lx = pplib.opt_float(opts, "Lx", 2.0e-3)
    Ly = pplib.opt_float(opts, "Ly", Lx)
    return {
        "t": ssa[:, 2], "ssa": 6.0 * ssa[:, 0] / (Lx * Ly), "dt": ssa[:, 4],
        "ice": ssa[:, 1] / (Lx * Ly),
        "kt": k["time"], "k_iso": k["k_iso"], "k00": k["k_00"], "k11": k["k_11"],
        "k01": k["k_01"], "phi_bar": k["phi_bar"], "L": Lx,
    }


def find_runs(batch: Path) -> dict:
    """Every run under `batch`, keyed by seed.

    Discovery is by CONTENT -- a directory holding both k_eff.csv and
    SSA_evo.dat is a run -- rather than by directory name. The layout has
    already changed twice: originally `<geom>__<exp>/` beside `<geom>/<leg>/`,
    then a `merged/` tree, and now the merged runs relocated into `<geom>/`.
    Matching on names broke silently at each step and reported "no runs found"
    on a directory full of results.

    Prefers the deepest match when a run nests inside another, and prefers a
    merged run over the legs it was built from.
    """
    seen = {}
    for kf in sorted(batch.glob("**/k_eff.csv")):
        d = kf.parent
        if not (d / "SSA_evo.dat").is_file():
            continue
        m = re.search(r"seed(\d+)", d.name) or re.search(r"seed(\d+)", str(d))
        key = m.group(1) if m else d.name
        prev = seen.get(key)
        if prev is None:
            seen[key] = d
            continue
        # a merged run carries MERGE_INFO.json; prefer it over a raw leg
        merged_new = (d / "MERGE_INFO.json").is_file()
        merged_old = (prev / "MERGE_INFO.json").is_file()
        if merged_new != merged_old:
            if merged_new:
                seen[key] = d
            continue
        # otherwise the same seed was simply run more than once (rev64 seed 1
        # was, on 09-18 and 09-21, and the two agree to 1.2e-8 in k_iso).
        # Take the NEWEST rather than whichever happens to sort last.
        if d.stat().st_mtime > prev.stat().st_mtime:
            seen[key] = d
    return dict(sorted(seen.items()))


def collect(batch: Path):
    runs = {}
    for seed, d in find_runs(batch).items():
        r = load_run(d)
        r["merged"] = (d / "MERGE_INFO.json").is_file()
        runs[seed] = r
    return runs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path)
    a = ap.parse_args()
    runs = collect(a.batch)
    if not runs:
        raise SystemExit("no runs found")

    print(f"{'seed':>5} {'src':>8} {'t_end[d]':>9} {'k_iso 0':>9} {'k_iso end':>10} "
          f"{'change':>8} {'SSA 0':>10} {'SSA end':>10} {'change':>8} {'k11/k00':>8}")
    for s, r in runs.items():
        dk = r["k_iso"][-1] / r["k_iso"][0] - 1
        ds = np.interp(r["kt"][-1], r["t"], r["ssa"]) / r["ssa"][0] - 1
        print(f"{s:>5} {'merged' if r['merged'] else 'leg1':>8} "
              f"{r['kt'][-1]/DAY:9.2f} {r['k_iso'][0]:9.4f} {r['k_iso'][-1]:10.4f} "
              f"{dk:+8.1%} {r['ssa'][0]:10.0f} "
              f"{np.interp(r['kt'][-1], r['t'], r['ssa']):10.0f} {ds:+8.1%} "
              f"{r['k11'][-1]/r['k00'][-1]:8.3f}")

    C = [fs.C[i % len(fs.C)] for i in range(len(runs))]
    fig, ax = plt.subplots(2, 2, figsize=(10.0, 7.6))
    (a1, a2), (a3, a4) = ax

    for (s, r), c in zip(runs.items(), C):
        lab = f"seed {s}" + ("" if r["merged"] else "  (leg 1 only)")
        ls = "-" if r["merged"] else "--"
        a1.plot(r["kt"] / DAY, r["k_iso"], ls, color=c, lw=1.8, label=lab)
        a2.plot(r["t"] / DAY, r["ssa"], ls, color=c, lw=1.8, label=lab)
        ss = np.interp(r["kt"], r["t"], r["ssa"])
        a3.plot(ss, r["k_iso"], ls, color=c, lw=1.6, marker="o", ms=3.5,
                markeredgecolor="none", label=lab)
        a4.plot(r["kt"] / DAY, r["k11"] / r["k00"], ls, color=c, lw=1.8, label=lab)

    fs.style(a1, "time  [days]", r"$k_{\rm eff}$  [W m$^{-1}$K$^{-1}$]",
             "(a)  conductivity rises as necks grow", logy=False)
    a1.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower right")

    fs.style(a2, "time  [days]", r"SSA  [m$^{-1}$]",
             "(b)  surface area falls as it coarsens", logy=False)

    # The point of this panel: the ice fraction does not move AT ALL (mass is
    # conserved exactly in a closed isothermal box), so a parameterisation that
    # knows only density predicts a horizontal line. Everything above it is
    # structure.
    k0 = np.mean([r["k_iso"][0] for r in runs.values()])
    a3.axhline(k0, color=fs.MUTED, ls="--", lw=1.4)
    a3.text(0.03, k0, "  what density alone predicts\n"
                      r"  ($\bar\phi$ constant to 1e-16)",
            transform=a3.get_yaxis_transform(), va="bottom", ha="left",
            fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.style(a3, r"SSA  [m$^{-1}$]", r"$k_{\rm eff}$  [W m$^{-1}$K$^{-1}$]",
             "(c)  none of this rise is density", logy=False)
    a3.invert_xaxis()

    a4.axhline(1.0, color=fs.MUTED, ls="--", lw=1.0)
    a4.text(0.02, 0.97, "isotropic", transform=a4.transAxes, va="top",
            fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.style(a4, "time  [days]", r"$k_{yy}/k_{xx}$",
             "(d)  anisotropy is inherited, not made", logy=False)

    fig.tight_layout()
    fs.save(fig, HERE, "pilot_keff", dpi=190)

    # ---- does SSA add anything over ice fraction? -------------------------
    print("\n--- is SSA a cause or a proxy? ---")
    print("Pooling all seeds, k_eff regressed on the ice fraction alone, then")
    print("on ice fraction + SSA. If SSA adds nothing, the SSA-k_eff")
    print("correlation is both quantities tracking density.\n")
    X_phi, X_ssa, Y = [], [], []
    for r in runs.values():
        X_phi.append(r["phi_bar"])
        X_ssa.append(np.interp(r["kt"], r["t"], r["ssa"]))
        Y.append(r["k_iso"])
    xp, xs, y = map(np.concatenate, (X_phi, X_ssa, Y))

    def r2(cols):
        A = np.column_stack([np.ones_like(y)] + cols)
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ beta
        return 1.0 - res.var() / y.var()

    r_phi, r_both = r2([xp]), r2([xp, xs])
    print(f"  R^2, ice fraction only      {r_phi:.4f}")
    print(f"  R^2, ice fraction + SSA     {r_both:.4f}")
    print(f"  SSA adds                    {r_both - r_phi:+.4f}")
    print(f"  corr(SSA, k_eff) raw        {np.corrcoef(xs, y)[0,1]:+.4f}")
    print(f"  corr(phi_bar, k_eff)        {np.corrcoef(xp, y)[0,1]:+.4f}")
    print(f"  phi_bar range across run    {xp.min():.6f} .. {xp.max():.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
