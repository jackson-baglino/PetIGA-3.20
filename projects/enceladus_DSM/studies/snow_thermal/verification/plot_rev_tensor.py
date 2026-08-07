#!/usr/bin/env python3
"""Plot every run of the REV sweep as the four components of the k_eff tensor.

One figure, four panels laid out the way the tensor is written:

    +-----------+-----------+
    |   k_00    |   k_01    |
    +-----------+-----------+
    |   k_10    |   k_11    |
    +-----------+-----------+

Each panel carries every window size as its own line, k_ij(t).

WHY THE PANELS ARE PAIRED ON A SHARED Y-AXIS. The diagonal components are
O(0.5) W/m/K and the off-diagonals are O(0.03) -- more than an order of
magnitude apart. Putting all four on one common scale would flatten k_01 and
k_10 into the zero line; giving all four independent scales would make the
off-diagonal noise look as big as the diagonal signal. So the diagonals share
one scale with each other and the off-diagonals share another, which is the
comparison that actually means something in each case (k_00 vs k_11 is the
anisotropy; k_01 vs k_10 is the symmetry check).

k_01 and k_10 come from two independent corrector solves and are never forced
equal, so their difference is a free error estimate on the whole k_eff
calculation. The printed table reports it per run.

This is deliberately a companion to plot_rev.py, not a replacement: that script
answers "where is the REV?" from k_iso, this one shows the full tensor the REV
verdict is drawn from.

Runs on the cluster against $SCRATCH directly, or locally against the tree that
fetch_rev_keff.sh brings down (the k_eff.csv files are a few kB; the fields are
not downloadable).

Usage:
    python3 plot_rev_tensor.py --results <dir> [--tag rev] [--out DIR]
"""
import argparse
import csv
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# Window size is an ordered magnitude, not a set of unrelated categories, so it
# is encoded on a sequential ramp rather than with categorical hues: small L
# light, large L dark. viridis is monotone in lightness and colour-vision-safe
# by construction; the ends are trimmed because its palest yellow disappears on
# a light surface and its darkest purple disappears on a dark one.
RAMP = {"light": (0.78, 0.06), "dark": (0.22, 0.92)}
THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e",
                  muted="#8a8985", grid="#e3e2dd"),
    "dark":  dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7",
                  muted="#8a8985", grid="#333330"),
}
COMPONENTS = [("k_00", 0, 0), ("k_01", 0, 1), ("k_10", 1, 0), ("k_11", 1, 1)]


def style(ax, t):
    ax.set_facecolor(t["surface"])
    ax.grid(True, which="both", color=t["grid"], linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(t["grid"])
    ax.tick_params(colors=t["ink2"], labelsize=9)
    ax.xaxis.label.set_color(t["ink2"])
    ax.yaxis.label.set_color(t["ink2"])
    ax.title.set_color(t["ink"])


def collect(results, tag):
    """Find every k_eff.csv under `results`, newest run per geometry.

    The tree is walked rather than globbed at a fixed depth because the sweep
    has been written under two different layouts -- <base>/<geom>/<run>/ from
    run_enceladus.sh, and <base>/rev_sweep/<geom>/<run>/ when SCRATCH pointed a
    level higher. Walking finds the CSVs either way.
    """
    found = {}
    for dirpath, _dirnames, filenames in os.walk(results):
        if "k_eff.csv" not in filenames:
            continue
        run = os.path.basename(dirpath)
        geom = os.path.basename(os.path.dirname(dirpath))
        if tag and tag not in run and tag not in geom:
            continue
        # Newest wins. Run folders are timestamp-led (YYYY-MM-DD__HH.MM.SS_...)
        # precisely so that a lexicographic max is also the chronological one.
        if geom not in found or run > os.path.basename(found[geom]):
            found[geom] = dirpath

    rows = []
    for geom, run in found.items():
        with open(os.path.join(run, "k_eff.csv")) as fh:
            d = [r for r in csv.DictReader(fh) if r.get("time")]
        if not d:
            print(f"  empty k_eff.csv in {run} -- skipped")
            continue
        m = re.search(r"L([0-9.]+)mm", geom)
        # The packing's own solid fraction, so the run can be checked against
        # the microstructure it was supposed to start from. Without this the
        # phi_bar column has nothing to be wrong relative to.
        solid = np.nan
        if m:
            mp = os.path.join(PROJ, "inputs", "packings_rev",
                              f"rev_L{m.group(1)}mm", "metadata.json")
            if os.path.isfile(mp):
                with open(mp) as fh:
                    solid = 1.0 - json.load(fh)["porosity_achieved"]
        rows.append({
            "geom": geom,
            "run": run,
            "L": float(m.group(1)) if m else np.nan,
            "label": f"L = {float(m.group(1)):g} mm" if m else geom,
            "t": np.array([float(r["time"]) for r in d]),
            "k": {c: np.array([float(r[c]) for r in d]) for c, _, _ in COMPONENTS},
            "k_iso": np.array([float(r["k_iso"]) for r in d]),
            "phi_bar": np.array([float(r["phi_bar"]) for r in d]),
            "its": np.array([int(r["ksp_its"]) for r in d]),
            "reason": np.array([int(r["ksp_reason"]) for r in d]),
            "solid_pack": solid,
        })
    # NaN sorts last, so unlabelled geometries fall to the end instead of
    # raising or scrambling the colour order.
    rows.sort(key=lambda r: (np.isnan(r["L"]), r["L"], r["geom"]))
    return rows


def make_figure(rows, theme, out):
    t = THEMES[theme]
    lo, hi = RAMP[theme]
    colors = [plt.cm.viridis(x) for x in np.linspace(lo, hi, len(rows))]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6))
    fig.patch.set_facecolor(t["surface"])

    for name, i, j in COMPONENTS:
        ax = axes[i][j]
        style(ax, t)
        for r, c in zip(rows, colors):
            th = r["t"] / 3600.0
            if th.size == 1:
                ax.plot(th, r["k"][name], "o", color=c, ms=8, label=r["label"])
            else:
                ax.plot(th, r["k"][name], "-", lw=2.0, color=c, label=r["label"])
        ax.set_ylabel(r"$k_{%d%d}$ [W m$^{-1}$ K$^{-1}$]" % (i, j))
        diag = i == j
        # color= is not optional: set_title re-applies rcParams['text.color']
        # (black) over whatever style() put on ax.title, which silently paints
        # the dark-theme titles black-on-black.
        ax.set_title(f"$k_{{{i}{j}}}$ — {'diagonal' if diag else 'off-diagonal'}",
                     fontsize=11, loc="left", color=t["ink"])
        if not diag:
            ax.axhline(0.0, color=t["muted"], ls=":", lw=1.0)
        if i == 1:
            ax.set_xlabel("time [h]")

    # Pair the scales: diagonals together, off-diagonals together. Done by hand
    # rather than with sharey= so the two pairs stay independent of each other.
    for group in ((axes[0][0], axes[1][1]), (axes[0][1], axes[1][0])):
        lims = [a.get_ylim() for a in group]
        span = (min(l[0] for l in lims), max(l[1] for l in lims))
        for a in group:
            a.set_ylim(span)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=min(len(rows), 6),
               frameon=False, fontsize=9.5, labelcolor=t["ink2"],
               bbox_to_anchor=(0.5, 0.005))

    fig.suptitle("REV sweep — the four components of $\\mathbf{k}_{eff}$\n"
                 "nested centred windows of one 3.5 mm master packing, "
                 "1 day at $-20\\,^\\circ$C",
                 fontsize=13, color=t["ink"], x=0.008, ha="left", y=0.995,
                 va="top")
    fig.tight_layout(rect=[0, 0.055, 1, 0.965])
    p = os.path.join(out, f"rev_keff_tensor_{theme}.png")
    fig.savefig(p, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="directory holding the sweep. On the cluster: "
                         "/resnick/scratch/jbaglino/enceladus_DSM/rev_sweep . "
                         "Locally: the tree fetch_rev_keff.sh brought down.")
    ap.add_argument("--tag", default="rev",
                    help="only use runs whose folder name contains this "
                         "(pass '' to take everything found)")
    ap.add_argument("--out", default=HERE)
    args = ap.parse_args()

    rows = collect(args.results, args.tag)
    if not rows:
        raise SystemExit(f"no k_eff.csv found under {args.results} "
                         f"with tag '{args.tag}'")

    print("THE TENSOR AT t_end")
    print(f"{'window':>12} {'n':>4} {'t_end[h]':>9} {'k_00':>11} {'k_11':>11} "
          f"{'k_01':>11} {'k_10':>11} {'|k01-k10|/k_iso':>16} {'k_iso':>11}")
    summary = []
    for r in rows:
        f = {c: r["k"][c][-1] for c, _, _ in COMPONENTS}
        asym = abs(f["k_01"] - f["k_10"]) / r["k_iso"][-1]
        print(f"{r['label']:>12} {r['t'].size:4d} {r['t'][-1]/3600:9.2f} "
              f"{f['k_00']:11.4e} {f['k_11']:11.4e} {f['k_01']:11.4e} "
              f"{f['k_10']:11.4e} {asym:16.2e} {r['k_iso'][-1]:11.4e}")
        summary.append([r["geom"], r["L"], r["t"][-1], f["k_00"], f["k_01"],
                        f["k_10"], f["k_11"], r["k_iso"][-1],
                        r["phi_bar"][0], r["phi_bar"][-1], r["solid_pack"],
                        r["k_iso"][0], asym, int(r["its"][-1]),
                        int(r["reason"][-1]), r["run"]])

    # THE ICE-FRACTION AUDIT. k_eff in this model goes roughly as phi_bar^4, so
    # a couple of percent of ice fraction is ~10% of k_eff -- more than the REV
    # tolerance. A k_eff(L) trend can therefore be entirely an artefact of
    # phi_bar not being what it should be, in two distinct ways this table
    # separates:
    #   pack->t0   the run did not start from the packing it was given
    #              (rasterisation, diffuse-interface bias, a bad IC)
    #   t0->end    ice was created or destroyed during the run
    #              (boundary conditions, mass conservation)
    # Only once BOTH are small is a k_eff(L) trend about L at all.
    print("\nTHE ICE-FRACTION AUDIT — k_eff ~ phi_bar^4, so read this first")
    print(f"{'window':>12} {'pack':>8} {'phi(t0)':>9} {'phi(end)':>9} "
          f"{'pack->t0':>10} {'t0->end':>9} {'k_iso drift':>12}")
    for r in rows:
        p0, pe, pk = r["phi_bar"][0], r["phi_bar"][-1], r["solid_pack"]
        ic = 100 * (p0 - pk) / pk if np.isfinite(pk) else np.nan
        run = 100 * (pe - p0) / p0
        kd = 100 * (r["k_iso"][-1] - r["k_iso"][0]) / r["k_iso"][0]
        print(f"{r['label']:>12} {pk:8.4f} {p0:9.4f} {pe:9.4f} "
              f"{ic:+9.2f}% {run:+8.2f}% {kd:+11.2f}%")

    ic = [abs(r["phi_bar"][0] - r["solid_pack"]) / r["solid_pack"]
          for r in rows if np.isfinite(r["solid_pack"])]
    if ic:
        note = "" if max(ic) < 0.02 else (
            "  <-- the ICs do not match the packings, so k_eff(L) "
            "is not about L until this is explained")
        print(f"\nworst pack->t0 ice-fraction error: {max(ic)*100:.1f}%{note}")

    worst = max(s[12] for s in summary)
    print(f"largest k_01 vs k_10 asymmetry: {worst:.2e} of k_iso "
          f"({'fine — the corrector solves agree' if worst < 1e-6 else 'the corrector solves are stopping short; check ksp_reason'})")

    csvp = os.path.join(args.out, "rev_keff_tensor_summary.csv")
    with open(csvp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["geometry", "L_mm", "t_end_s", "k_00", "k_01", "k_10",
                    "k_11", "k_iso_end", "phi_bar_t0", "phi_bar_end",
                    "solid_frac_packing", "k_iso_t0", "asym_over_k_iso",
                    "ksp_its_end", "ksp_reason_end", "run_dir"])
        w.writerows(summary)
    print("wrote", csvp)

    for theme in ("light", "dark"):
        print("wrote", make_figure(rows, theme, args.out))


if __name__ == "__main__":
    main()
