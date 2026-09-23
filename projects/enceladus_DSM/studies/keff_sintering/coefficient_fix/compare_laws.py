#!/usr/bin/env python3
"""Compare k_eff under the band-interpolation laws, on a real packing.

    venv_enceladus/bin/python studies/keff_sintering/coefficient_fix/compare_laws.py <batch_dir>

WHAT THIS IS FOR, AND WHY THE LADDERS DO NOT ALREADY ANSWER IT
--------------------------------------------------------------
`studies/keff_sharp_limit/` verifies the tensor law on two geometries with
closed-form answers: a planar slab (tensor k_perp exact to 4e-7) and an
isolated disk array (tensor error +1.2% -> +0.006%, observed order ~2.5).
Neither has a CONTACT.

The packing does. Its grains are seated in exact tangency, so at t = 0 the
conduction path between two grains is a single point. That is precisely where
the arithmetic law does its worst: the diffuse band bridges the contact with
partially-conducting material, inflating k_eff(0) for a structure that ought
to barely conduct, and thereby suppressing the RELATIVE rise -- the campaign's
headline. The disk ladder's bias at the pilot's eps/R = 0.02 is about +8%, and
that is a LOWER BOUND on the packing, not an estimate of it. This script
measures the real thing.

INPUTS. One directory per run, each holding the in-line arith `k_eff.csv` from
the original run plus whatever `k_eff_<law>.csv` files
`scripts/HPC/submit_keff_replay.sh` has produced. Discovery is by content, as
in analyze_pilot.py, because the layout has moved more than once.

LEGS. A run killed on walltime and resumed has its snapshots split across two
directories, and a per-leg replay therefore writes two CSVs. They are
concatenated here by time, with leg-1 rows at or after the join dropped -- the
same rule postprocess/merge_restart_legs.py applies to the in-line CSV.

Writes compare_laws.csv and compare_laws.png next to this file.
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
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "preprocess"))
import figstyle as fs                                      # noqa: E402

DAY = 86400.0

# Plot order and colour, coarsest-biased first so the figure reads as a
# correction. `arith` is the legacy law every run before 2026-09-23 used.
LAWS = ("arith", "tensor", "sharp")
COLOR = {"arith": "#c0392b", "tensor": "#1f6fb4", "sharp": "#2e8b57"}
LABEL = {"arith": "arith (legacy)", "tensor": "tensor", "sharp": "sharp"}


def law_of(path: Path) -> str:
    """`k_eff.csv` is arith; `k_eff_<law>.csv` names its own law."""
    m = re.fullmatch(r"k_eff(?:_(\w+))?\.csv", path.name)
    if not m:
        return ""
    return m.group(1) or "arith"


def seed_of(d: Path) -> str:
    m = re.search(r"seed(\d+)", d.name) or re.search(r"seed(\d+)", str(d))
    return m.group(1) if m else d.name


def collect(batch: Path) -> dict:
    """{seed: {law: [ (mtime, recarray), ... ]}} -- one entry per leg."""
    out: dict = {}
    for csv in sorted(batch.glob("**/k_eff*.csv")):
        law = law_of(csv)
        if not law:
            continue
        try:
            a = np.atleast_1d(np.genfromtxt(csv, delimiter=",", names=True))
        except Exception as exc:                       # noqa: BLE001
            print(f"  ⚠ unreadable, skipping: {csv} ({exc})")
            continue
        if a.size == 0 or "k_iso" not in (a.dtype.names or ()):
            continue
        out.setdefault(seed_of(csv.parent), {}).setdefault(law, []).append(a)
    return dict(sorted(out.items()))


def merge_legs(parts: list) -> np.ndarray:
    """Concatenate legs by time, dropping the earlier leg's superseded tail.

    A resume restarts from the last good snapshot, so the leg it replaces has
    rows at and after that time which the new leg recomputes. Keeping both
    would put two points at one time and bend the curve.
    """
    parts = sorted(parts, key=lambda a: float(np.min(a["time"])))
    kept = [parts[0]]
    for nxt in parts[1:]:
        join = float(np.min(nxt["time"]))
        kept[-1] = kept[-1][kept[-1]["time"] < join]
        kept.append(nxt)
    merged = np.concatenate(kept)
    return merged[np.argsort(merged["time"])]


# The initial condition is an analytic sum of tanh profiles, not an equilibrated
# phase field, so the first hours of every run are the field relaxing onto its
# equilibrium profile rather than sintering. Measured on pilot seed 1: the
# log-log slope of SSA(t) runs 0.00 -> -0.03 -> -0.075 over the first ~8 h and
# is flat at ~-0.077 thereafter. Anything read at t = 0 is therefore a property
# of the initial condition, not of the model.
BASELINE_DAYS = 1.0


def rise(a: np.ndarray, baseline_days: float = BASELINE_DAYS) -> tuple:
    """(k at baseline, k at t_end, % rise between them, actual baseline time).

    The baseline is the first sample at or after `baseline_days`, NOT t = 0.
    The choice moves the headline by several points -- on the pilot's arith
    curves the ensemble rise is +19.6% from t=0, +17.7% from 0.5 d, +15.9%
    from 1 d -- so it has to be stated, not defaulted into.
    """
    t = np.asarray(a["time"], dtype=float)
    k = np.asarray(a["k_iso"], dtype=float)
    idx = np.flatnonzero(t >= baseline_days * DAY)
    i = int(idx[0]) if idx.size else 0
    return float(k[i]), float(k[-1]), 100.0 * (k[-1] / k[i] - 1.0), float(t[i])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path, help="downloaded batch directory to scan")
    ap.add_argument("--out", type=Path, default=HERE)
    ap.add_argument("--dpi", type=int, default=160)
    ap.add_argument("--baseline-days", type=float, default=BASELINE_DAYS,
                    help="measure the rise from the first sample at or after "
                         "this time, not from t=0 (the IC has not relaxed at t=0)")
    args = ap.parse_args()

    if not args.batch.is_dir():
        print(f"not a directory: {args.batch}", file=sys.stderr)
        return 2

    runs = collect(args.batch)
    if not runs:
        print(f"no k_eff*.csv found under {args.batch}", file=sys.stderr)
        return 1

    data: dict = {}
    for seed, bylaw in runs.items():
        data[seed] = {law: merge_legs(parts) for law, parts in bylaw.items()}

    present = sorted({law for d in data.values() for law in d},
                     key=lambda l: LAWS.index(l) if l in LAWS else 99)
    print(f"seeds: {', '.join(data)}      laws present: {', '.join(present)}")
    if present == ["arith"]:
        print("\nOnly the legacy arith law is present -- the replay has not been\n"
              "run, or its CSVs have not been downloaded. Submit it with:\n"
              "  ./scripts/HPC/submit_keff_replay.sh --dry-run \\\n"
              "      --roots <cluster batch parent> <cluster geom dir>/")

    # ---- table -----------------------------------------------------------
    rows = []
    b = args.baseline_days
    print(f"\nrise measured from the first sample at or after t = {b:g} d "
          f"(NOT t=0: the IC has not relaxed there)")
    print(f"\n{'seed':>5} {'law':>8} {'t_b[d]':>7} {'k(t_b)':>9} {'k(end)':>9} "
          f"{'rise':>9} {'t_end[d]':>9} {'samples':>8} {'ksp_its':>8}")
    for seed in data:
        for law in present:
            a = data[seed].get(law)
            if a is None:
                continue
            kb, kN, r, tb = rise(a, b)
            tend = float(a["time"][-1]) / DAY
            its = float(np.median(a["ksp_its"])) if "ksp_its" in a.dtype.names else float("nan")
            print(f"{seed:>5} {law:>8} {tb/DAY:7.2f} {kb:9.4f} {kN:9.4f} {r:+8.1f}% "
                  f"{tend:9.2f} {len(a):8d} {its:8.0f}")
            rows.append((seed, law, kb, kN, r, tend, len(a), its, tb / DAY))

    # ---- ensemble and the two things worth checking ----------------------
    print()
    ens = {}
    for law in present:
        rs = [row[4] for row in rows if row[1] == law]
        if rs:
            ens[law] = (float(np.mean(rs)), float(np.std(rs, ddof=1)) if len(rs) > 1 else 0.0)
            sd, n = ens[law][1], len(rs)
            print(f"  {law:>8}: rise {ens[law][0]:+.1f}%  sd {sd:.1f} pts  "
                  f"SEM {sd/max(np.sqrt(n),1):.1f}  (n={n})")

    if "arith" in ens and "tensor" in ens:
        print(f"\n  the correction: {ens['arith'][0]:+.1f}% -> {ens['tensor'][0]:+.1f}%  "
              f"(x{ens['tensor'][0]/ens['arith'][0]:.2f} on the headline)")
    if "tensor" in ens and "sharp" in ens:
        # The cross-check. These are completely different routes to removing
        # the same bias -- one analytic, one by thresholding the geometry --
        # so agreement is the evidence that either is right, and disagreement
        # is a result in its own right.
        gap = [(s, 100.0 * (float(data[s]["sharp"]["k_iso"][-1]
                            / data[s]["tensor"]["k_iso"][-1]) - 1.0))
               for s in data if "sharp" in data[s] and "tensor" in data[s]]
        if gap:
            worst = max(abs(g) for _, g in gap)
            print(f"  tensor vs sharp at t_end: max |diff| {worst:.2f}% "
                  f"across {len(gap)} seed(s)" + ("  ✓" if worst < 5.0 else "  ⚠ >5%"))

    # ---- CSV -------------------------------------------------------------
    outcsv = args.out / "compare_laws.csv"
    with outcsv.open("w") as fh:
        fh.write("seed,law,baseline_days,k_iso_baseline,k_iso_end,rise_pct,"
                 "t_end_days,n_samples,ksp_its_median\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[8]:.4f},{r[2]:.6e},{r[3]:.6e},{r[4]:.4f},"
                     f"{r[5]:.4f},{r[6]},{r[7]:.0f}\n")

    # ---- figure ----------------------------------------------------------
    # Panel (a) is absolute, panel (b) normalised to each curve's own t = 0.
    # Both are needed: the law moves the LEVEL and the RISE by different
    # factors, and the campaign reports the rise.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for law in present:
        for i, seed in enumerate(data):
            a = data[seed].get(law)
            if a is None:
                continue
            t = a["time"] / DAY
            kb, _, _, tb = rise(a, b)
            axes[0].plot(t, a["k_iso"], color=COLOR.get(law, "k"), lw=1.3, alpha=0.85,
                         label=LABEL.get(law, law) if i == 0 else None)
            # Normalised to the BASELINE sample, not to t = 0. Points before it
            # are drawn dotted: they are the initial condition relaxing onto an
            # equilibrium profile, not sintering, and they are not part of the
            # reported rise.
            pre = t < tb / DAY
            axes[1].plot(t[pre], a["k_iso"][pre] / kb, color=COLOR.get(law, "k"),
                         lw=1.0, alpha=0.5, ls=":")
            axes[1].plot(t[~pre], a["k_iso"][~pre] / kb, color=COLOR.get(law, "k"),
                         lw=1.3, alpha=0.85,
                         label=LABEL.get(law, law) if i == 0 else None)
    axes[0].set_ylabel(r"$k_{\mathrm{eff}}$  [W m$^{-1}$ K$^{-1}$]")
    axes[1].set_ylabel(rf"$k_{{\mathrm{{eff}}}}(t)\,/\,k_{{\mathrm{{eff}}}}(t_0)$,  $t_0={b:g}$ d")
    axes[1].axhline(1.0, color="0.6", lw=0.8, ls=":")
    axes[1].axvline(b, color="0.6", lw=0.8, ls="--")
    for ax, ttl in zip(axes, ("(a) absolute",
                              "(b) relative to each law's own baseline "
                              "(dotted = IC relaxation, excluded)")):
        ax.set_xlabel("time [days]")
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.25, lw=0.6)
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle("k_eff under each band-interpolation law — pilot packing, one curve per seed",
                 fontsize=11)
    fig.tight_layout()
    out = args.out / "compare_laws.png"
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print(f"\nwrote {outcsv}\n      {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
