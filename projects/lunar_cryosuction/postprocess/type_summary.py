#!/usr/bin/env python3
"""One CSV and one figure per simulation type (simulation plan tool T6).

A manifest lists the runs of a type, one per line:

    # type  name   run folder (restart legs joined with '+')
    4       4A     /path/to/run
    4       4B     /path/to/legA+/path/to/legB

For every run this collects (computing what is missing):

    water_budget  ice left, ice lost, vapor change, outflow, closure error,
                  recapture fraction, time to lose 50 % and 90 % of the ice (M11)
    ice_bodies    bodies at the end, events, final face positions and speeds
    resolution    largest Pi after the start-up window (2 tau_sub) and the
                  largest phi seen on a reservoir face (resolution.dat)

and writes <out>/type<T>_summary.csv and <out>/type<T>_summary.png (ice
fraction left, upper-face depth and closure error against time, one line per
run). The figure is checked with figcheck before it is saved; problems are
printed and make the exit status 1.

ice_bodies.csv in a run folder is reused when it is newer than every snapshot;
otherwise ice_bodies is run (with --every, default 1) and its CSV written into
--cache (default: the output folder), never into the run folder.

Usage:
    python3 postprocess/type_summary.py --manifest runs_type4.txt --out <folder>
        [--every 2] [--vertical] [--no-recapture]
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ice_bodies  # noqa: E402
import pplib  # noqa: E402
import water_budget  # noqa: E402

COLUMNS = ["type", "name", "t_end_s", "ice_left_frac", "ice_lost_kg_per_m", "vapor_change_kg_per_m",
           "outflow_kg_per_m", "closure_rel", "recapture", "t50_s", "t90_s", "bodies_end", "events",
           "face_lo_end_m", "v_face_lo_end_m_per_s", "face_hi_end_m", "v_face_hi_end_m_per_s",
           "pi_max", "phi_reservoir_max", "folder"]


def read_manifest(path):
    runs = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split(None, 2)
            if len(parts) != 3:
                raise SystemExit(f"{path}:{ln}: expected 'type name folder', got {line!r}")
            legs = [os.path.expanduser(p) for p in parts[2].split("+")]
            for leg in legs:
                if not os.path.isdir(leg):
                    raise SystemExit(f"{path}:{ln}: no such folder {leg}")
            runs.append((parts[0], parts[1], legs))
    if not runs:
        raise SystemExit(f"{path}: no runs listed")
    return runs


def time_to_lose(t, lost, total, frac):
    """First time the ice lost reaches frac of the initial ice (NaN if never)."""
    if total <= 0:
        return np.nan
    target = frac * total
    k = np.flatnonzero(lost >= target)
    if k.size == 0:
        return np.nan
    k = k[0]
    if k == 0:
        return t[0]
    f = (target - lost[k - 1]) / (lost[k] - lost[k - 1])
    return t[k - 1] + f * (t[k] - t[k - 1])


def resolution_summary(run_dir):
    path = os.path.join(run_dir, "resolution.dat")
    if not os.path.isfile(path):
        return np.nan, np.nan
    a = np.genfromtxt(path, comments="#", ndmin=2, invalid_raise=False)
    if a.size == 0:
        return np.nan, np.nan
    tau = np.nan
    outp = os.path.join(run_dir, "outp.txt")
    if os.path.isfile(outp):
        with open(outp, errors="replace") as fh:
            m = re.search(r"tau_sub\s*=\s*([0-9.eE+-]+)\s*s", fh.read())
            if m:
                tau = float(m.group(1))
    skip = 2.0 * tau if np.isfinite(tau) else 0.0
    late = a[a[:, 0] - a[0, 0] >= skip]
    pi = late[:, 4].max() if late.size else np.nan
    return float(pi), float(a[:, 5].max())


def bodies_for(legs, cache, every):
    """ice_bodies rows for a run: reuse a fresh CSV or compute."""
    csv_path = os.path.join(legs[0], "ice_bodies.csv")
    snaps = [os.path.join(d, f) for d in legs for f in os.listdir(d)
             if re.fullmatch(r"sol_\d+\.dat", f)]
    newest = max((os.path.getmtime(f) for f in snaps), default=0.0)
    if len(legs) == 1 and os.path.isfile(csv_path) and os.path.getmtime(csv_path) >= newest:
        with open(csv_path) as fh:
            rows = [r for r in csv.DictReader(fh)]
        for r in rows:
            for k, v in r.items():
                try:
                    r[k] = float(v)
                except (TypeError, ValueError):
                    pass
        events = []
        ev = os.path.join(legs[0], "ice_bodies_events.txt")
        if os.path.isfile(ev):
            events = [ln for ln in open(ev) if not ln.startswith("#")]
        return rows, len(events)
    rows, contacts, events, _ = ice_bodies.analyze(legs, every=every)
    ice_bodies.add_rates(rows)
    os.makedirs(cache, exist_ok=True)
    tag = os.path.basename(os.path.normpath(legs[0]))
    ice_bodies.write_csv(os.path.join(cache, f"ice_bodies_{tag}.csv"), rows, ice_bodies.FIELDS)
    return rows, len(events)


def summarize(typ, name, legs, cache, every, do_recapture):
    b = water_budget.budget(legs)
    rec = water_budget.recapture(legs, "auto", 1) if do_recapture else None
    rows, n_events = bodies_for(legs, cache, every)
    m0 = b["m_ice"][0]
    lost = b["lost"]
    end = {}
    if rows:
        t_last = max(r["t_s"] for r in rows)
        last = [r for r in rows if r["t_s"] == t_last]
        # The deepest-reaching body's faces (largest area if tied) stand for the run.
        main = max(last, key=lambda r: r["area_m2"])
        end = dict(bodies_end=len(last), face_lo_end_m=main.get("face_lo"),
                   v_face_lo_end_m_per_s=main.get("v_face_lo"), face_hi_end_m=main.get("face_hi"),
                   v_face_hi_end_m_per_s=main.get("v_face_hi"))
    pi, phi_res = resolution_summary(legs[-1])
    out = dict(type=typ, name=name, t_end_s=b["t"][-1],
               ice_left_frac=b["m_ice"][-1] / m0 if m0 else np.nan,
               ice_lost_kg_per_m=lost[-1], vapor_change_kg_per_m=b["vap"][-1] - b["vap"][0],
               outflow_kg_per_m=b["out"][-1],
               closure_rel=np.abs(b["closure"]).max() / max(np.abs(lost).max(), 1e-300),
               recapture=(rec[1] / rec[0]) if rec and rec[0] > 0 else np.nan,
               t50_s=time_to_lose(b["t"], lost, m0, 0.5), t90_s=time_to_lose(b["t"], lost, m0, 0.9),
               events=n_events, pi_max=pi, phi_reservoir_max=phi_res, folder="+".join(legs), **end)
    return out, b, rows


def figure(typ, results, path, vertical):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.8), constrained_layout=True)
    tmax = max(b["t"][-1] for _, b, _ in results)
    unit = pplib.auto_time_unit(tmax)
    for s, b, rows in results:
        t = pplib.in_time_unit(b["t"], unit)
        ln, = axs[0].plot(t, b["m_ice"] / b["m_ice"][0], lw=1.4, label=s["name"])
        main = {}
        for r in rows:
            if np.isfinite(r.get("face_lo", np.nan)):
                if r["t_s"] not in main or r["area_m2"] > main[r["t_s"]]["area_m2"]:
                    main[r["t_s"]] = r
        if main:
            tt = sorted(main)
            axs[1].plot(pplib.in_time_unit(tt, unit), [main[k]["face_lo"] * 1e6 for k in tt],
                        lw=1.4, color=ln.get_color())
        ref = max(np.abs(b["lost"]).max(), 1e-300)
        axs[2].plot(t, b["closure"] / ref, lw=1.2, color=ln.get_color())
    axs[0].set_ylabel("ice left (fraction)")
    axs[1].set_ylabel("upper face x [µm]")
    if vertical:
        axs[1].invert_yaxis()
        axs[1].set_ylabel("upper face depth x [µm]")
    axs[2].set_ylabel("closure / ice lost")
    for ax in axs:
        ax.set_xlabel(f"time [{unit}]")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.legend(loc="outside lower center", ncol=min(6, len(results)), frameon=False, fontsize=8)
    fig.suptitle(f"Simulation type {typ}", fontsize=10)
    problems = check_figure(fig, min_pt=7.0)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return problems


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cache", default=None)
    p.add_argument("--every", type=int, default=1)
    p.add_argument("--vertical", action="store_true")
    p.add_argument("--no-recapture", action="store_true")
    a = p.parse_args(argv)
    runs = read_manifest(a.manifest)
    os.makedirs(a.out, exist_ok=True)
    cache = a.cache or a.out
    by_type = {}
    for typ, name, legs in runs:
        by_type.setdefault(typ, []).append(summarize(typ, name, legs, cache, a.every, not a.no_recapture))
    status = 0
    for typ, results in by_type.items():
        path = os.path.join(a.out, f"type{typ}_summary.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
            w.writeheader()
            for s, _, _ in results:
                w.writerow({k: (f"{v:.6e}" if isinstance(v, (float, np.floating)) else v) for k, v in s.items()})
        probs = figure(typ, results, os.path.join(a.out, f"type{typ}_summary.png"), a.vertical)
        print(f"type {typ}: {len(results)} runs -> {path}")
        for pr in probs:
            print(f"  figure check: {pr}")
            status = 1
    return status


if __name__ == "__main__":
    raise SystemExit(main())
