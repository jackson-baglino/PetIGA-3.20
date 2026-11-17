#!/usr/bin/env python3
"""analyze_type5.py: lens position and area in the type 5 runs, against the reduced model.

For each run folder (a batch folder is expanded to its run folders whose geometry starts with
t5_cell or t5_pore), reads ice_bodies.csv (postprocess/ice_bodies.py, run by
run_postprocess.sh): the largest body per step, its phi = 0.5 crossings on the pore axis
(face_lo, face_hi) and its area. The lens center is the mean of the two crossings, given as the
distance from the throat (xi = x - Lx/2; positive = cold side, x > Lx/2, since x = 0 is the warm
end in every gradient run). The model prediction comes from type5_model.py with the run's own
settings: the periodic cell (5P) for -periodic_x 1 runs, the ends model (5G, 5S) otherwise, with
h taken from the run's -rhovfix_lo (h = frac / (rho_vs(T_warm)/rho_vs(T0)), rhov_bc.fracs) and
G from -grad_temp0. Water closure is read from water_budget.csv when present.

Writes <out>/type5_summary.csv and <out>/type5_lens.png (center and area against time, run and
model), and prints one line per run.

    venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_type5.py <batch or run folders...> [--out DIR]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "postprocess"))
from pplib import opt_float, read_opts  # noqa: E402
from rhov_bc import fracs  # noqa: E402
from type5_model import DAY, EndsPore, PeriodicCell, consts, side_words  # noqa: E402


def run_folders(paths):
    out = []
    for p in map(Path, paths):
        if (p / "ice_bodies.csv").exists() or (p / "SSA_evo.dat").exists():
            out.append(p)
        else:
            out += sorted(q for q in p.iterdir() if q.is_dir() and q.name.startswith(("t5_cell", "t5_pore")))
    return out


def settings(run: Path):
    """dict: periodic, wt (m; = W for a straight pore), W, Lx, G (K/m, warm at x = 0), h, T0."""
    o = read_opts(str(run))
    lx = float(opt_float(o, "-Lx"))
    W = float(opt_float(o, "-Ly"))
    T0 = float(opt_float(o, "-temp", -20.0))
    periodic = int(float(opt_float(o, "-periodic_x", 0) or 0)) == 1
    geo = next((f.name for f in run.glob("t5_*.opts") if "_T-20_" not in f.name), run.name)
    m = re.search(r"_wt(\d+(?:\.\d+)?)um", geo)
    wt = float(m.group(1)) * 1e-6 if m else W
    g0 = o.get("-grad_temp0", "0,0,0")
    G = -float(str(g0).split(",")[0])
    h = None
    if not periodic and int(float(opt_float(o, "-flag_BC_rhovfix", 0) or 0)) == 1:
        lo = float(opt_float(o, "-rhovfix_lo"))
        unit = fracs(T0, -G, lx, "h-face", 1.0)[0]
        h = lo / unit
    return {"periodic": periodic, "wt": wt, "W": W, "Lx": lx, "G": G, "h": h, "T0": T0, "geo": geo}


def lens_track(run: Path, lx: float):
    """arrays t [s], xi [m] (center from the throat), A/A0, from the largest body per step."""
    f = run / "ice_bodies.csv"
    if not f.exists():
        raise FileNotFoundError(f"{f} missing: run bash {run}/postprocess/run_postprocess.sh {run} first")
    best = {}
    with f.open() as fh:
        for r in csv.DictReader(fh):
            s = int(r["step"])
            a = float(r["area_m2"])
            if s not in best or a > best[s][1]:
                best[s] = (float(r["t_s"]), a, float(r["face_lo"]), float(r["face_hi"]))
    rows = [best[s] for s in sorted(best)]
    t = np.array([r[0] for r in rows])
    A = np.array([r[1] for r in rows])
    xc = 0.5 * (np.array([r[2] for r in rows]) + np.array([r[3] for r in rows]))
    return t, xc - 0.5 * lx, A / A[0]


def model_track(st, lens0, start0, days):
    k = consts(st["T0"], 1e-2)
    if st["periodic"]:
        cell = PeriodicCell(st["wt"], k, W=st["W"], L=st["Lx"], lens=lens0, start=start0)
        traj = cell.run(days)
        return np.array([q[0] for q in traj]), np.array([q[1] for q in traj]), np.ones(len(traj))
    pore = EndsPore(st["wt"], k, W=st["W"], L=st["Lx"], lens=lens0, start=start0)
    r = pore.run(st["h"], st["G"], days, path=True)
    tr = r["traj"]
    return np.array([q[0] for q in tr]), np.array([q[1] for q in tr]), np.array([q[2] for q in tr])


def water_closure(run: Path):
    """largest |closure| over the run relative to the initial ice mass (water_budget.csv,
    postprocess/water_budget.py); nan when the file is missing."""
    f = run / "water_budget.csv"
    if not f.exists():
        return float("nan")
    with f.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return float("nan")
    m0 = float(rows[0]["ice_kg_per_m"]) + float(rows[0]["vapor_kg_per_m"])
    return max(abs(float(r["closure_kg_per_m"])) for r in rows) / m0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--out", default=None, help="output folder (default: the first path)")
    ap.add_argument("--lens", type=float, default=1.0e-4, help="initial lens length on the axis [m]")
    ap.add_argument("--start", type=float, default=7.0e-5, help="initial center, cold side of the throat [m]")
    a = ap.parse_args()
    runs = run_folders(a.paths)
    if not runs:
        sys.exit("analyze_type5: no type 5 run folders found")
    out = Path(a.out or a.paths[0])
    rows, tracks = [], []
    for run in runs:
        st = settings(run)
        try:
            t, xi, ar = lens_track(run, st["Lx"])
        except FileNotFoundError as e:
            print(f"SKIP {run.name}: {e}")
            continue
        days = t[-1] / DAY
        tm, xm, am = model_track(st, a.lens, a.start, max(days, 1e-3))
        x_model_end = float(np.interp(t[-1], tm, xm))
        kind = "5P" if st["periodic"] else ("5S" if st["wt"] >= st["W"] else "5G")
        row = {"run": run.name, "kind": kind, "wt_um": st["wt"] * 1e6, "G_K_per_m": st["G"],
               "h_end": "" if st["h"] is None else f"{st['h']:.10f}", "days": f"{days:.2f}",
               "xi_start_um": f"{xi[0] * 1e6:.2f}", "xi_end_um": f"{xi[-1] * 1e6:.2f}",
               "xi_model_end_um": f"{x_model_end * 1e6:.2f}", "area_end_rel": f"{ar[-1] - 1:+.4f}",
               "area_maxdev": f"{np.max(np.abs(ar - 1)):.4f}", "water_closure": f"{water_closure(run):.2e}"}
        rows.append(row)
        tracks.append((row, t, xi, ar, tm, xm, am))
        print(f"{kind} wt {st['wt'] * 1e6:5.1f} um, G {st['G']:5.2f} K/m, {days:6.1f} d: lens center "
              f"{side_words(xi[-1])} (model {side_words(x_model_end)}); area {ar[-1] - 1:+.2%}  {run.name}")
    if not rows:
        return 1
    out.mkdir(parents=True, exist_ok=True)
    with (out / "type5_summary.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 8})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 6.0), sharex=True)
    for k, (row, t, xi, ar, tm, xm, am) in enumerate(tracks):
        c = plt.cm.tab10(k % 10)
        lab = f"{row['kind']} wt {row['wt_um']:.0f}, G {row['G_K_per_m']:.2f}"
        ax1.plot(t / DAY, xi * 1e6, color=c, lw=1.2, label=lab)
        ax1.plot(tm / DAY, xm * 1e6, color=c, lw=0.9, ls="--")
        ax2.plot(t / DAY, (ar - 1) * 100, color=c, lw=1.2)
        if row["kind"] != "5P":
            ax2.plot(tm / DAY, (am - 1) * 100, color=c, lw=0.9, ls="--")
    ax1.axhline(0, color="0.5", lw=0.6)
    ax1.set_ylabel("lens center from the throat [µm]\n(+ cold side)")
    ax1.set_title("solid: simulation; dashed: type5_model", fontsize=8)
    ax1.legend(fontsize=7, ncol=2, frameon=False)
    ax2.set_ylabel("lens area change [%]")
    ax2.set_xlabel("time [days]")
    fig.tight_layout()
    fig.savefig(out / "type5_lens.png", dpi=200)
    print(f"wrote {out / 'type5_summary.csv'} and {out / 'type5_lens.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
