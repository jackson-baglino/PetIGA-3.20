#!/usr/bin/env python3
"""Ice bodies in every snapshot: area, center of mass, extent, faces on the
pore axis, contact lines, and their velocities and accelerations
(simulation plan tool T1; measurements M1, M4, M5, M6).

WHAT A "BODY" IS

A body is a connected region of phi > 0.5 (4-connected: two cores that touch
only at a corner are two bodies). Every node of the grid is then given to the
nearest core, so the diffuse band around a core belongs to it, and the body's
area is the integral of phi over its share of the grid. The areas of all
bodies therefore add up to the integral of phi over the domain, which is
tot_ice in SSA_evo.dat up to the sampling error (checked at 6 digits on the
stage-0c channel run). Cores smaller than --min-area are not bodies (noise or
a body about to vanish); their share goes to the nearest real body.

WHAT IS MEASURED PER BODY AND SNAPSHOT

  area_m2         integral of phi over the body's share [m^2 per unit depth]
  area_core_m2    area of phi > 0.5 (the sharp-interface area)
  com_x, com_y    center of mass, weights phi [m]
  xmin..ymax      extent of the phi > 0.5 core (interpolated crossings on the
                  axis line, node positions across it) [m]
  face_lo, face_hi  phi = 0.5 crossings on the pore axis line y = --axis-y
                  (default mid-height) that bound this body, refined with the
                  tanh fit (pplib.refine_tanh); NaN when the body does not
                  reach the axis line
  n_contacts      phi = 0.5 crossings on the wall faces (-wall_faces) that
                  belong to this body; their positions go to the contacts CSV

Velocities and accelerations of com_x, com_y, face_lo, face_hi and area come
from second-order finite differences on the (uneven) snapshot times
(np.gradient with the times as coordinates). An acceleration needs three
snapshots; with fewer it is NaN.

TRACKING

Bodies are matched between consecutive snapshots by the overlap of their
cores; the body that keeps most of an old core keeps its id. A core that
overlaps no previous core is new; one that overlaps two is a merge; two cores
sharing one parent is a split. These events are written to the events file,
and a body id never jumps from one physical body to another: after a merge
the merged body gets a new id, so its velocity history starts fresh.

FIXED-GAP SAMPLING (--at)

Velocities change as a run proceeds, so comparing runs at "the second half of
the run" mixes different states (the reason for tool T8). --at COL=VALUE
reports every quantity at the moments when COL of a body crosses VALUE,
linearly interpolated in time between the two snapshots that bracket the
crossing. Example: --at face_lo=2e-4 gives each face speed at the moment the
upper face passes x = 200 um; --at gap=1e-4 uses the gap between two bodies
on the axis (bodies 0 and 1 by default; see --gap-bodies).

OUTPUT (in <run>/ unless --out)

  ice_bodies.csv           one row per (snapshot, body)
  ice_bodies_contacts.csv  one row per contact-line crossing
  ice_bodies_events.txt    births, deaths, merges and splits
  ice_bodies_at.csv        with --at
  plots/ice_bodies.png     area and face positions vs time (with --plot)

Usage:
    python3 postprocess/ice_bodies.py --dir <run> [--dir <restart leg> ...]
        [--source auto|nurbs|vts] [--axis-y Y] [--at face_lo=2e-4] [--plot]
        [--vertical]
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pplib  # noqa: E402

FIELDS = ["step", "t_s", "body", "area_m2", "area_core_m2", "com_x", "com_y",
          "xmin", "xmax", "ymin", "ymax", "face_lo", "face_hi", "n_contacts"]
RATE_COLS = ["com_x", "com_y", "face_lo", "face_hi", "area_m2"]


# ---------------------------------------------------------------------------
# One snapshot
# ---------------------------------------------------------------------------
def label_cores(phi, level=0.5):
    """Connected components of phi > level (4-connected). Returns (labels, n)."""
    from scipy import ndimage
    lab, n = ndimage.label(phi > level)
    return lab, n


def assign_nodes(lab):
    """Give every node to the nearest core (index-space distance)."""
    from scipy import ndimage
    if lab.max() == 0:
        return lab.copy()
    _, (jj, ii) = ndimage.distance_transform_edt(lab == 0, return_indices=True)
    return lab[jj, ii]


def axis_line(X, Y, y0):
    """(x, j, w): for each column i, the x and the interpolation weights that
    put a field on the line y = y0, column by column (Y monotone in j).
    Columns that do not span y0 get NaN."""
    ny, nx = Y.shape
    j = np.empty(nx, dtype=np.intp)
    w = np.full(nx, np.nan)
    x = np.full(nx, np.nan)
    for i in range(nx):
        col = Y[:, i]
        if not (col.min() <= y0 <= col.max()):
            j[i] = 0
            continue
        k = int(np.clip(np.searchsorted(col, y0) - 1, 0, ny - 2))
        dy = col[k + 1] - col[k]
        j[i] = k
        w[i] = 0.0 if dy == 0 else (y0 - col[k]) / dy
        x[i] = (1 - w[i]) * X[k, i] + w[i] * X[k + 1, i]
    return x, j, w


def on_line(F, j, w):
    i = np.arange(F.shape[1])
    jj = np.clip(j, 0, F.shape[0] - 2)
    return (1 - w) * F[jj, i] + w * F[jj + 1, i]


def wall_curves(X, Y, walls):
    """[(name, x, y), ...] node coordinates along each wall face, as index
    slices of the grid (axis 0 = x faces are columns, axis 1 = y faces rows)."""
    out = []
    for axis, side in walls:
        name = "xyz"[axis] + str(side)
        if axis == 0:
            sl = (slice(None), 0 if side == 0 else -1)
        elif axis == 1:
            sl = (0 if side == 0 else -1, slice(None))
        else:
            continue
        out.append((name, sl))
    return out


def measure(snap, eps, y_axis, walls, min_area):
    """Per-body measurements of one Snapshot. Returns (rows, contacts, core_labels)."""
    X, Y, phi = snap.X, snap.Y, snap.phi
    w = pplib.node_weights(X, Y)
    lab, n = label_cores(phi)
    # Drop cores below the area threshold before assigning the diffuse band.
    keep = np.zeros(n + 1, dtype=bool)
    if n:
        core_area = np.bincount(lab.ravel(), weights=w.ravel(), minlength=n + 1)
        keep[1:] = core_area[1:] >= min_area
    relabel = np.zeros(n + 1, dtype=np.intp)
    relabel[np.flatnonzero(keep)] = np.arange(1, keep.sum() + 1)
    lab = relabel[lab]
    n = int(keep.sum())
    share = assign_nodes(lab)

    # Axis line and its crossings. Each crossing belongs to the body that owns
    # the node on its ice side (the higher-phi neighbor) on the nearer row.
    xl, jl, wl = axis_line(X, Y, y_axis)
    ok = np.flatnonzero(np.isfinite(xl))
    xs = xl[ok]
    ps = on_line(phi, jl, np.nan_to_num(wl))[ok]
    rows_near = np.clip(jl[ok] + (np.nan_to_num(wl[ok]) >= 0.5), 0, phi.shape[0] - 1)
    own = share[rows_near, ok]
    cross = pplib.crossings(xs, ps)
    owner = []
    for c in cross:
        k = int(np.clip(np.searchsorted(xs, c), 1, xs.size - 1))
        owner.append(int(own[k] if ps[k] > ps[k - 1] else own[k - 1]))
    if eps and np.isfinite(eps):
        cross = np.array([pplib.refine_tanh(xs, ps, c, eps) for c in cross])

    rows = []
    for b in range(1, n + 1):
        m = share == b
        core = lab == b
        wb = w * phi * m
        A = wb.sum()
        mine = sorted(c for c, o in zip(cross, owner) if o == b)
        face_lo = face_hi = np.nan
        if len(mine) >= 2:
            face_lo, face_hi = mine[0], mine[-1]
        elif len(mine) == 1:
            # One crossing: the body runs into a domain end on the axis. The
            # crossing is the face on the side where phi falls toward 0.
            c = mine[0]
            k = int(np.clip(np.searchsorted(xs, c), 1, xs.size - 1))
            if ps[k] < ps[k - 1]:
                face_hi = c
            else:
                face_lo = c
        rows.append({
            "body": b, "area_m2": A, "area_core_m2": float((w * core).sum()),
            "com_x": (wb * X).sum() / A if A > 0 else np.nan,
            "com_y": (wb * Y).sum() / A if A > 0 else np.nan,
            "xmin": float(X[core].min()), "xmax": float(X[core].max()),
            "ymin": float(Y[core].min()), "ymax": float(Y[core].max()),
            "face_lo": face_lo, "face_hi": face_hi, "n_contacts": 0,
            "n_axis_crossings": len(mine),
        })

    contacts = []
    for name, sl in wall_curves(X, Y, walls):
        cx, cy, cp, cs = X[sl], Y[sl], phi[sl], share[sl]
        s = np.r_[0.0, np.cumsum(np.hypot(np.diff(cx), np.diff(cy)))]   # arc length
        for c in pplib.crossings(s, cp):
            k = int(np.clip(np.searchsorted(s, c), 1, s.size - 1))
            b = int(cs[k] if cp[k] > cp[k - 1] else cs[k - 1])   # the ice side
            f = (c - s[k - 1]) / (s[k] - s[k - 1])
            contacts.append({"step": snap.step, "t_s": snap.t, "body": b, "wall": name,
                             "s": c, "x": cx[k - 1] + f * (cx[k] - cx[k - 1]),
                             "y": cy[k - 1] + f * (cy[k] - cy[k - 1])})
            if 1 <= b <= n:
                rows[b - 1]["n_contacts"] += 1
    return rows, contacts, lab


# ---------------------------------------------------------------------------
# Tracking
# ---------------------------------------------------------------------------
def match(prev_lab, prev_ids, lab, next_id):
    """Map the bodies of `lab` to persistent ids. Returns (ids, next_id, events)."""
    n = int(lab.max())
    ids = [None] * (n + 1)
    events = []
    if prev_lab is None:
        for b in range(1, n + 1):
            ids[b] = next_id
            next_id += 1
        return ids, next_id, events
    m = int(prev_lab.max())
    ov = np.zeros((m + 1, n + 1))
    both = (prev_lab > 0) & (lab > 0)
    np.add.at(ov, (prev_lab[both], lab[both]), 1)
    parents = {b: [a for a in range(1, m + 1) if ov[a, b] > 0] for b in range(1, n + 1)}
    children = {a: [b for b in range(1, n + 1) if ov[a, b] > 0] for a in range(1, m + 1)}
    for b in range(1, n + 1):
        ps = parents[b]
        if not ps:
            ids[b] = next_id
            events.append(("birth", next_id, []))
            next_id += 1
        elif len(ps) > 1:
            ids[b] = next_id
            events.append(("merge", next_id, [prev_ids[a] for a in ps]))
            next_id += 1
        else:
            a = ps[0]
            if len(children[a]) == 1:
                ids[b] = prev_ids[a]
            else:
                # Split: the child holding most of the parent keeps the id.
                best = max(children[a], key=lambda c: ov[a, c])
                if b == best:
                    ids[b] = prev_ids[a]
                else:
                    ids[b] = next_id
                    events.append(("split", next_id, [prev_ids[a]]))
                    next_id += 1
    for a in range(1, m + 1):
        if not children[a]:
            events.append(("death", prev_ids[a], []))
    return ids, next_id, events


# ---------------------------------------------------------------------------
# Rates and fixed-value sampling
# ---------------------------------------------------------------------------
def add_rates(rows):
    """Add d<col>/dt and d2<col>/dt2 for each body's time series in place."""
    by = {}
    for r in rows:
        by.setdefault(r["body"], []).append(r)
    for series in by.values():
        series.sort(key=lambda r: r["t_s"])
        t = np.array([r["t_s"] for r in series])
        for col in RATE_COLS:
            y = np.array([r[col] for r in series], dtype=float)
            v = np.full(t.size, np.nan)
            a = np.full(t.size, np.nan)
            good = np.isfinite(y)
            # Rates only over unbroken runs of finite values, never across a gap.
            start = 0
            while start < t.size:
                if not good[start]:
                    start += 1
                    continue
                end = start
                while end + 1 < t.size and good[end + 1]:
                    end += 1
                seg = slice(start, end + 1)
                if end - start + 1 >= 2:
                    edge = 2 if end - start + 1 >= 3 else 1
                    v[seg] = np.gradient(y[seg], t[seg], edge_order=edge)
                if end - start + 1 >= 3:
                    a[seg] = np.gradient(v[seg], t[seg], edge_order=2)
                start = end + 1
            for r, vi, ai in zip(series, v, a):
                r["v_" + col] = vi
                r["a_" + col] = ai
    return rows


def add_gap(rows, pair):
    """Gap on the axis between two bodies: face_lo of the deeper minus
    face_hi of the shallower, per snapshot. Written into both bodies' rows."""
    by_t = {}
    for r in rows:
        by_t.setdefault(r["t_s"], {})[r["body"]] = r
    for bodies in by_t.values():
        g = np.nan
        if pair[0] in bodies and pair[1] in bodies:
            a, b = bodies[pair[0]], bodies[pair[1]]
            if a["com_x"] > b["com_x"]:
                a, b = b, a
            g = b["face_lo"] - a["face_hi"]
        for r in bodies.values():
            r["gap"] = g
    return rows


def sample_at(rows, col, value):
    """Rows interpolated in time where `col` of a body crosses `value`."""
    out = []
    by = {}
    for r in rows:
        by.setdefault(r["body"], []).append(r)
    for body, series in by.items():
        series.sort(key=lambda r: r["t_s"])
        y = np.array([r.get(col, np.nan) for r in series], dtype=float)
        for k in range(len(series) - 1):
            y0, y1 = y[k], y[k + 1]
            if not (np.isfinite(y0) and np.isfinite(y1)) or y0 == y1:
                continue
            if (y0 - value) * (y1 - value) > 0:
                continue
            if y1 == value and k + 2 < len(series):
                continue                    # counted as the next interval's start
            f = (value - y0) / (y1 - y0)
            r0, r1 = series[k], series[k + 1]
            row = {"body": body, "crossing_of": col, "value": value}
            for key in r0:
                if key in ("body", "step"):
                    continue
                a0, a1 = r0.get(key, np.nan), r1.get(key, np.nan)
                try:
                    row[key] = (1 - f) * float(a0) + f * float(a1)
                except (TypeError, ValueError):
                    row[key] = a0
            row["step_before"], row["step_after"] = r0["step"], r1["step"]
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def analyze(run_dirs, source="auto", n_per_elem=2, axis_y=None, min_area=None, every=1):
    opts = pplib.read_opts(run_dirs[0])
    eps = pplib.opt_float(opts, "-eps", np.nan)
    walls = pplib.face_list(opts, "-wall_faces")
    rows, contacts, events = [], [], []
    prev_lab, prev_ids, next_id = None, None, 0
    for snap in pplib.iter_snapshots(run_dirs, source=source, n_per_elem=n_per_elem, every=every):
        y0 = axis_y if axis_y is not None else 0.5 * (snap.Y.min() + snap.Y.max())
        thr = min_area if min_area is not None else (np.pi * (2 * eps) ** 2 if np.isfinite(eps) else 0.0)
        r, c, lab = measure(snap, eps, y0, walls, thr)
        ids, next_id, ev = match(prev_lab, prev_ids, lab, next_id)
        for row in r:
            row.update(step=snap.step, t_s=snap.t, body=ids[row["body"]])
        for row in c:
            if 1 <= row["body"] < len(ids):
                row["body"] = ids[row["body"]]
        for kind, bid, parents in ev:
            events.append((snap.step, snap.t, kind, bid, parents))
        rows += r
        contacts += c
        prev_lab, prev_ids = lab, ids
    return rows, contacts, events, opts


def write_csv(path, rows, first=None):
    if not rows:
        with open(path, "w") as fh:
            fh.write("# no rows\n")
        return
    keys = list(first or [])
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        wtr = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        wtr.writeheader()
        for r in rows:
            wtr.writerow({k: (f"{v:.10e}" if isinstance(v, float) else v) for k, v in r.items()})


def plot(rows, path, vertical):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from figcheck import check_figure
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 3.0), constrained_layout=True)
    tmax = max(r["t_s"] for r in rows)
    unit = pplib.auto_time_unit(tmax)
    for b in sorted({r["body"] for r in rows}):
        s = sorted((r for r in rows if r["body"] == b), key=lambda r: r["t_s"])
        t = pplib.in_time_unit([r["t_s"] for r in s], unit)
        a1.plot(t, [r["area_m2"] * 1e12 for r in s], lw=1.5, label=f"body {b}")
        for col, ls in (("face_lo", "-"), ("face_hi", "--")):
            a2.plot(t, [r[col] * 1e6 for r in s], ls=ls, lw=1.2, color=a1.lines[-1].get_color())
    a1.set_xlabel(f"time [{unit}]")
    a1.set_ylabel("ice area [µm²]")
    a2.set_xlabel(f"time [{unit}]")
    a2.set_ylabel("face position on the axis x [µm]")
    if vertical:
        a2.invert_yaxis()
        a2.set_ylabel("face depth x [µm]")
    a1.legend(frameon=False, fontsize=8)
    for ax in (a1, a2):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.savefig(path, dpi=200)
    problems = check_figure(fig, min_pt=7.0)
    plt.close(fig)
    return problems


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", action="append", required=True, help="run folder (repeat for a restart chain)")
    p.add_argument("--source", choices=("auto", "nurbs", "vts"), default="auto")
    p.add_argument("--n-per-elem", type=int, default=2)
    p.add_argument("--every", type=int, default=1, help="use every n-th snapshot")
    p.add_argument("--axis-y", type=float, default=None, help="pore axis line y [m] (default mid-height)")
    p.add_argument("--min-area", type=float, default=None,
                   help="smallest core that counts as a body [m^2] (default pi (2 eps)^2)")
    p.add_argument("--at", action="append", default=[], metavar="COL=VALUE")
    p.add_argument("--gap-bodies", default="0,1", help="ids of the two bodies whose axis gap is 'gap'")
    p.add_argument("--out", default=None, help="output folder (default: the first --dir)")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--vertical", action="store_true", help="plot depth downward (plan tool T7)")
    a = p.parse_args(argv)

    rows, contacts, events, opts = analyze(a.dir, a.source, a.n_per_elem, a.axis_y, a.min_area, a.every)
    if not rows:
        print("ice_bodies: no ice bodies in any snapshot", file=sys.stderr)
    add_rates(rows)
    pair = tuple(int(v) for v in a.gap_bodies.split(","))
    add_gap(rows, pair)
    out = a.out or a.dir[0]
    os.makedirs(out, exist_ok=True)
    write_csv(os.path.join(out, "ice_bodies.csv"), rows, FIELDS)
    write_csv(os.path.join(out, "ice_bodies_contacts.csv"), contacts)
    with open(os.path.join(out, "ice_bodies_events.txt"), "w") as fh:
        fh.write("# step t_s event body parents\n")
        for st, t, kind, bid, parents in events:
            fh.write(f"{st} {t:.6e} {kind} {bid} {','.join(map(str, parents)) or '-'}\n")
    at_rows = []
    for spec in a.at:
        col, _, val = spec.partition("=")
        at_rows += sample_at(rows, col.strip(), float(val))
    if a.at:
        write_csv(os.path.join(out, "ice_bodies_at.csv"), at_rows, ["body", "crossing_of", "value", "t_s"])
    n_snap = len({r["step"] for r in rows})
    print(f"ice_bodies: {n_snap} snapshots, {len({r['body'] for r in rows})} bodies, "
          f"{len(events)} events, eps {pplib.opt_float(opts, '-eps', float('nan')):.4e} m, "
          f"walls {opts.get('-wall_faces', '-')}; wrote {out}/ice_bodies*.csv")
    if a.plot and rows:
        os.makedirs(os.path.join(out, "plots"), exist_ok=True)
        probs = plot(rows, os.path.join(out, "plots", "ice_bodies.png"), a.vertical)
        for pr in probs:
            print("  figure check:", pr, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
