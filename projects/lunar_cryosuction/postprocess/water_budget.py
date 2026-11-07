#!/usr/bin/env python3
"""Water budget of a run: where the ice went (simulation plan tool T2;
measurements M2 and M3).

THE BALANCE

All water in the domain is ice or vapor. Per unit depth in 2D (kg/m):

    ice mass      M_ice(t) = rho_ice * tot_ice(t)
    vapor stored  M_vap(t) = tot_rhov(t) = integral of rho_v (1 - phi)
    outflow       O(t)     = integral from t0 to t of sum over pinned faces of J_vap
                             (bflux.dat; positive out of the domain)

and conservation says

    M_ice(t) + M_vap(t) + O(t) = M_ice(t0) + M_vap(t0).

The closure error is the left side minus the right side, reported relative to
the ice lost, rho_ice (tot_ice(t0) - tot_ice(t)). The solver integrates
tot_ice and tot_rhov with its own quadrature every step (SSA_evo.dat); bflux.dat
has the face fluxes every accepted step. O(t) is the trapezoid rule over the
bflux rows, so its error is second order in the step; a run whose flux changes
fast between rows shows it as closure error. A closed box has no bflux.dat
(no face pinned); then O = 0 and the check is that M_ice + M_vap stays flat.

RECAPTURE (M3)

Of the ice that sublimates, how much is redeposited inside the pore instead
of leaving? From consecutive snapshots k, k+1 the local ice change is
dphi = phi_{k+1} - phi_k. Summed over the run,

    sublimated = rho_ice * sum_k integral of max(-dphi, 0)
    deposited  = rho_ice * sum_k integral of max(+dphi, 0)
    recapture  = deposited / sublimated

A lens translating without loss gives recapture = 1; a body losing ice with
no redeposition gives 0. Changes that reverse between two snapshots cancel,
so this is a LOWER bound on both sums that tightens as snapshots get closer;
the script reports it at full cadence and at half cadence so the size of that
bias is visible (if the two differ by more than a few percent, save snapshots
more often).

HEAT

Q_in(t) = -integral of sum over faces of Q_heat (bflux.dat; positive into the
domain). Reported per kg of water that left, next to the latent heat of
sublimation 2.83e6 J/kg (for reference only; the model's latent heat is in
the solver opts).

OUTPUT

    water_budget.csv    per SSA_evo.dat row: t, ice, vapor, outflow, closure
    water_budget.txt    the summary printed below

Usage:
    python3 postprocess/water_budget.py --dir <run> [--dir <restart leg> ...]
        [--no-recapture] [--source auto|nurbs|vts]
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pplib  # noqa: E402

L_SUB_REF = 2.83e6   # J/kg near 250 K (reference value only)


def load_bflux(run_dir):
    """(t, J_total, Q_total, {face: J}) from bflux.dat, or None when absent."""
    path = os.path.join(run_dir, "bflux.dat")
    if not os.path.isfile(path):
        return None
    names = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("# t step"):
                names = line[1:].split()[2:]
                break
    a = np.genfromtxt(path, comments="#", invalid_raise=False, ndmin=2)
    if a.size == 0:
        return None
    a = a[np.isfinite(a).all(axis=1)]
    # A restart repeats steps; keep the last row of each step, in time order.
    _, last = np.unique(a[::-1, 1], return_index=True)
    a = a[np.sort(len(a) - 1 - last)]
    a = a[np.argsort(a[:, 0], kind="stable")]
    jcols = [k for k, n in enumerate(names) if n.startswith("J_vap")]
    qcols = [k for k, n in enumerate(names) if n.startswith("Q_heat")]
    t = a[:, 0]
    J = a[:, [2 + k for k in jcols]].sum(axis=1) if jcols else np.zeros_like(t)
    Q = a[:, [2 + k for k in qcols]].sum(axis=1) if qcols else np.zeros_like(t)
    faces = {names[k][len("J_vap_"):]: a[:, 2 + k] for k in jcols}
    return t, J, Q, faces


def chain(run_dirs, loader):
    """Concatenate a per-folder time series over a restart chain; a later leg
    wins where times overlap."""
    parts = [p for p in (loader(d) for d in run_dirs) if p is not None]
    if not parts:
        return None
    out = parts[0]
    for nxt in parts[1:]:
        t0 = nxt[0][0]
        keep = out[0] < t0
        out = tuple(np.r_[o[keep], n] if isinstance(o, np.ndarray) else
                    {k: np.r_[o[k][keep], n[k]] for k in o} for o, n in zip(out, nxt))
    return out


def load_ssa_cols(run_dir):
    a = pplib.load_ssa(run_dir)
    if a is None:
        return None
    if a.shape[1] <= pplib.RHOV:
        raise SystemExit(f"{run_dir}: SSA_evo.dat has no tot_rhov column (run predates it); "
                         "the vapor term of the budget cannot be formed")
    order = np.argsort(a[:, pplib.TIME], kind="stable")
    a = a[order]
    return a[:, pplib.TIME], a[:, pplib.ICE], a[:, pplib.RHOV]


def cumulative(t, y):
    """Cumulative trapezoid integral, 0 at t[0]."""
    return np.r_[0.0, np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))]


def recapture(run_dirs, source, every):
    """(sublimated, deposited) in kg/m over the snapshot sequence."""
    prev = None
    w = None
    sub = dep = 0.0
    n = 0
    for s in pplib.iter_snapshots(run_dirs, source=source, every=every):
        if w is None:
            w = pplib.node_weights(s.X, s.Y)
        phi = np.clip(s.phi, 0.0, 1.0)
        if prev is not None:
            d = phi - prev
            sub += float((w * np.clip(-d, 0, None)).sum())
            dep += float((w * np.clip(d, 0, None)).sum())
        prev = phi
        n += 1
    return pplib.RHO_ICE * sub, pplib.RHO_ICE * dep, n


def budget(run_dirs, rho_ice=None):
    opts = pplib.read_opts(run_dirs[0])
    rho_ice = rho_ice or pplib.opt_float(opts, "-rho_ice", pplib.RHO_ICE)
    ssa = chain(run_dirs, load_ssa_cols)
    if ssa is None:
        raise SystemExit("no SSA_evo.dat in " + ", ".join(run_dirs))
    t, ice, vap = ssa
    bf = chain(run_dirs, load_bflux)
    if bf is None:
        out = np.zeros_like(t)
        heat = np.zeros_like(t)
        faces = {}
    else:
        tb, J, Q, faces = bf
        if tb.size < 2:
            out = np.zeros_like(t)
            heat = np.zeros_like(t)
        else:
            Oc = cumulative(tb, J)
            Hc = cumulative(tb, -Q)
            # Put the outflow on the SSA times. Before the first bflux row the
            # outflow is 0; after the last it is held (the run ended there).
            out = np.interp(t, tb, Oc, left=0.0, right=Oc[-1])
            heat = np.interp(t, tb, Hc, left=0.0, right=Hc[-1])
            if tb[0] > t[0] + 1e-12 * max(1.0, abs(t[-1])):
                print(f"water_budget: bflux.dat starts at t = {tb[0]:.4e} s, after SSA_evo.dat "
                      f"({t[0]:.4e} s); outflow before it is taken as 0", file=sys.stderr)
    m_ice = rho_ice * ice
    total0 = m_ice[0] + vap[0]
    closure = m_ice + vap + out - total0
    lost = m_ice[0] - m_ice
    return dict(t=t, ice=ice, m_ice=m_ice, vap=vap, out=out, heat=heat, closure=closure,
                lost=lost, faces=faces, rho_ice=rho_ice, has_bflux=bf is not None)


def summary(b, rec=None, rec_half=None):
    t, lost, out, vap, clo = b["t"], b["lost"], b["out"], b["vap"], b["closure"]
    ref = max(np.abs(lost).max(), 1e-300)
    lines = [
        f"rows {t.size}, t {t[0]:.4e} to {t[-1]:.4e} s ({(t[-1] - t[0]) / 86400:.3f} d), "
        f"rho_ice {b['rho_ice']:g} kg/m3, bflux.dat {'yes' if b['has_bflux'] else 'no (closed box: outflow 0)'}",
        f"ice at start {b['m_ice'][0]:.6e} kg/m, at end {b['m_ice'][-1]:.6e} kg/m "
        f"(fraction left {b['m_ice'][-1] / b['m_ice'][0] if b['m_ice'][0] else float('nan'):.5f})",
        f"ice lost to sublimation (net) {lost[-1]:.6e} kg/m",
        f"vapor held in the pore: start {vap[0]:.6e}, end {vap[-1]:.6e}, change {vap[-1] - vap[0]:+.6e} kg/m",
        f"vapor that left through the faces {out[-1]:.6e} kg/m",
        f"closure error (ice + vapor + outflow - initial): end {clo[-1]:+.3e} kg/m, "
        f"max |.| {np.abs(clo).max():.3e} kg/m = {np.abs(clo).max() / ref:.3e} of the ice lost",
    ]
    for f, J in b["faces"].items():
        lines.append(f"  face {f}: mean J {J.mean():+.4e} kg/(m s) over its rows")
    if out[-1] > 0:
        lines.append(f"heat drawn in through the faces {b['heat'][-1]:+.4e} J/m = "
                     f"{b['heat'][-1] / out[-1]:+.4e} J per kg that left (L_sub ref {L_SUB_REF:.2e})")
    if rec is not None:
        sub, dep, n = rec
        r = dep / sub if sub > 0 else float("nan")
        lines.append(f"recapture: sublimated {sub:.6e}, redeposited {dep:.6e} kg/m over {n} snapshots "
                     f"-> recapture fraction {r:.4f} (lower bound)")
        if rec_half is not None and rec_half[0] > 0:
            rh = rec_half[1] / rec_half[0]
            lines.append(f"  at half cadence ({rec_half[2]} snapshots): {rh:.4f}; "
                         f"cadence bias about {abs(r - rh):.4f}")
    return lines


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", action="append", required=True)
    p.add_argument("--no-recapture", action="store_true", help="skip the snapshot pass")
    p.add_argument("--source", choices=("auto", "nurbs", "vts"), default="auto")
    p.add_argument("--out", default=None)
    a = p.parse_args(argv)
    b = budget(a.dir)
    rec = rec_half = None
    if not a.no_recapture:
        try:
            rec = recapture(a.dir, a.source, 1)
            if rec[2] >= 5:
                rec_half = recapture(a.dir, a.source, 2)
        except FileNotFoundError as e:
            print(f"water_budget: recapture skipped ({e})", file=sys.stderr)
    lines = summary(b, rec, rec_half)
    out = a.out or a.dir[0]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "water_budget.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    hdr = "t_s,ice_m2,ice_kg_per_m,vapor_kg_per_m,outflow_kg_per_m,heat_in_J_per_m,closure_kg_per_m"
    np.savetxt(os.path.join(out, "water_budget.csv"),
               np.column_stack([b["t"], b["ice"], b["m_ice"], b["vap"], b["out"], b["heat"], b["closure"]]),
               delimiter=",", header=hdr, comments="", fmt="%.10e")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
