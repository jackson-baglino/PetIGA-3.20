#!/usr/bin/env python3
"""Why do packings fan out in k_eff against SSA? Two checks on the exponent.

    venv_enceladus/bin/python studies/keff_sintering/master_curve/exponent_scatter.py <campaign dir>

p is each packing's exponent in k_eff/k_r = (SSA/SSA_r)^p, fitted from the
reference age theta_r = 30 on (as the state-law figure).

1. DOMAIN SIZE. phi 0.325, -20 C, ungated packings at L/R = 20, 30, 40, 56,
   80: does the packing-to-packing scatter of p, and of the level k_r, shrink
   as the cell grows? Sampling of a finite cell predicts scatter ~ (L/R)^-1.
2. CONNECTIVITY. For the 25 production packings (-10 C runs), p against
     eta    = k_r / k_HS, the conductivity at theta_r over the 2D
              Hashin-Shtrikman upper bound for that ice fraction: how much of
              the best possible conduction the ice network delivers;
     aniso  = min(k_xx, k_yy) / max(k_xx, k_yy) at theta_r;
     z_band = contacts per grain of the initial packing (master_packings.csv).
   Correlations are given raw and with each porosity's mean removed, i.e.
   between packings at fixed porosity.
Writes exponent_scatter.txt next to this file.
"""
import collections, re, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402

K_I, K_A, TH_R = 2.29, 0.02, 30.0


def fit(d):
    r = load(d)
    if r is None:
        return None
    tau = read_tau_sub(d); th = r["t"] / tau
    if th[-1] < 1.5 * TH_R:
        return None
    kr, sr = np.interp(TH_R, th, r["kiso"]), np.interp(TH_R, th, r["ssa"])
    w = th >= TH_R
    X, Y = np.log10(r["ssa"][w] / sr), np.log10(r["kiso"][w] / kr)
    c = np.genfromtxt(d / "k_eff.csv", delimiter=",", names=True); thc = c["time"] / tau
    f = float(np.interp(TH_R, thc, c["phi_bar"]))                 # ice fraction
    hs = K_I + (1 - f) * (K_A - K_I) * 2 * K_I / (2 * K_I + f * (K_A - K_I))   # 2D HS upper bound
    kxx, kyy = np.interp(TH_R, thc, c["k_00"]), np.interp(TH_R, thc, c["k_11"])
    return dict(p=float(np.sum(X * Y) / np.sum(X * X)), kr=float(kr), eta=float(kr / hs),
                aniso=float(min(kxx, kyy) / max(kxx, kyy)))


def main():
    L = Path(sys.argv[1]); out = []
    P = lambda s="": (print(s), out.append(s))
    P("1. Domain size (phi 0.325, -20 C, ungated packings)")
    by = collections.defaultdict(list)
    for d in sorted(L.glob("packing_2D_phi0.325_Rave50um_LR*_eps1000nm_perxy_T-20__*")):
        m = re.search(r"LR(\d+)_seed(\d+)", d.name)
        if int(m.group(2)) < 1000:
            continue
        e = fit(d)
        if e:
            by[int(m.group(1))].append(e)
    P("   L/R   n   exponent mean    sd     level k_r/k_ice   sd/mean")
    for lr in sorted(by):
        p = np.array([e["p"] for e in by[lr]]); k = np.array([e["kr"] for e in by[lr]]) / K_I
        P(f"   {lr:3d}  {len(p):2d}     {p.mean():+.3f}      {p.std(ddof=1):.3f}      {k.mean():.3f}         {100 * k.std(ddof=1) / k.mean():4.1f} %")
    lr = np.array(sorted(by))
    for lab, v in (("exponent", [np.std([e["p"] for e in by[l]], ddof=1) for l in lr]),
                   ("level   ", [np.std([e["kr"] for e in by[l]], ddof=1) / np.mean([e["kr"] for e in by[l]]) for l in lr])):
        P(f"   scatter of the {lab} ~ (L/R)^{np.polyfit(np.log(lr), np.log(v), 1)[0]:+.2f}")
    P()
    P("2. Connectivity (25 production packings, -10 C runs)")
    z = {}
    mp = np.genfromtxt(HERE / "master_packings.csv", delimiter=",", names=True)
    for s, zb in zip(mp["seed"], mp["z_band"]):
        z[int(s)] = zb
    rows = []
    for d in sorted(L.glob("packing_2D_phi0.*_Rave50um_LR40_seed[12][0-9][0-9][0-9]_L2mm_eps1000nm_perxy_T-10__snow_T-10_h1.00_30d")):
        m = re.search(r"phi([\d.]+)_.*seed(\d+)", d.name); e = fit(d)
        rows.append((float(m.group(1)), e["p"], e["eta"], e["aniso"], z.get(int(m.group(2)), np.nan)))
    A = np.array(rows)
    P("   phi    exponent        eta = k_r/k_HS    min/max k     z_band")
    for phi in sorted(set(A[:, 0])):
        s = A[A[:, 0] == phi]
        P(f"   {phi:.3f}  {s[:, 1].mean():+.2f} ± {s[:, 1].std(ddof=1):.2f}    {s[:, 2].mean():.2f} ± {s[:, 2].std(ddof=1):.2f}      {s[:, 3].mean():.2f} ± {s[:, 3].std(ddof=1):.2f}   {s[:, 4].mean():.2f}")
    cc = lambda a, b: np.corrcoef(a, b)[0, 1]
    for name, sel in (("all 25 packings   ", np.ones(len(A), bool)), ("phi <= 0.375 (15) ", A[:, 0] <= 0.376)):
        S = A[sel]; Rm = S.copy()
        for phi in set(S[:, 0]):
            k = S[:, 0] == phi; Rm[k, 1:] = S[k, 1:] - S[k, 1:].mean(0)
        P(f"   {name} corr(exponent, .) raw / at fixed porosity:  eta {cc(S[:, 1], S[:, 2]):+.2f} / {cc(Rm[:, 1], Rm[:, 2]):+.2f}"
          f"   min/max k {cc(S[:, 1], S[:, 3]):+.2f} / {cc(Rm[:, 1], Rm[:, 3]):+.2f}   z_band {cc(S[:, 1], S[:, 4]):+.2f} / {cc(Rm[:, 1], Rm[:, 4]):+.2f}")
        X = np.column_stack([np.ones(len(S)), S[:, 2], S[:, 0]]); c, *_ = np.linalg.lstsq(X, S[:, 1], rcond=None)
        P(f"      exponent = {c[0]:+.2f} {c[1]:+.2f} eta {c[2]:+.2f} phi;  rms residual {np.sqrt(np.mean((S[:, 1] - X @ c) ** 2)):.3f}"
          f"  (sd of the exponent {S[:, 1].std(ddof=1):.3f}; about porosity means {Rm[:, 1].std(ddof=1):.3f})")
    (HERE / "exponent_scatter.txt").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
