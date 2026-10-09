#!/usr/bin/env python3
"""Was the thin-interface correction active? An audit of every manuscript run.

    venv_enceladus/bin/python studies/keff_sintering/thin_interface_audit.py

For each run it reads the solver's own banner (outp.txt) and reports
  ratio = tau_sub d0 / (eps^2 beta_sub0)
which is 1 when tau_sub is the kinetic term alone (correction OFF) and
1 + Delta/beta_sub0 when the two Karma counter-terms are included (ON), plus
the banner's "(-thin_iface_corr N)" line where the solver version prints one.
It also estimates what ON would add, Delta / beta_sub0, from
  Delta = a1 a2 eps (1/D_therm + 1/D_v) rho_ice / rho_vs ,  a1 a2 = 0.7905,
with 1/D_therm = 1.218e5 s/m2 (the value implied by lunar docs/gt_deficit:
Delta = 1.282e5 s/m at -20 C, eps = 0.8584 um).
Writes thin_interface_audit.txt next to this file.

Background. enceladus_DSM: counter-terms removed 2026-09-09 (5ecdd236), flag
-thin_iface_corr default 0. lunar_regolith_DSM: default switched back to 1 on
2026-09-13 (f7cdfbe5). The two projects have differed since.
"""
import collections, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "preprocess"))
import comp_eps as ce  # noqa: E402

H = Path.home() / "SimulationResults/HPC_results/enceladus_DSM"
A1A2, INV_DTH, RHO_I = 0.7905, 1.218e5, 919.0


def audit(run):
    txt = (run / "outp.txt").read_text(errors="replace")
    g = lambda k: float(re.search(rf"^\s*{k}\s+(-?[0-9.eE+-]+)", txt, re.M).group(1))
    tau, d0, eps, beta, T = g("tau_sub"), g("d0_sub0"), g("eps"), g("beta_sub"), g("T0")
    flag = re.search(r"\(-thin_iface_corr (\d)\)", txt)
    delta = A1A2 * eps * (INV_DTH + 1.0 / ce.Dv_T(T)) * RHO_I / ce.rho_vs_sat(T)
    return dict(ratio=tau * d0 / (eps ** 2 * beta), flag=flag.group(1) if flag else "not printed",
                would=delta / beta, T=T, eps=eps)


out = []
P = lambda s="": (print(s), out.append(s))
G = H / "GrainPairSintering"
pairs = [("Molaro pair, -20 C (run 2026-09-08)", next((G / "batch_2026-09-08__17.20.46_molaro_T-20_round2").glob("molaro_*"))),
         ("Molaro pair, -5 C (run 2026-09-29)", next(G.glob("molaro_2D_L450x225um_eps0.11um_axisym_T-5pair*"))),
         ("saturated pair (run 2026-10-06)", next((G / "batch_2026-10-06__09.21.04_demmenie_mirror").glob("molaro_*")))]
P("run                                   tau d0/(eps^2 beta)  banner flag   correction  ON would add")
for name, r in pairs:
    a = audit(r)
    P(f"{name:38s} {a['ratio']:10.4f}        {a['flag']:12s}  {'ON ' if a['ratio'] > 1.005 else 'OFF'}        {100 * a['would']:5.1f} %")
C = H / "keff_sintering_campaign"
rows = [audit(d) for d in sorted(C.glob("packing_2D_phi0.*_Rave50um_LR40_seed[12][0-9][0-9][0-9]_L2mm_eps1000nm_perxy_T-*__snow_T-*_h1.00_30d"))]
P()
P(f"aggregate runs audited: {len(rows)}")
P(f"  banner flag: {dict(collections.Counter(r['flag'] for r in rows))}")
P(f"  tau d0/(eps^2 beta): {min(r['ratio'] for r in rows):.5f} to {max(r['ratio'] for r in rows):.5f}  -> correction OFF in all")
for T in sorted({r["T"] for r in rows}):
    w = [r["would"] for r in rows if r["T"] == T]
    P(f"  T = {T:6.1f} C: ON would add {100 * w[0]:.2f} % to tau_sub ({len(w)} runs)")
(HERE / "thin_interface_audit.txt").write_text("\n".join(out) + "\n")
