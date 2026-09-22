#!/usr/bin/env python3
"""Red-flag scan: is this pilot trustworthy, and is the domain still an REV?

    venv_enceladus/bin/python studies/keff_sintering/health_check.py <batch_dir>

Two separate questions, deliberately kept apart:

NUMERICAL HEALTH -- did the solver do what it claims?
  mass conservation, phase-field bounds, k_eff linear-solve convergence,
  off-diagonal symmetry of the tensor, rejected steps, guard trips.

REPRESENTATIVENESS -- is a 2 mm box still an REV at the END of the run?
  Coarsening grows the grains, so L/R_ave FALLS during the run. A domain
  sized as an REV at t=0 need not be one at t_final, and the operational
  test (Kanit et al.) is whether independent realisations agree: if the
  seed-to-seed scatter in k_eff GROWS with time, the box is losing its
  representativeness as the structure coarsens into it.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "postprocess"))
import pplib                                                # noqa: E402

DAY = 86400.0
OK, BAD = "  ok  ", " FLAG "


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
        # a merged run carries MERGE_INFO.json; prefer it over a raw leg
        if prev is None or ((d / "MERGE_INFO.json").is_file()
                            and not (prev / "MERGE_INFO.json").is_file()):
            seen[key] = d
    return dict(sorted(seen.items()))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path)
    ap.add_argument("--L", type=float, default=2.0e-3)
    ap.add_argument("--LR0", type=float, default=40.0, help="L/R_ave at t=0")
    a = ap.parse_args()

    runs = {}
    for seed, d in find_runs(a.batch).items():
        ssa = pplib.load_ssa(str(d))
        k = np.atleast_1d(np.genfromtxt(d / "k_eff.csv", delimiter=",", names=True))
        runs[seed] = {"ssa_raw": ssa, "k": k, "dir": d,
                      "ssa": 6.0 * ssa[:, 0] / a.L**2, "t": ssa[:, 2],
                      "ice": ssa[:, 1], "mass": ssa[:, 7], "dt": ssa[:, 4]}
    if not runs:
        raise SystemExit("no merged runs found")

    print("=" * 72)
    print("NUMERICAL HEALTH")
    print("=" * 72)
    for s, r in runs.items():
        k = r["k"]
        drift = abs(r["mass"][-1] / r["mass"][0] - 1)
        ice_d = abs(r["ice"][-1] / r["ice"][0] - 1)
        offsym = np.max(np.abs(k["k_01"] - k["k_10"]) / np.abs(k["k_00"]))
        offmag = np.max(np.abs(k["k_01"]) / k["k_00"])
        bad_ksp = int(np.sum(k["ksp_reason"] <= 0))
        its = k["ksp_its"]
        print(f"\n  seed {s}")
        # 1e-6, not 1e-9. The pilot happened to conserve to exactly 0.0, and
        # a threshold set from that flags ordinary floating-point accumulation:
        # the L/R_ave = 64 runs drift 1.5e-7 relative over 371 steps, which is
        # 0.00002% and is not a fault. 1e-6 still catches a real leak.
        print(f"   {OK if drift   < 1e-6 else BAD} mass drift            {drift:.3e}")
        print(f"   {OK if ice_d   < 1e-6 else BAD} ice-fraction drift    {ice_d:.3e}")
        print(f"   {OK if bad_ksp == 0   else BAD} k_eff solves failed   {bad_ksp} of {len(k)}")
        print(f"   {OK if offsym  < 1e-6 else BAD} tensor symmetry       "
              f"max|k01-k10|/k00 = {offsym:.2e}")
        print(f"        off-diagonal magnitude max|k01|/k00 = {offmag:.3f}"
              f"   (structure, not error)")
        print(f"        k_eff KSP iterations {int(its.min())}..{int(its.max())}")
        # warnings the solver itself raised
        outp = r["dir"] / "outp.txt"
        if outp.is_file():
            txt = outp.read_text(errors="ignore")
            cfl = len(re.findall(r"Interface-CFL violated", txt))
            guard = len(re.findall(r"PHASE GUARD TRIPPED", txt))
            oob = len(re.findall(r"out of bounds", txt))
            # A FEW CFL rollbacks are the system working, not failing: they
            # mean -dtmax was set slightly above what the interface allows and
            # the limiter caught it, which is exactly the division of labour
            # (dtmax proposes, the measured phase rate disposes). Flag only if
            # they are frequent enough to be wasting real time, or if the
            # phase guard -- which has no such excuse -- ever tripped.
            nstep = len(r["t"])
            frac = cfl / max(nstep, 1)
            print(f"   {OK if frac < 0.05 else BAD} CFL rollbacks         "
                  f"{cfl} in {nstep} steps ({frac:.1%}) — the limiter working")
            print(f"   {OK if guard+oob == 0 else BAD} phase guard / bounds  "
                  f"{guard} trips, {oob} out-of-bounds")

    print("\n" + "=" * 72)
    print("REPRESENTATIVENESS  (is a 2 mm box still an REV at t_final?)")
    print("=" * 72)

    # Grain size from SSA. For a fixed morphology SSA ~ 1/R, so the coarsening
    # factor is SSA(0)/SSA(t) and L/R_ave falls by the same factor. This is a
    # proxy -- it assumes the shape does not change, which sintering violates --
    # so it is a LOWER bound on the coarsening and an UPPER bound on L/R_ave.
    print("\n  coarsening (from SSA, assuming SSA ~ 1/R):")
    print(f"  {'seed':>5} {'SSA0':>8} {'SSAend':>8} {'g=R/R0':>8} {'L/R_ave end':>12}")
    for s, r in runs.items():
        g = r["ssa"][0] / r["ssa"][-1]
        print(f"  {s:>5} {r['ssa'][0]:8.0f} {r['ssa'][-1]:8.0f} {g:8.3f} "
              f"{a.LR0/g:12.1f}")

    # The operational REV test: scatter across independent realisations.
    print("\n  seed-to-seed scatter in k_eff vs time"
          "  (the operational REV criterion):")
    tgrid = np.linspace(0, min(r["k"]["time"].max() for r in runs.values()), 25)
    print(f"  {'t[d]':>7} {'mean k_eff':>11} {'sd':>9} {'CV':>7} {'L/R_ave':>8}")
    rows = []
    for tt in tgrid[[0, 2, 6, 12, 18, 24]]:
        vals = np.array([np.interp(tt, r["k"]["time"], r["k"]["k_iso"])
                         for r in runs.values()])
        gs = np.mean([r["ssa"][0] / np.interp(tt, r["t"], r["ssa"])
                      for r in runs.values()])
        cv = vals.std(ddof=1) / vals.mean()
        rows.append((tt / DAY, vals.mean(), vals.std(ddof=1), cv, a.LR0 / gs))
        print(f"  {tt/DAY:7.2f} {vals.mean():11.4f} {vals.std(ddof=1):9.4f} "
              f"{cv:7.1%} {a.LR0/gs:8.1f}")

    cv0, cv1 = rows[0][3], rows[-1][3]
    print(f"\n  CV grows {cv0:.1%} -> {cv1:.1%} over the run.")
    if cv1 > 1.5 * cv0:
        print("   FLAG  The realisations are DIVERGING. A box that averages well")
        print("         at t=0 averages worse once the structure has coarsened")
        print("         into it -- exactly the L/R_ave drop above. Size the")
        print("         production domain for t_final, not t=0.")
    print(f"\n  with 4 seeds the standard error of the mean is "
          f"{cv1/np.sqrt(len(runs)):.1%} at t_final;")
    print(f"  a claimed difference between conditions must clear that.")

    # early-time behaviour: is the first day physics or IC relaxation?
    print("\n  early transient (is day 1 sintering, or the IC relaxing?):")
    for s, r in list(runs.items())[:1]:
        m = (r["t"] > 2 * DAY)
        p = np.polyfit(np.log(r["t"][m]), np.log(r["ssa"][m]), 1)
        pred1 = np.exp(np.polyval(p, np.log(1 * DAY)))
        act1 = np.interp(1 * DAY, r["t"], r["ssa"])
        print(f"    seed {s}: SSA ~ t^{p[0]:+.4f} fitted on t > 2 d")
        print(f"      at t = 1 d that law predicts {pred1:.0f}, actual {act1:.0f}"
              f"  ({act1/pred1-1:+.1%})")
        print("      a large mismatch means the first day is the initial")
        print("      condition equilibrating, not sintering, and should be")
        print("      excluded from fits and from the t=0 reference.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
