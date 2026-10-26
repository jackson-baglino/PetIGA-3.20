"""When does the eps-dependent early transient end?

Uses the paired interface-width runs (phi 0.325, -20 C, seeds 1701-1703,
eps = 1 um production and eps = 2 um). For a baseline time t_b, the rise
k(30 d)/k(t_b) - 1 is computed for both widths; the transient is over once the
two rises agree. The baseline is reported in days and in theta = t/tau_sub of
each width, read from each run's outp.txt.

    python relaxation_end.py <keff_sintering_campaign dir>
"""
import csv, re, sys
from pathlib import Path
import numpy as np

C = Path(sys.argv[1])
DAY = 86400.0

def load(d):
    t, k = [], []
    with open(d / "k_eff.csv") as f:
        for r in csv.DictReader(f):
            t.append(float(r["time"])); k.append(float(r["k_iso"]))
    tau = float(re.search(r"tau_sub\s+([0-9.eE+-]+)\s*s", (d / "outp.txt").read_text()).group(1))
    t, k = np.array(t), np.array(k)
    keep = t >= 1.0                       # drop the analytic IC frame
    return t[keep], k[keep], tau

def kat(t, k, x):
    return np.interp(np.log(x), np.log(t), k)

base_days = [0.25, 0.5, 1, 1.5, 2, 2.7, 3, 4, 5, 8, 12]
print("seed  eps  tau_sub[s]  t_end[d]")
runs = {}
for s in (1701, 1702, 1703):
    for tag in ("eps1000nm", "eps2.00um"):
        d = next(C.glob(f"packing_*seed{s}_L2mm_{tag}_perxy_T-20__*"))
        runs[s, tag] = load(d)
        t, k, tau = runs[s, tag]
        print(f"{s}  {tag}  {tau:.1f}  {t[-1]/DAY:.2f}")
tau1 = runs[1701, "eps1000nm"][2]; tau2 = runs[1701, "eps2.00um"][2]

print("\nRise k(t_end)/k(t_b) - 1 [%], eps 1 um / eps 2 um, and the difference [points]")
print(f"{'t_b [d]':>7} {'theta(1um)':>10} {'theta(2um)':>10} | " +
      " | ".join(f"seed {s}" for s in (1701, 1702, 1703)) + " | mean diff")
for b in base_days:
    tb = b * DAY; diffs = []; cells = []
    for s in (1701, 1702, 1703):
        r = []
        for tag in ("eps1000nm", "eps2.00um"):
            t, k, _ = runs[s, tag]
            tend = min(runs[s, "eps1000nm"][0][-1], runs[s, "eps2.00um"][0][-1])
            r.append(100 * (kat(t, k, tend) / kat(t, k, tb) - 1))
        diffs.append(r[1] - r[0]); cells.append(f"{r[0]:5.1f}/{r[1]:5.1f}")
    print(f"{b:7.2f} {tb/tau1:10.1f} {tb/tau2:10.1f} | " + " | ".join(cells) + f" | {np.mean(diffs):+5.1f}")
