#!/usr/bin/env python3
"""Collect the disk-array ice-fraction sweep from a downloaded batch folder.

    venv_enceladus/bin/python studies/keff_sharp_limit/disk_sweep/collect_sweep.py <batch_dir>

Reads each <geom>__<exp>__<label>/ run folder, takes its k_eff row (k_eff.csv
for the arithmetic law, k_eff_tensor.csv for the tensor law), joins it with
sweep_manifest.csv, and writes keff_disk_sweep.csv here.

Checks, each printed and counted, none silently skipped:
  * every manifest label has a run folder with a k_eff row;
  * k_00 = k_11 and k_01 = k_10 = 0 (square symmetry), to 1e-3;
  * phi_bar equals the exact diffuse-disk fraction f(1 + pi^2 eps^2 / 3R^2)
    to 1e-4, i.e. the solver built the disk that was asked for. The formula
    is for an ISOLATED disk; the IC is one disk per cell, so its tanh tail is
    cut at the cell edge. That costs O(exp(-(L/2 - R)/eps)) of the ice: 5.6e-5
    relative at R = 400 um, eps/R = 0.04 (gap/2 = 6 eps), 7e-7 at R = 350 um,
    eps/R = 0.04, and below 4e-8 in every other run.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "disk"))
import keff_disk_analytic as ana                          # noqa: E402

L = 1.0e-3


def last_row(path: Path) -> dict | None:
    if not path.is_file():
        return None
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    return rows[-1] if rows else None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch_dir", type=Path)
    a = ap.parse_args()

    with (HERE / "sweep_manifest.csv").open() as fh:
        manifest = list(csv.DictReader(fh))

    out, problems = [], []
    for m in manifest:
        hits = sorted(a.batch_dir.glob(f"*__{m['label']}"))
        if len(hits) != 1:
            problems.append(f"{m['label']}: {len(hits)} run folders")
            continue
        name = "k_eff.csv" if m["law"] == "arithmetic" else "k_eff_tensor.csv"
        row = last_row(hits[0] / name)
        if row is None:
            problems.append(f"{m['label']}: no {name}")
            continue
        R, eps = float(m["R"]), float(m["eps"])
        f = ana.area_fraction(R, L)
        k00, k11 = float(row["k_00"]), float(row["k_11"])
        k01, k10 = float(row["k_01"]), float(row["k_10"])
        k = 0.5 * (k00 + k11)
        ks = ana.k_sharp(f)
        phi_exact = f * (1.0 + math.pi**2 * eps**2 / (3.0 * R**2))
        phi_bar = float(row["phi_bar"])
        if abs(k00 - k11) / k00 > 1e-3 or max(abs(k01), abs(k10)) / k00 > 1e-3:
            problems.append(f"{m['label']}: not square-symmetric ({k00:.6g}, {k11:.6g})")
        if abs(phi_bar / phi_exact - 1.0) > 1e-4:
            problems.append(f"{m['label']}: phi_bar {phi_bar:.8g} vs exact {phi_exact:.8g}")
        out.append(dict(law=m["law"], R=R, f=f, eps_over_R=float(m["eps_over_R"]),
                        eps=eps, N=int(m["N"]), k=k, k_sharp=ks, rel_err=k / ks - 1.0,
                        phi_bar=phi_bar, ksp_its=row.get("ksp_its", ""),
                        run_dir=hits[0].name))

    out.sort(key=lambda r: (r["law"], r["eps_over_R"], r["R"]))
    dst = HERE / "keff_disk_sweep.csv"
    with dst.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]) if out else ["law"],
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    print(f"{len(out)}/{len(manifest)} runs collected -> {dst}")
    print(f"\n  {'law':<10} {'eps/R':>6} " + " ".join(f"f={r:>5.3f}" for r in
          sorted({r['f'] for r in out})))
    for law in ("arithmetic", "tensor"):
        for ratio in sorted({r["eps_over_R"] for r in out}, reverse=True):
            cells = [r for r in out if r["law"] == law and r["eps_over_R"] == ratio]
            print(f"  {law:<10} {ratio:6g} " +
                  " ".join(f"{r['rel_err']:+7.2%}" for r in cells))
    if problems:
        print(f"\n{len(problems)} problem(s):")
        for p in problems:
            print("  " + p)
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
