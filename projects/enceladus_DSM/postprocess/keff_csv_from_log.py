#!/usr/bin/env python3
"""Rebuild a k_eff CSV from a solver log.

    venv_enceladus/bin/python postprocess/keff_csv_from_log.py <run_dir> [...]

WHY THIS EXISTS. `-keff_replay` writes its CSV next to the run being REPLAYED,
not into the job's own output folder (src/keff.c:377-381) -- which is the right
place to look for it later, but means a replay batch can be downloaded without
its results: the numbers are on $SCRATCH in the source directory, while what
came back is the job folder. The per-sample values are in the log either way,
because KeffSample prints every one of them.

Parses the blocks

    [keff] step 191    t = 2.873318e+04 s   phi_bar = 0.675798   k_iso = 4.716e-01 W/m/K   (200 its, 96.62 s)
           k = [ 4.913376583e-01  2.468149017e-02 ]
               [ 2.468149017e-02  4.518856950e-01 ]

and writes k_eff_<law>.csv with the schema KeffCSVAppend uses, so the result is
indistinguishable from the file the solver would have written. The law comes
from the startup banner ("band interpolation: tensor ..."); without it the file
is named k_eff_from_log.csv rather than guessed at.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

STEP_RE = re.compile(
    r"\[keff\]\s+step\s+(\d+)\s+t\s*=\s*([-\d.eE+]+)\s*s\s+"
    r"phi_bar\s*=\s*([-\d.eE+]+)\s+k_iso\s*=\s*([-\d.eE+]+)\s*W/m/K"
    r"\s*\((\d+)\s+its,\s*([\d.]+)\s*s\)")
ROW_RE  = re.compile(r"\[\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s*\]")
LAW_RE  = re.compile(r"band interpolation:\s*(\w+)")


def parse(text: str):
    law = None
    m = LAW_RE.search(text)
    if m:
        law = m.group(1)
    lines = text.splitlines()
    out = []
    for i, ln in enumerate(lines):
        m = STEP_RE.search(ln)
        if not m:
            continue
        step, t, phi, kiso, its, wall = m.groups()
        # The two tensor rows follow immediately; absent in 3D or if truncated.
        rows = []
        for j in (i + 1, i + 2):
            if j < len(lines):
                r = ROW_RE.search(lines[j])
                if r:
                    rows.append((float(r.group(1)), float(r.group(2))))
        if len(rows) == 2:
            k00, k01 = rows[0]
            k10, k11 = rows[1]
        else:
            k00 = k11 = float(kiso)
            k01 = k10 = 0.0
        out.append((int(step), float(t), k00, k01, k10, k11,
                    float(phi), float(kiso), int(its), 2, float(wall)))
    return law, out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--glob", default="outp*.txt")
    args = ap.parse_args()

    n_ok = 0
    for d in args.dirs:
        logs = sorted(d.glob(args.glob)) if d.is_dir() else [d]
        # Prefer the per-job leg file: outp.txt is cumulative and, across a
        # resume, holds more than this job's samples.
        legs = [p for p in logs if p.name.startswith("outp_job")]
        logs = legs or logs
        law, rows = None, []
        for p in logs:
            lw, rs = parse(p.read_text(errors="ignore"))
            law = law or lw
            rows.extend(rs)
        if not rows:
            print(f"  no [keff] samples in {d}")
            continue
        rows.sort(key=lambda r: r[1])
        name = f"k_eff_{law}.csv" if law and law != "arithmetic" else \
               ("k_eff.csv" if law == "arithmetic" else "k_eff_from_log.csv")
        dest = (d if d.is_dir() else d.parent) / name
        with dest.open("w") as fh:
            fh.write("step,time,k_00,k_01,k_10,k_11,phi_bar,k_iso,"
                     "ksp_its,ksp_reason,wall_s\n")
            for r in rows:
                fh.write(f"{r[0]},{r[1]:.12e},{r[2]:.12e},{r[3]:.12e},{r[4]:.12e},"
                         f"{r[5]:.12e},{r[6]:.12e},{r[7]:.12e},{r[8]},{r[9]},{r[10]:.3f}\n")
        print(f"  {dest}  ({len(rows)} samples, law={law})")
        n_ok += 1
    return 0 if n_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
