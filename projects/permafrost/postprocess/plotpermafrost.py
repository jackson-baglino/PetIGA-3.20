#!/usr/bin/env python3
"""
plotpermafrost.py  —  Convert PetIGA sol_*.dat output to VTK files and write
a ParaView Data (.pvd) collection index so that ParaView understands the
simulation time at each snapshot.

Exports four DOFs plus one derived field:
  IcePhase     — φ_i          (DOF 0)
  Temperature  — T            (DOF 1)
  VaporDensity — ρ_v          (DOF 2)
  SedPhase     — φ_s          (DOF 3)
  AirPhase     — 1 − φ_i − φ_s  (derived, clipped to [0, 1])

VTK files are written to ./vtkOut/.  Existing files are skipped unless
--force is given.  After conversion, permafrost.pvd is written next to
vtkOut/ so that opening it in ParaView gives a time-aware dataset.

Time values are read from outp.txt (the monitor table).  If outp.txt is
absent, the PVD file uses the step index as a proxy for time.

Usage
-----
  python plotpermafrost.py                 # run from the output directory
  python plotpermafrost.py --dir /path/to/run
  python plotpermafrost.py --dir /path/to/run --force
"""

import argparse
import glob
import os
import re

import numpy as np
from igakit.io import PetIGA, VTK


# ---------------------------------------------------------------------------
# Time-map: step index → simulation time, parsed from outp.txt
# ---------------------------------------------------------------------------

# outp.txt has two kinds of pipe-delimited rows:
#   domain rows (9 fields):  STEP | TIME | DT | TOT_ICE | ... | TRIPL_JUNC
#   SNES rows  (13 fields):  it   | fnorm | n0 | r0 | ...
# We want only domain rows.
_OUTP_ROW_RE   = re.compile(r"^\s*(\d+)\s*\|(.+)$")
_DOMAIN_NFIELDS = 9


def _load_time_map(outp_path: str) -> dict:
    """Return {step_index: simulation_time_s} from outp.txt."""
    time_map = {}
    if not os.path.isfile(outp_path):
        return time_map
    with open(outp_path) as fh:
        for line in fh:
            m = _OUTP_ROW_RE.match(line)
            if not m:
                continue
            fields = m.group(2).split("|")
            if len(fields) != _DOMAIN_NFIELDS:  # skip SNES/Newton rows
                continue
            try:
                step = int(m.group(1))
                t    = float(fields[0].strip())   # TIME [s]
                time_map[step] = t
            except ValueError:
                continue
    return time_map


# ---------------------------------------------------------------------------
# PVD writer
# ---------------------------------------------------------------------------

def _write_pvd(pvd_path: str, entries: list) -> None:
    """
    Write a ParaView Data (.pvd) XML collection file.

    entries : list of (timestep_float, vtk_path_relative_to_pvd)
    """
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        '  <Collection>',
    ]
    for t, rel_path in sorted(entries):
        lines.append(
            f'    <DataSet timestep="{t:.6e}" part="0" file="{rel_path}"/>'
        )
    lines += ['  </Collection>', '</VTKFile>']
    with open(pvd_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  PVD collection written: {pvd_path}")


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------

def convert(run_dir: str = ".", iga_file: str = "igasol.dat",
            force: bool = False):
    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        raise FileNotFoundError(f"IGA geometry file not found: {iga_path}")

    out_dir = os.path.join(run_dir, "vtkOut")
    os.makedirs(out_dir, exist_ok=True)

    nrb = PetIGA().read(iga_path)

    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol*.dat")))
    if not sol_files:
        print(f"No sol*.dat files found in '{run_dir}'")
        return

    # Build step → time mapping from outp.txt (best-effort; falls back to
    # step index if outp.txt is absent or a step is missing).
    time_map = _load_time_map(os.path.join(run_dir, "outp.txt"))
    if not time_map:
        print("  Warning: outp.txt not found or empty — PVD will use step "
              "index as time proxy.")

    pvd_entries = []   # (time, relative_vtk_path) for every VTK file

    for infile in sol_files:
        name    = os.path.splitext(os.path.basename(infile))[0]  # "sol_00042"
        number  = name.split("l")[1]                              # "_00042"
        step    = int(name.split("_")[1])                         # 42
        outfile = os.path.join(out_dir, f"solV{number}.vtk")

        # Relative path from run_dir (where the PVD will live) to the VTK file
        rel_vtk = os.path.join("vtkOut", f"solV{number}.vtk")

        # Time for this snapshot
        t = time_map.get(step, float(step))   # fall back to step index
        pvd_entries.append((t, rel_vtk))

        if not force and os.path.isfile(outfile):
            print(f"  Skipping (exists): {outfile}")
            continue

        sol = PetIGA().read_vec(infile, nrb)

        scalars = {
            "IcePhase":    0,
            "Temperature": 1,
            "VaporDensity": 2,
            "SedPhase":    3,
        }
        # Only export DOFs present in the solution vector
        if sol.ndim < 2 or sol.shape[-1] < 4:
            scalars.pop("SedPhase")

        # Append AirPhase = 1 - φ_i - φ_s as a derived field
        if sol.ndim >= 2 and sol.shape[-1] >= 4:
            air = np.clip(1.0 - sol[..., 0:1] - sol[..., 3:4], 0.0, 1.0)
            sol = np.concatenate([sol, air], axis=-1)
            scalars["AirPhase"] = sol.shape[-1] - 1

        VTK().write(outfile, nrb, fields=sol, scalars=scalars)
        print(f"  Written: {outfile}")

    # Always (re)write the PVD so it covers all snapshots, including any that
    # were already present and skipped above.
    if pvd_entries:
        pvd_path = os.path.join(run_dir, "permafrost.pvd")
        _write_pvd(pvd_path, pvd_entries)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert PetIGA output to VTK and write a PVD time index."
    )
    p.add_argument("--dir",   default=".",          help="Run directory (default: .)")
    p.add_argument("--iga",   default="igasol.dat", help="IGA geometry file")
    p.add_argument("--force", action="store_true",  help="Overwrite existing VTK files")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(run_dir=args.dir, iga_file=args.iga, force=args.force)
