#!/usr/bin/env python3
"""
plotDSM.py — Convert PetIGA sol_*.dat snapshots to VTK/VTS files and write a
ParaView PVD collection file.

Reads the IGA geometry from igasol.dat and converts each sol_NNNNN.dat
snapshot to a VTK file under vtkOut/. After all snapshots are processed,
writes dsm.pvd (a ParaView data collection) so all snapshots can be opened
at once as a time series.

Usage
-----
  # From inside a run folder:
  python3 postprocess/plotDSM.py

  # Specify a different run directory:
  python3 postprocess/plotDSM.py --dir /path/to/run

  # Force regeneration of existing VTK files:
  python3 postprocess/plotDSM.py --force

Fields written per snapshot (indices match DOF order):
  IcePhase    — col 0 — φ_i ∈ [0, 1]
  Temperature — col 1 — T   [°C relative to T₀]
  VaporDensity — col 2 — ρ_v [kg/m³]
"""

import argparse
import glob
import os
import re
import sys

try:
    from igakit.io import PetIGA, VTK
except ImportError:
    sys.exit("ERROR: igakit is not installed.\nInstall with:  pip install igakit")


# ---- PVD writer ----

_PVD_HEADER = """\
<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
"""
_PVD_FOOTER = """\
  </Collection>
</VTKFile>
"""


def write_pvd(pvd_path: str, entries: list):
    """Write a ParaView PVD collection.

    entries: list of (time, vtk_rel_path) tuples, sorted by time.
    """
    with open(pvd_path, "w") as f:
        f.write(_PVD_HEADER)
        for t, relpath in entries:
            f.write(f'    <DataSet timestep="{t:.6e}" group="" part="0" file="{relpath}"/>\n')
        f.write(_PVD_FOOTER)
    print(f"PVD collection written: {pvd_path}")


def parse_step(filename: str) -> int:
    """Extract the integer step number from sol_NNNNN.dat."""
    m = re.search(r"(\d+)", os.path.basename(filename))
    return int(m.group(1)) if m else -1


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir",   default=".",
                   help="Run directory containing igasol.dat and sol_*.dat (default: .)")
    p.add_argument("--force", action="store_true",
                   help="Regenerate VTK files even if they already exist")
    p.add_argument("--iga",   default="igasol.dat",
                   help="IGA geometry file name (default: igasol.dat)")
    p.add_argument("--vtk-dir", default="vtkOut",
                   help="Sub-directory for VTK output (default: vtkOut)")
    p.add_argument("--pvd",   default="dsm.pvd",
                   help="PVD collection file name in --dir (default: dsm.pvd)")
    p.add_argument("--ssa",   default="SSA_evo.dat",
                   help="SSA_evo.dat for time-axis data (default: SSA_evo.dat)")
    args = p.parse_args()

    run_dir = os.path.abspath(args.dir)
    iga_path = os.path.join(run_dir, args.iga)
    vtk_dir  = os.path.join(run_dir, args.vtk_dir)
    pvd_path = os.path.join(run_dir, args.pvd)
    ssa_path = os.path.join(run_dir, args.ssa)

    if not os.path.isfile(iga_path):
        sys.exit(f"ERROR: IGA geometry file not found: {iga_path}")

    os.makedirs(vtk_dir, exist_ok=True)

    print(f"Reading IGA geometry: {iga_path}")
    nrb = PetIGA().read(iga_path)

    # ---- Build step→time mapping from SSA_evo.dat ----
    step_to_time = {}
    if os.path.isfile(ssa_path):
        import numpy as np
        try:
            ssa = np.genfromtxt(ssa_path, dtype=float, comments="#",
                                 invalid_raise=False)
            if ssa.ndim == 1:
                ssa = ssa[np.newaxis, :]
            if ssa.shape[1] >= 4:
                for row in ssa:
                    if not any(map(lambda x: x != x, row[:4])):  # no NaN
                        step_to_time[int(row[3])] = float(row[2])
        except Exception:
            pass

    # ---- Sort solution files by step number ----
    sol_files = sorted(
        glob.glob(os.path.join(run_dir, "sol_*.dat")),
        key=parse_step
    )
    if not sol_files:
        sol_files = sorted(
            glob.glob(os.path.join(run_dir, "sol*.dat")),
            key=parse_step
        )
    if not sol_files:
        sys.exit(f"ERROR: No sol_*.dat files found in {run_dir}")

    print(f"Found {len(sol_files)} snapshot(s)")

    pvd_entries = []

    for sol_file in sol_files:
        step = parse_step(sol_file)
        base = os.path.splitext(os.path.basename(sol_file))[0]
        vtk_name = f"solV{step:05d}.vtk"
        vtk_path = os.path.join(vtk_dir, vtk_name)
        vtk_rel  = os.path.join(args.vtk_dir, vtk_name)

        if not args.force and os.path.isfile(vtk_path):
            print(f"  Skipping {vtk_name} (already exists; use --force to regenerate)")
        else:
            print(f"  Writing {vtk_name} ...")
            sol = PetIGA().read_vec(sol_file, nrb)
            VTK().write(
                vtk_path, nrb,
                fields=sol,
                scalars={
                    "IcePhase":    0,
                    "Temperature": 1,
                    "VaporDensity": 2,
                },
            )

        t = step_to_time.get(step, float(step))
        pvd_entries.append((t, vtk_rel))

    pvd_entries.sort(key=lambda e: e[0])
    write_pvd(pvd_path, pvd_entries)
    print(f"\nDone. Open {pvd_path} in ParaView to view the time series.")


if __name__ == "__main__":
    main()
