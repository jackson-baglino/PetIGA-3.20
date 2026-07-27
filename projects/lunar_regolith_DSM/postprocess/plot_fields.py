#!/usr/bin/env python3
"""
plotpermafrost.py  —  Convert PetIGA sol_*.dat output to XML VTK structured-grid
files (.vts) and write a ParaView Data (.pvd) collection so that ParaView
understands the simulation time at each snapshot.

The PVD reader in ParaView (vtkXMLCollectionReader) requires XML-format VTK
datasets.  Legacy binary .vtk files are NOT supported inside a PVD collection.
This script writes proper VTK XML Structured Grid (.vts) files, which the PVD
reader can open without errors.

Exports four DOFs plus one derived field:
  IcePhase     — φ_i               (DOF 0)
  Temperature  — T                 (DOF 1)
  VaporDensity — ρ_v               (DOF 2)
  SedPhase     — φ_s               (DOF 3)
  AirPhase     — 1 − φ_i − φ_s    (derived, NOT clipped — see note below)

Note on AirPhase: previously this field was np.clip()-ed to [0, 1] before
write. That hid small out-of-bound excursions (e.g. phi_air = -1e-3 from
the B-spline overshoot at sharp ice-sed interfaces) and made the .vts
inconsistent with the unclipped BOUNDS check printed by the simulation
monitor. The clip has been removed so the .vts faithfully reflects the
actual field; use ParaView's colormap range controls if you want a
visually clean display.

VTK files are written to ./vtkOut/ with extension .vts.  Existing files are
skipped unless --force is given.  After conversion, permafrost.pvd is written
next to vtkOut/ so that opening it in ParaView gives a time-aware dataset.

Time values are read from outp.txt (the monitor table).  If outp.txt is absent,
the PVD file uses the step index as a proxy for time.

Usage
-----
  python plotpermafrost.py                 # run from the output directory
  python plotpermafrost.py --dir /path/to/run
  python plotpermafrost.py --dir /path/to/run --force
"""

import argparse
import base64
import glob
import os
import re
import struct

import numpy as np
from igakit.io import PetIGA


# ---------------------------------------------------------------------------
# Time-map: step index → simulation time, parsed from outp.txt
# ---------------------------------------------------------------------------

# outp.txt has two kinds of pipe-delimited rows:
#   domain rows (8 fields): STEP | TIME | DT | TOT_ICE | TOT_AIR | TEMP |
#                           TOT_RHOV | I-A INTERF | TOTAL_MASS
#   SNES rows  (10 fields): it   | fnorm | n0 | r0 | s0 | n1 | r1 | s1 | n2 | r2 | s2
# Both start with "<small int> |", so we disambiguate by field count, not by
# the leading index alone -- _DOMAIN_NFIELDS must track monitoring.c's actual
# column count exactly (it previously said 10, a stale value from when the
# domain row also had TOT_SED/TRIPL_JUNC columns; that made this regex match
# SNES rows instead, since those happen to also have 10 fields after the
# leading index -- silently corrupting the time map for nearly every step
# with the SNES iteration's fnorm value instead of the real simulated time).
_OUTP_ROW_RE    = re.compile(r"^\s*(\d+)\s*\|(.+)$")
_DOMAIN_NFIELDS = 8  # TIME, DT, TOT_ICE, TOT_AIR, TEMP, TOT_RHOV, I-A INTERF, TOTAL_MASS


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
            if len(fields) != _DOMAIN_NFIELDS:
                continue
            try:
                step = int(m.group(1))
                t    = float(fields[0].strip())
                time_map[step] = t
            except ValueError:
                continue
    return time_map


# ---------------------------------------------------------------------------
# VTK XML writer helpers
# ---------------------------------------------------------------------------

def _b64_block(arr: np.ndarray) -> str:
    """
    Encode a numpy array as a base64 binary block for VTK XML format="binary".

    Layout: [UInt64 byte-count header][raw float64 data], base64 encoded.
    Matches header_type="UInt64" byte_order="LittleEndian" in the VTKFile tag.
    """
    data = np.ascontiguousarray(arr.ravel(), dtype='<f8').tobytes()
    header = struct.pack('<Q', len(data))   # 8-byte little-endian uint64
    return base64.b64encode(header + data).decode('ascii')


def _vtk_scalar(arr: np.ndarray) -> np.ndarray:
    """
    Reorder a scalar field from igakit's (Nx[, Ny[, Nz]]) layout to VTK's
    point linearisation where x (first index) varies fastest.

    igakit stores points as arr[ix, iy] in C order, so a plain .ravel()
    gives ix-slow, iy-fast — the opposite of VTK's expected ix-fast order.
    Reversing the axis order before flattening corrects this.
    """
    ndim = arr.ndim
    if ndim == 1:
        return arr.ravel()
    axes = list(range(ndim - 1, -1, -1))   # e.g. [1,0] for 2-D, [2,1,0] for 3-D
    return np.ascontiguousarray(arr.transpose(axes)).ravel()


def _vtk_coords(nrb) -> np.ndarray:
    """
    Return physical coordinates as (N_total, 3) in VTK point order (x fastest).

    VTK requires points stored with x varying fastest.  igakit's nrb.points
    has shape (Nx[, Ny[, Nz]], sdim) in C order, so the spatial axes must be
    reversed before flattening.  The component axis (last) is kept last.
    """
    grid_shape = nrb.points.shape[:-1]
    sdim       = nrb.points.shape[-1]
    dim        = len(grid_shape)

    coords = np.zeros((*grid_shape, 3))
    coords[..., :sdim] = nrb.points[..., :sdim]

    # Reverse only the spatial axes; keep component axis last.
    spatial_axes = list(range(dim - 1, -1, -1))  # e.g. [1,0] for 2-D
    axes = spatial_axes + [dim]                   # keep component last
    return np.ascontiguousarray(coords.transpose(axes)).reshape(-1, 3)


def _write_vts(outfile: str, nrb, sol: np.ndarray) -> None:
    """
    Write a VTK XML Structured Grid (.vts) file.

    Parameters
    ----------
    outfile : path to write
    nrb     : igakit NURBS geometry object (provides grid shape and coordinates)
    sol     : solution array of shape (*grid_shape, ndof);
              DOFs: 0=phi_i, 1=T, 2=rho_v, 3=phi_s
    """
    grid_shape = nrb.points.shape[:-1]   # (Nx,) or (Nx, Ny) or (Nx, Ny, Nz)
    dim        = len(grid_shape)
    ndof       = sol.shape[-1] if sol.ndim > dim else 1

    Nx = grid_shape[0]
    Ny = grid_shape[1] if dim >= 2 else 1
    Nz = grid_shape[2] if dim >= 3 else 1

    # --- Physical coordinates (always 3-component for VTK) -----------------
    coords_vtk = _vtk_coords(nrb)    # (N_total, 3) in VTK point order

    # --- Field scalars -------------------------------------------------------
    fields = {}
    if ndof >= 1:
        fields["IcePhase"]     = sol[..., 0]
    if ndof >= 2:
        fields["Temperature"]  = sol[..., 1]
    if ndof >= 3:
        fields["VaporDensity"] = sol[..., 2]
    if ndof >= 4:
        fields["SedPhase"]     = sol[..., 3]
        fields["AirPhase"]     = 1.0 - sol[..., 0] - sol[..., 3]

    # Reorder each scalar field to VTK point order
    fields_vtk = {name: _vtk_scalar(arr) for name, arr in fields.items()}

    # --- XML assembly -------------------------------------------------------
    extent = f"0 {Nx-1} 0 {Ny-1} 0 {Nz-1}"

    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="StructuredGrid" version="0.1" '
        'byte_order="LittleEndian" header_type="UInt64">',
        f'  <StructuredGrid WholeExtent="{extent}">',
        f'    <Piece Extent="{extent}">',
        '      <Points>',
        '        <DataArray type="Float64" NumberOfComponents="3" format="binary">',
        f'          {_b64_block(coords_vtk)}',
        '        </DataArray>',
        '      </Points>',
        '      <PointData>',
    ]
    for name, arr in fields_vtk.items():
        lines.append(
            f'        <DataArray type="Float64" Name="{name}" format="binary">'
            f'{_b64_block(arr)}</DataArray>'
        )
    lines += [
        '      </PointData>',
        '    </Piece>',
        '  </StructuredGrid>',
        '</VTKFile>',
    ]

    with open(outfile, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


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

    # Build step → time mapping from outp.txt.
    time_map = _load_time_map(os.path.join(run_dir, "outp.txt"))
    if not time_map:
        print("  Warning: outp.txt not found or empty — PVD will use step "
              "index as time proxy.")

    pvd_entries = []   # (time, relative_vts_path)

    for infile in sol_files:
        name    = os.path.splitext(os.path.basename(infile))[0]  # "sol_00042"
        number  = name.split("l")[1]                              # "_00042"
        step    = int(name.split("_")[1])                         # 42
        outfile = os.path.join(out_dir, f"solV{number}.vts")      # XML format
        rel_vts = os.path.join("vtkOut", f"solV{number}.vts")

        t = time_map.get(step, float(step))
        pvd_entries.append((t, rel_vts))

        if not force and os.path.isfile(outfile):
            print(f"  Skipping (exists): {outfile}")
            continue

        sol = PetIGA().read_vec(infile, nrb)
        _write_vts(outfile, nrb, sol)
        print(f"  Written: {outfile}")

    if pvd_entries:
        pvd_path = os.path.join(run_dir, "permafrost.pvd")
        _write_pvd(pvd_path, pvd_entries)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Convert PetIGA output to VTK XML (.vts) and write a "
                    "PVD time-collection index for ParaView."
    )
    p.add_argument("--dir",   default=".",          help="Run directory (default: .)")
    p.add_argument("--iga",   default="igasol.dat", help="IGA geometry file")
    p.add_argument("--force", action="store_true",  help="Overwrite existing VTK files")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(run_dir=args.dir, iga_file=args.iga, force=args.force)
