#!/usr/bin/env python3
"""
plot_fields.py  —  Convert PetIGA sol_*.dat output to XML VTK structured-grid
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
skipped unless --force is given.  After conversion, pf.pvd is written
next to vtkOut/ so that opening it in ParaView gives a time-aware dataset.

pf.pvd is ALWAYS rewritten, and is built by globbing the .vts files that
actually exist rather than from the sol_*.dat list -- so re-running after more
timesteps have landed always extends the collection, and an entry never points
at a snapshot that was not written.  A sol_*.dat that cannot be read (routine
while a job is still writing its newest file) is warned about and skipped
instead of aborting the run.

Time values are read from outp.txt (the monitor table).  If outp.txt is absent,
the PVD file uses the step index as a proxy for time.

Usage
-----
  python plot_fields.py                 # run from the output directory
  python plot_fields.py --dir /path/to/run
  python plot_fields.py --dir /path/to/run --force
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


def _vtk_coords_from(points: np.ndarray) -> np.ndarray:
    """
    Return physical coordinates as (N_total, 3) in VTK point order (x fastest).

    VTK requires points stored with x varying fastest.  igakit's nrb.points
    has shape (Nx[, Ny[, Nz]], sdim) in C order, so the spatial axes must be
    reversed before flattening.  The component axis (last) is kept last.
    """
    grid_shape = points.shape[:-1]
    sdim       = points.shape[-1]
    dim        = len(grid_shape)

    coords = np.zeros((*grid_shape, 3))
    coords[..., :sdim] = points[..., :sdim]

    # Reverse only the spatial axes; keep component axis last.
    spatial_axes = list(range(dim - 1, -1, -1))  # e.g. [1,0] for 2-D
    axes = spatial_axes + [dim]                   # keep component last
    return np.ascontiguousarray(coords.transpose(axes)).reshape(-1, 3)


DEFAULT_MAX_AXIS = 1200   # points per axis in the written .vts


def _axis_index(n: int, stride: int) -> np.ndarray:
    """Sub-sample indices along one axis, ALWAYS keeping both endpoints.

    Plain `arr[::k]` keeps index 0 but usually drops the last point, which
    truncates the domain and biases any volume integral taken from the .vts.
    """
    if stride <= 1 or n <= 2:
        return np.arange(n)
    n_out = max(2, int(np.ceil((n - 1) / stride)) + 1)
    return np.unique(np.linspace(0, n - 1, n_out).round().astype(int))


def _stride_for(grid_shape, max_axis: int) -> int:
    """Uniform stride so no axis exceeds `max_axis` points. 0 disables."""
    if max_axis <= 0:
        return 1
    return max(1, int(np.ceil(max(grid_shape) / float(max_axis))))


def _write_vts(outfile: str, nrb, sol: np.ndarray, stride: int = 1) -> None:
    """
    Write a VTK XML Structured Grid (.vts) file.

    Parameters
    ----------
    outfile : path to write
    nrb     : igakit NURBS geometry object (provides grid shape and coordinates)
    sol     : solution array of shape (*grid_shape, ndof);
              DOFs: 0=phi_i, 1=T, 2=rho_v, 3=phi_s
    """
    full_shape = nrb.points.shape[:-1]   # (Nx,) or (Nx, Ny) or (Nx, Ny, Nz)
    dim        = len(full_shape)
    ndof       = sol.shape[-1] if sol.ndim > dim else 1

    # Sub-sample every axis by the same stride so the aspect ratio and the
    # physical coordinates stay honest.
    idx = [_axis_index(n, stride) for n in full_shape]
    pts = nrb.points[np.ix_(*idx)] if stride > 1 else nrb.points
    sol = sol[np.ix_(*idx)] if stride > 1 else sol
    grid_shape = pts.shape[:-1]

    Nx = grid_shape[0]
    Ny = grid_shape[1] if dim >= 2 else 1
    Nz = grid_shape[2] if dim >= 3 else 1

    # --- Physical coordinates (always 3-component for VTK) -----------------
    coords_vtk = _vtk_coords_from(pts)    # (N_total, 3) in VTK point order

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
            force: bool = False, max_axis: int = DEFAULT_MAX_AXIS,
            stride: int = 0):
    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        raise FileNotFoundError(f"IGA geometry file not found: {iga_path}")

    out_dir = os.path.join(run_dir, "vtkOut")
    os.makedirs(out_dir, exist_ok=True)

    nrb = PetIGA().read(iga_path)

    # Resolution cap. These meshes are sized for the SOLVER, not for viewing:
    # a 5394 x 2697 run writes ~870 MB per snapshot, 60 GB for 69 of them.
    # Sub-sampling to a viewable resolution is lossless for anything you do in
    # ParaView and, at the default cap, changes the neck and grain-volume
    # measurements by far less than their own convergence error. Pass
    # --max-axis 0 (or --stride 1) for the untouched grid.
    full_shape = nrb.points.shape[:-1]
    k = stride if stride > 0 else _stride_for(full_shape, max_axis)
    if k > 1:
        out_shape = tuple(len(_axis_index(n, k)) for n in full_shape)
        print(f"  Resolution cap: stride {k} -> "
              f"{' x '.join(map(str, full_shape))} becomes "
              f"{' x '.join(map(str, out_shape))} "
              f"({np.prod(full_shape)/np.prod(out_shape):.0f}x smaller). "
              f"Use --max-axis 0 for the full grid.")

    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol*.dat")))
    if not sol_files:
        print(f"No sol*.dat files found in '{run_dir}'")
        return

    # Build step → time mapping from outp.txt.
    time_map = _load_time_map(os.path.join(run_dir, "outp.txt"))
    if not time_map:
        print("  Warning: outp.txt not found or empty — PVD will use step "
              "index as time proxy.")

    n_fail = 0
    for infile in sol_files:
        name    = os.path.splitext(os.path.basename(infile))[0]  # "sol_00042"
        number  = name.split("l")[1]                              # "_00042"
        outfile = os.path.join(out_dir, f"solV{number}.vts")      # XML format

        if not force and os.path.isfile(outfile):
            print(f"  Skipping (exists): {outfile}")
            continue

        # One unreadable sol_*.dat must not abort the run. While a job is still
        # writing, the newest file is routinely half-flushed; letting that
        # exception escape used to kill the script before the .pvd was rewritten
        # (it is written after this loop), which is how a stale collection
        # survives a re-run.
        try:
            sol = PetIGA().read_vec(infile, nrb)
            _write_vts(outfile, nrb, sol, stride=k)
        except Exception as exc:                                 # noqa: BLE001
            n_fail += 1
            print(f"  WARNING: skipping {os.path.basename(infile)} ({exc}) — "
                  f"truncated or still being written?")
            continue
        print(f"  Written: {outfile}")

    # Build the collection from the .vts files that ACTUALLY EXIST, not from
    # the sol_*.dat list, and always rewrite it.
    #
    # The old version appended an entry per sol_*.dat before attempting the
    # conversion, so the .pvd could advertise snapshots that were never
    # written; and because it only ran when that loop completed, a single bad
    # file left the previous .pvd in place. Re-running after more timesteps
    # landed then appeared to do nothing -- ParaView kept showing however many
    # steps existed the first time.
    #
    # Globbing vtkOut/ also picks up snapshots whose sol_*.dat has since been
    # cleaned up, so the collection stays complete.
    pvd_entries = []
    for vts in sorted(glob.glob(os.path.join(out_dir, "solV_*.vts"))):
        stem = os.path.splitext(os.path.basename(vts))[0]        # "solV_00042"
        try:
            step = int(stem.split("_")[1])
        except (IndexError, ValueError):
            continue
        pvd_entries.append((time_map.get(step, float(step)),
                            os.path.join("vtkOut", os.path.basename(vts))))

    if pvd_entries:
        _write_pvd(os.path.join(run_dir, "pf.pvd"), pvd_entries)
        print(f"  PVD lists {len(pvd_entries)} snapshot(s)"
              + (f"; {n_fail} sol file(s) skipped" if n_fail else ""))
    else:
        print("  No .vts files present — pf.pvd not written.")


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
    p.add_argument("--max-axis", type=int, default=DEFAULT_MAX_AXIS,
                   help="Cap points per axis in the written .vts; 0 = no cap. "
                        "Solver meshes are far finer than any display needs.")
    p.add_argument("--stride", type=int, default=0,
                   help="Explicit sub-sampling stride (overrides --max-axis)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(run_dir=args.dir, iga_file=args.iga, force=args.force,
            max_axis=args.max_axis, stride=args.stride)
