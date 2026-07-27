#!/usr/bin/env python3
"""
plot_permafrost_highres.py — Convert PetIGA sol_*.dat output to VTK XML
structured-grid (.vts) files like plotpermafrost.py, but evaluate the TRUE
NURBS solution field at a dense set of points per element instead of writing
the raw control points.

Why: plotpermafrost.py writes nrb.points (the control net) directly as VTK
nodes, one-to-one with the solution DOFs. ParaView then draws straight lines
between those nodes. But the actual IGA solution is a smooth degree-(P,P),
C^{P-1} field *within* each element -- control points are not interpolation
nodes, so connecting them linearly throws away exactly the curvature that
the basis functions are there to represent, and shows up as visible
faceting on contours even when the underlying solution is fine.

This script instead evaluates the real field via igakit's NURBS `fields=`
mechanism (the same basis functions used for the geometry map), sampling
`--n-per-elem` points per element in each direction, and writes that dense,
smoothly-interpolated grid instead. No solver re-run, no mesh change --
this is pure post-processing of the existing sol_*.dat output.

Cost note: points grow ~ (n_per_elem)^2 relative to the element grid, which
is already close to the control-point grid. n_per_elem=4 (default) gives
visibly smooth contours for a P=2/C1 field without excessive file size;
push higher only if you still see faceting.

Usage
-----
  python plot_permafrost_highres.py                 # run from the output directory
  python plot_permafrost_highres.py --dir /path/to/run --n-per-elem 4
  python plot_permafrost_highres.py --dir /path/to/run --force
"""

import argparse
import glob
import os

import numpy as np
from igakit.io import PetIGA

from plotpermafrost import (
    _b64_block,
    _vtk_scalar,
    _load_time_map,
    _write_pvd,
)


def _dense_uv(nrb, n_per_elem: int):
    """Per-axis dense parametric sample points: n_per_elem points per
    element (true knot-span breaks, not control points), deduplicated at
    shared element boundaries -- same convention as
    preprocess/build_geometry_multi_grain.py's write_vtk()."""
    uv = []
    for axis in range(nrb.dim):
        breaks = nrb.breaks(axis)
        pts = np.unique(np.concatenate([
            np.linspace(breaks[i], breaks[i + 1], n_per_elem, endpoint=False)
            for i in range(len(breaks) - 1)
        ] + [[breaks[-1]]]))
        uv.append(pts)
    return uv


def _vtk_coords_dense(C: np.ndarray) -> np.ndarray:
    """Same x-fastest VTK point reordering as plotpermafrost.py's
    _vtk_coords(), but operating on a precomputed dense coordinate array
    (shape (nu, nv[, nw], sdim)) instead of nrb.points."""
    grid_shape = C.shape[:-1]
    sdim       = C.shape[-1]
    dim        = len(grid_shape)

    coords = np.zeros((*grid_shape, 3))
    coords[..., :sdim] = C[..., :sdim]

    spatial_axes = list(range(dim - 1, -1, -1))
    axes = spatial_axes + [dim]
    return np.ascontiguousarray(coords.transpose(axes)).reshape(-1, 3)


def _write_vts_dense(outfile: str, C: np.ndarray, F: np.ndarray) -> None:
    """Write a dense-sampled VTK XML Structured Grid (.vts) file.

    C : dense physical coordinates, shape (nu, nv[, nw], 3)
    F : dense solution fields,      shape (nu, nv[, nw], ndof)
    """
    grid_shape = C.shape[:-1]
    dim        = len(grid_shape)
    ndof       = F.shape[-1]

    Nx = grid_shape[0]
    Ny = grid_shape[1] if dim >= 2 else 1
    Nz = grid_shape[2] if dim >= 3 else 1

    coords_vtk = _vtk_coords_dense(C)

    fields = {}
    if ndof >= 1:
        fields["IcePhase"]     = F[..., 0]
    if ndof >= 2:
        fields["Temperature"]  = F[..., 1]
    if ndof >= 3:
        fields["VaporDensity"] = F[..., 2]
    if ndof >= 4:
        fields["SedPhase"]     = F[..., 3]
        fields["AirPhase"]     = 1.0 - F[..., 0] - F[..., 3]

    fields_vtk = {name: _vtk_scalar(arr) for name, arr in fields.items()}

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


def convert(run_dir: str = ".", iga_file: str = "igasol.dat",
            force: bool = False, n_per_elem: int = 4, steps=None):
    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        raise FileNotFoundError(f"IGA geometry file not found: {iga_path}")

    out_dir = os.path.join(run_dir, "vtkOut_highres")
    os.makedirs(out_dir, exist_ok=True)

    nrb = PetIGA().read(iga_path)
    u, v = _dense_uv(nrb, n_per_elem)
    print(f"  Dense grid: {len(u)} x {len(v)} points "
          f"({n_per_elem} per element, vs {nrb.points.shape[0]} x "
          f"{nrb.points.shape[1]} control points)")

    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol*.dat")))
    if not sol_files:
        print(f"No sol*.dat files found in '{run_dir}'")
        return
    if steps is not None:
        steps = set(steps)
        sol_files = [f for f in sol_files
                     if int(os.path.splitext(os.path.basename(f))[0].split("_")[1]) in steps]
        if not sol_files:
            print(f"No sol*.dat files match --steps {sorted(steps)}")
            return

    time_map = _load_time_map(os.path.join(run_dir, "outp.txt"))
    if not time_map:
        print("  Warning: outp.txt not found or empty — PVD will use step "
              "index as time proxy.")

    pvd_entries = []

    for infile in sol_files:
        name    = os.path.splitext(os.path.basename(infile))[0]
        number  = name.split("l")[1]
        step    = int(name.split("_")[1])
        outfile = os.path.join(out_dir, f"solV{number}.vts")
        rel_vts = os.path.join("vtkOut_highres", f"solV{number}.vts")

        t = time_map.get(step, float(step))
        pvd_entries.append((t, rel_vts))

        if not force and os.path.isfile(outfile):
            print(f"  Skipping (exists): {outfile}")
            continue

        sol = PetIGA().read_vec(infile, nrb)
        C, F = nrb(u, v, fields=sol)
        _write_vts_dense(outfile, C, F)
        print(f"  Written: {outfile}")

    if pvd_entries:
        pvd_path = os.path.join(run_dir, "permafrost_highres.pvd")
        _write_pvd(pvd_path, pvd_entries)


def parse_args():
    p = argparse.ArgumentParser(
        description="Convert PetIGA output to VTK XML (.vts) at a dense, "
                    "true-NURBS-evaluated resolution for smooth contours, "
                    "and write a PVD time-collection index for ParaView."
    )
    p.add_argument("--dir",        default=".",          help="Run directory (default: .)")
    p.add_argument("--iga",        default="igasol.dat", help="IGA geometry file")
    p.add_argument("--force",      action="store_true",  help="Overwrite existing VTK files")
    p.add_argument("--n-per-elem", type=int, default=4,
                    help="Sample points per element per direction (default: 4)")
    p.add_argument("--steps", type=int, nargs="+", default=None,
                    help="Convert only these step indices (default: all). "
                         "Dense output is large (~16x a control-point .vts "
                         "per step at the default --n-per-elem) -- prefer "
                         "this for spot-checking a few snapshots rather "
                         "than converting a full time series.")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    convert(run_dir=args.dir, iga_file=args.iga, force=args.force,
            n_per_elem=args.n_per_elem, steps=args.steps)
