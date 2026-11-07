"""Synthetic run folders for the postprocessing tests.

Every file a run folder can hold is written here in the solver's own format
(the format strings are copied from src/monitoring.c and plot_fields.py), so a
test can build a run whose correct answer is known exactly:

    run = make_run(tmp, geom="g", exp="e")                 # folder + staged opts
    write_ssa(run, t, ice, rhov=...)                       # SSA_evo.dat
    write_bflux(run, t, {"x0": (J, Q)})                    # bflux.dat
    write_resolution(run, t, vn, eps=..., beta0=..., d0=...)
    write_outp(run, steps, t)                              # monitor table
    write_iga_run(run, Lx, Ly, nx, ny, snaps)              # igasol.dat + sol_*.dat
    write_vts_run(run, X, Y, snaps)                        # vtkOut/*.vts + pf.pvd

Analytic shapes with the model's own profile phi = (1 + tanh(d / (2 eps)))/2,
where d is the signed distance into the ice: disc, lens (slab across a
channel), and their sums.

Nothing here imports the tools under test.
"""
from __future__ import annotations

import base64
import os
import struct

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
POSTPROCESS = os.path.dirname(HERE)
PROJECT = os.path.dirname(POSTPROCESS)
PY = os.path.join(PROJECT, "venv_cryosuction", "bin", "python")

# The run this project keeps locally for integration tests (stage 0c v2check,
# channel 375 x 125 um, open face x0 at h = 0.99, 22 snapshots). Tests that use
# it skip when it is absent.
REAL_RUN = os.path.expanduser(
    "~/SimulationResults/lunar_cryosuction/scratch/channel_2D_H125um_eps0.39um_deadend/"
    "2026-11-05__07.57.41_cryo_T-20_h0.990_iso_20ks_ac1e-2_v2check")


def real_run_available() -> bool:
    return os.path.isfile(os.path.join(REAL_RUN, "igasol.dat"))


# ---------------------------------------------------------------------------
# Shapes
# ---------------------------------------------------------------------------
def profile(d, eps):
    """phi for signed distance d (positive inside the ice)."""
    return 0.5 * (1.0 + np.tanh(np.asarray(d, dtype=float) / (2.0 * eps)))


def disc(X, Y, cx, cy, R, eps):
    return profile(R - np.hypot(X - cx, Y - cy), eps)


def lens(X, a, b, eps):
    """Slab a < x < b spanning the channel (a plug)."""
    return profile(np.minimum(X - a, b - X), eps)


def union(*phis):
    return np.clip(np.sum(phis, axis=0), 0.0, 1.0)


# ---------------------------------------------------------------------------
# Folder and opts
# ---------------------------------------------------------------------------
SOLVER_OPTS_DEFAULT = {"-periodic": "1", "-flag_BC_rhovfix": "0", "-humidity": "1.0",
                       "-phi_open_bc": "1", "-dim": "2"}


def write_opts(path, opts: dict, header="# synthetic"):
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for k, v in opts.items():
            fh.write(f"{k} {v}\n" if v != "" else f"{k}\n")


def make_run(root, geom="geomA", exp="expA", layout="local", solver=None, geometry=None,
             experiment=None, tag=""):
    """Create a run folder in either layout and stage three opts files.

    local: <root>/<geom>/2026-01-01__00.00.00_<exp>[_tag]/
    hpc:   <root>/<geom>__<exp>/
    """
    if layout == "local":
        run = os.path.join(root, geom, f"2026-01-01__00.00.00_{exp}" + (f"_{tag}" if tag else ""))
    elif layout == "hpc":
        run = os.path.join(root, f"{geom}__{exp}")
    else:
        raise ValueError(layout)
    os.makedirs(run, exist_ok=True)
    write_opts(os.path.join(run, "solver.opts"), {**SOLVER_OPTS_DEFAULT, **(solver or {})})
    write_opts(os.path.join(run, f"{geom}.opts"), {"-periodic": "0", **(geometry or {})})
    write_opts(os.path.join(run, f"{exp}.opts"),
               {"-flag_BC_rhovfix": "1", "-t_final": "1e4", **(experiment or {})})
    return run


# ---------------------------------------------------------------------------
# Solver text outputs (formats from src/monitoring.c)
# ---------------------------------------------------------------------------
def write_ssa(run, t, ice, rhov=None, dt=None, air=None, rho_ice=919.0, ssa=None, steps=None,
              ncol=8):
    """SSA_evo.dat: ssa ice t step dt air rhov mass, one row per step, %e."""
    t = np.asarray(t, dtype=float)
    n = t.size
    ice = np.broadcast_to(np.asarray(ice, dtype=float), (n,))
    rhov = np.zeros(n) if rhov is None else np.broadcast_to(np.asarray(rhov, dtype=float), (n,))
    dt = np.r_[np.diff(t), np.diff(t)[-1:] if n > 1 else [1.0]] if dt is None else dt
    air = np.zeros(n) if air is None else air
    ssa = np.zeros(n) if ssa is None else ssa
    steps = np.arange(n) if steps is None else steps
    with open(os.path.join(run, "SSA_evo.dat"), "w") as fh:
        for k in range(n):
            row = [ssa[k], ice[k], t[k], steps[k], dt[k], air[k], rhov[k], rho_ice * ice[k] + rhov[k]]
            row = row[:ncol]
            fh.write(" ".join(f"{int(v)}" if j == 3 else f"{v:e}" for j, v in enumerate(row)) + "\n")


def write_bflux(run, t, faces: dict, steps=None, axisym=False):
    """bflux.dat. faces: {"x0": (J, Q), ...} with J, Q arrays (or scalars), in
    the solver's column order x0 x1 y0 y1."""
    t = np.asarray(t, dtype=float)
    steps = np.arange(t.size) if steps is None else steps
    order = [f for f in ("x0", "x1", "y0", "y1", "z0", "z1") if f in faces]
    with open(os.path.join(run, "bflux.dat"), "w") as fh:
        fh.write("# outward physical flux (unscaled D_v, k); vapor kg/s, heat W; "
                 + ("axisymmetric, full 2*pi*r surface" if axisym else "per unit depth in 2D")
                 + "\n# t step" + "".join(f" J_vap_{f} Q_heat_{f}" for f in order) + "\n")
        cols = {f: [np.broadcast_to(np.asarray(v, dtype=float), t.shape) for v in faces[f]]
                for f in order}
        for k in range(t.size):
            fh.write(f"{t[k]:.10e} {int(steps[k])}"
                     + "".join(f" {cols[f][0][k]:.10e} {cols[f][1][k]:.10e}" for f in order) + "\n")


def write_resolution(run, t, vn, eps=3.9358e-7, beta0=7.9408e5, d0=1.01517e-9, phi_face=None,
                     steps=None, dt=None):
    t = np.asarray(t, dtype=float)
    vn = np.broadcast_to(np.asarray(vn, dtype=float), t.shape)
    phi_face = np.zeros(t.size) if phi_face is None else np.broadcast_to(phi_face, t.shape)
    steps = np.arange(1, t.size + 1) if steps is None else steps
    dt = np.full(t.size, 1.0) if dt is None else dt
    with open(os.path.join(run, "resolution.dat"), "w") as fh:
        fh.write("# K&P Eq. 45 check: Pi = beta0 v_n eps/d0 (<= 0.5 design, > 1 violated); "
                 f"eps {eps:.6e} m, beta0 {beta0:.6e} s/m, d0 {d0:.6e} m\n"
                 "# t step dt vn_max Pi phi_max_reservoir_face\n")
        for k in range(t.size):
            pi = beta0 * vn[k] * eps / d0
            fh.write(f"{t[k]:.10e} {int(steps[k])} {dt[k]:.6e} {vn[k]:.6e} {pi:.6e} {phi_face[k]:.6e}\n")


def write_outp(run, steps, t, tau_sub=130.49, extra=None, every_header=20):
    """outp.txt: header lines the tools grep for and the 9-column monitor table."""
    with open(os.path.join(run, "outp.txt"), "w") as fh:
        if extra:
            fh.write(f"Extra opts   : {extra}\n")
        fh.write(f"   tau_sub  =  {tau_sub:.4e} s\n   dtmax/tau_sub = 2.0002   (ok, <= 5.0)\n")
        hdr = ("    STEP |     TIME [s] |    DT [s] |    TOT_ICE |    TOT_AIR |      TEMP |  TOT_RHOV "
               "| I-A INTERF | TOTAL_MASS\n")
        for k, (s, tt) in enumerate(zip(steps, t)):
            if k % every_header == 0:
                fh.write(hdr)
            fh.write(f"  {int(s):4d} |  {tt:.5e} | 1.000e-04 |  2.156e-08 |  2.531e-08 | -9.375e-07 "
                     f"| 2.148e-11 |  1.839e-11 |  1.982e-05\n")


def write_solver_evo(run, steps, t, newton=1, krylov=30, rej=0):
    with open(os.path.join(run, "solver_evo.dat"), "w") as fh:
        fh.write("# step t dt newton_its krylov_its krylov_per_newton rejections\n")
        for s, tt in zip(steps, t):
            fh.write(f"{int(s)} {tt:e} {1.0:e} {newton} {krylov} {krylov / max(newton, 1):.2f} {rej}\n")


# ---------------------------------------------------------------------------
# Snapshots: VTS (plot_fields.py format) and igakit (solver format)
# ---------------------------------------------------------------------------
def _b64(arr):
    a = np.ascontiguousarray(arr, dtype="<f8").ravel()
    raw = struct.pack("<Q", a.nbytes) + a.tobytes()
    return base64.b64encode(raw).decode()


def write_vts(path, X, Y, fields: dict):
    """VTK XML StructuredGrid with inline base64 Float64 arrays and a UInt64
    header, the layout pplib.read_vts decodes. X, Y and fields are (ny, nx)."""
    ny, nx = X.shape
    pts = np.stack([X, Y, np.zeros_like(X)], axis=-1)
    with open(path, "w") as fh:
        fh.write('<?xml version="1.0"?>\n<VTKFile type="StructuredGrid" version="0.1" '
                 'byte_order="LittleEndian" header_type="UInt64">\n'
                 f'<StructuredGrid WholeExtent="0 {nx - 1} 0 {ny - 1} 0 0">\n'
                 f'<Piece Extent="0 {nx - 1} 0 {ny - 1} 0 0">\n<PointData>\n')
        for name, F in fields.items():
            fh.write(f'<DataArray type="Float64" Name="{name}" format="binary">{_b64(F)}</DataArray>\n')
        fh.write('</PointData>\n<Points>\n<DataArray type="Float64" NumberOfComponents="3" '
                 f'format="binary">{_b64(pts)}</DataArray>\n</Points>\n</Piece>\n'
                 '</StructuredGrid>\n</VTKFile>\n')


def write_vts_run(run, X, Y, snaps, write_pvd=True):
    """snaps: list of (step, t, phi, T, rhov) on the (ny, nx) grid X, Y."""
    vdir = os.path.join(run, "vtkOut")
    os.makedirs(vdir, exist_ok=True)
    entries = []
    for step, t, phi, T, rv in snaps:
        rel = os.path.join("vtkOut", f"solV_{step:05d}.vts")
        write_vts(os.path.join(run, rel), X, Y, {"IcePhase": phi, "Temperature": T,
                                                 "VaporDensity": rv})
        entries.append((t, rel))
    if write_pvd:
        with open(os.path.join(run, "pf.pvd"), "w") as fh:
            fh.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                     'byte_order="LittleEndian">\n  <Collection>\n')
            for t, rel in entries:
                fh.write(f'    <DataSet timestep="{t:.6e}" part="0" file="{rel}"/>\n')
            fh.write("  </Collection>\n</VTKFile>\n")


def iga_rectangle(Lx, Ly, nx, ny, p=2, mapping=None):
    """igakit NURBS of [0, Lx] x [0, Ly], nx by ny elements, degree p, C^(p-1).

    mapping(X, Y) -> (X', Y') optionally moves the control points (a mapped
    geometry, like the -geom_file meshes)."""
    from igakit.nurbs import NURBS

    def knots(n, L):
        return np.r_[[0.0] * p, np.linspace(0.0, L, n + 1), [L] * p]

    U, V = knots(nx, Lx), knots(ny, Ly)
    gu = np.array([U[i + 1:i + p + 1].mean() for i in range(nx + p)])     # Greville points
    gv = np.array([V[j + 1:j + p + 1].mean() for j in range(ny + p)])
    GX, GY = np.meshgrid(gu, gv, indexing="ij")
    if mapping is not None:
        GX, GY = mapping(GX, GY)
    C = np.zeros((nx + p, ny + p, 4))
    C[..., 0], C[..., 1], C[..., 3] = GX, GY, 1.0
    return NURBS([U, V], C)


def write_iga_run(run, nrb, snaps):
    """Write igasol.dat and sol_<step>.dat.

    snaps: list of (step, f) where f(X, Y) -> (phi, T, rhov) is evaluated at the
    control points (Greville abscissae). Linear fields are reproduced exactly;
    smooth fields to O(h^2)."""
    from igakit.io import PetIGA
    io = PetIGA()
    io.write(os.path.join(run, "igasol.dat"), nrb)
    GX, GY = nrb.points[..., 0], nrb.points[..., 1]
    for step, f in snaps:
        phi, T, rv = f(GX, GY)
        A = np.stack([np.broadcast_to(phi, GX.shape), np.broadcast_to(T, GX.shape),
                      np.broadcast_to(rv, GX.shape)], axis=-1)
        io.write_vec(os.path.join(run, f"sol_{step:05d}.dat"), A, nrb)


def grid(Lx, Ly, nx, ny):
    """Uniform (ny, nx) node grid of [0, Lx] x [0, Ly]."""
    return np.meshgrid(np.linspace(0, Lx, nx), np.linspace(0, Ly, ny))
