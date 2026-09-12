"""Shared postprocessing helpers.

Every function here previously existed as two to five copy-pasted variants
across the scripts in this directory, and every one of them had drifted:

    load_ssa        5 copies, 5 different contracts for the SAME file
    rho_vs          3 copies (and all three were wrong -- see below)
    read_vts        3 copies
    step_times      3 copies
    auto_time_unit  2 copies

Import from here instead of re-deriving. Running a script directly puts this
directory on sys.path, so ``from pplib import rho_vs`` just works.
"""

from __future__ import annotations

import base64
import os
import re
import struct
import xml.etree.ElementTree as ET

import numpy as np


# ---------------------------------------------------------------------------
# SSA_evo.dat
# ---------------------------------------------------------------------------
# Column layout written by monitoring.c. The last three are only present in
# runs produced after those quadrature integrals were added; older files stop
# at DT, and files older still stop at STEP.
SSA, ICE, TIME, STEP, DT, AIR, RHOV, MASS = range(8)

SSA_COLUMNS = """\
  0 sub_interf/eps   ice-air interface density proxy
  1 tot_ice          integrated ice volume
  2 t                [s]
  3 step
  4 dt               [s]   (NaN for legacy 4-column files)
  5 tot_air          only in runs newer than the monitoring.c quadrature change
  6 tot_rhov               "
  7 tot_mass               "
"""


def load_ssa(path: str, min_cols: int = 4):
    """Load SSA_evo.dat, returning an (N, >=5) array or None.

    `path` may be the file itself or the run directory containing it.

    A legacy 4-column file is padded with a NaN dt column so callers can index
    DT unconditionally. Rows are dropped only when one of the first four
    columns is NaN -- never on a NaN in dt or in the newer integral columns,
    which is what an earlier variant of this function got wrong: it filtered on
    `isnan(row).any()` and so discarded *every* row of a padded legacy file.

    Rows are then deduplicated by step number, keeping the last occurrence, to
    guard against repeated monitor calls after a timestep retry.
    """
    if os.path.isdir(path):
        path = os.path.join(path, "SSA_evo.dat")
    if not os.path.isfile(path):
        return None

    try:
        data = np.genfromtxt(path, dtype=float, comments="#", invalid_raise=False)
    except Exception:
        return None

    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.size == 0 or data.shape[1] < min_cols:
        return None

    if data.shape[1] == 4:
        data = np.hstack([data, np.full((len(data), 1), np.nan)])

    data = data[~np.isnan(data[:, :4]).any(axis=1)]
    if len(data) == 0:
        return None

    steps = data[:, STEP].astype(int)
    _, last = np.unique(steps[::-1], return_index=True)
    return data[np.sort(len(steps) - 1 - last)]


# ---------------------------------------------------------------------------
# Staged .opts files
# ---------------------------------------------------------------------------
def read_opts(run_dir: str) -> dict:
    """Merge every ``*.opts`` staged in a run folder into a {flag: value} dict.

    The run scripts copy solver / geometry / experiment opts into the run
    folder, so the settings a run actually used are recoverable from the folder
    alone -- no source tree, no guessing from the folder name.

    Later files win, mirroring how PETSc processes repeated ``-options_file``
    flags. That ordering is only approximated here (files are read in sorted
    name order) so do not rely on it to resolve a flag that is genuinely set
    twice to different values; every flag read by this module is set in exactly
    one of the three.

    Keys keep their leading dash. Valueless boolean flags map to "".
    """
    opts: dict[str, str] = {}
    if not os.path.isdir(run_dir):
        return opts
    for fn in sorted(os.listdir(run_dir)):
        if not fn.endswith(".opts"):
            continue
        try:
            with open(os.path.join(run_dir, fn), errors="replace") as fh:
                for line in fh:
                    line = line.split("#", 1)[0].strip()
                    if not line.startswith("-"):
                        continue
                    parts = line.split(None, 1)
                    opts[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
        except OSError:
            continue
    return opts


def opt_float(opts: dict, key: str, default=None):
    """Read a numeric flag out of read_opts(), or `default` if absent/bad."""
    try:
        return float(opts[key])
    except (KeyError, ValueError, TypeError):
        return default


def grain_radii(run_dir: str):
    """Radii [m] from the packing this run used, or None.

    Resolves ``-grains_file`` out of the staged opts. The path recorded there
    is relative to the PROJECT root, not the run folder, so the file usually
    does not travel with the results -- a copy staged next to the opts is
    preferred when present.

    Rows are ``x y r`` (or legacy ``x y z r``), with an optional leading
    two-field ``Lx Ly`` header, matching ReadGrainsFromFile in
    src/initial_conditions.c. Periodic edge images are included in the file and
    are NOT filtered out here: they are real grain area inside the domain.
    """
    opts = read_opts(run_dir)
    rel = opts.get("-grains_file", "").strip()
    if not rel:
        return None

    candidates = [os.path.join(run_dir, os.path.basename(rel)), rel]
    path = next((p for p in candidates if os.path.isfile(p)), None)
    if path is None:
        return None

    radii = []
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split()
            if len(f) == 2:          # 'Lx Ly' domain header
                continue
            try:
                vals = [float(v) for v in f]
            except ValueError:
                continue
            if len(vals) == 3:
                radii.append(vals[2])
            elif len(vals) == 4:
                radii.append(vals[3])
    return np.array(radii) if radii else None


# ---------------------------------------------------------------------------
# Material properties
# ---------------------------------------------------------------------------
# Empirical saturation-pressure coefficients, copied from RhoVS_I in
# src/material_properties.c. Keep them in step with that function: postprocess
# reporting a different rho_vs than the solver used is worse than useless.
_K = (-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918)
_PATM = 101325.0
_BB = 0.62
RHO_AIR_DEFAULT = 1.341          # -rho_air default in enceladus_main.c


def rho_vs(T_C, rho_air: float = RHO_AIR_DEFAULT):
    """Saturation vapour density over ice [kg/m^3] at T_C [deg C].

    This mirrors the solver's RhoVS_I exactly.

    The three previous copies of this function all used

        3.25e-3 * exp(-6150 / T_K)

    which is roughly 1e10 times too small -- the prefactor should have been
    ~3e7 for that Clausius-Clapeyron form, so the exponent looks like a typo.
    Because supersaturation is computed as (rhov - rho_vs) / rho_vs, the error
    did not cancel: a saturated run, which should report 0, reported ~9e9.
    Any supersaturation figure produced before 2026-08-07 is wrong.
    """
    T_K = np.asarray(T_C, dtype=float) + 273.15
    Pvs = np.exp(_K[0] / T_K + _K[1] + _K[2] * T_K + _K[3] * T_K ** 2
                 + _K[4] * T_K ** 3 + _K[5] * np.log(T_K))
    return rho_air * _BB * Pvs / (_PATM - Pvs)


def supersaturation(rhov, T_C, rho_air: float = RHO_AIR_DEFAULT):
    """(rhov - rho_vs) / rho_vs, zero where rho_vs is non-positive."""
    rvs = rho_vs(T_C, rho_air)
    return np.where(rvs > 0, (np.asarray(rhov) - rvs) / rvs, 0.0)


# ---------------------------------------------------------------------------
# VTS snapshots
# ---------------------------------------------------------------------------
def _decode(da) -> np.ndarray:
    """Decode a base64 appended-format DataArray element."""
    raw = base64.b64decode("".join(da.text.split()))
    n = struct.unpack("<Q", raw[:8])[0]
    return np.frombuffer(raw[8:8 + n], dtype=np.float64)


def read_vts(fn, want=None):
    """Read a solV_*.vts snapshot.

    Returns (fields, X, Y) where X and Y are (ny, nx) meshgrids. Pass `want`
    to keep only named point-data arrays; the default reads all of them.

    Callers that want 1-D coordinate vectors should slice: X[0, :], Y[:, 0].
    """
    root = ET.parse(fn).getroot()
    grid = root.find(".//StructuredGrid")
    ext = [int(v) for v in grid.get("WholeExtent").split()]
    nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1

    pts = _decode(root.find(".//Points/DataArray")).reshape(ny, nx, 3)

    fields = {}
    for da in root.findall(".//PointData/DataArray"):
        name = da.get("Name")
        if want is None or name in want:
            fields[name] = _decode(da).reshape(ny, nx)

    return fields, pts[:, :, 0], pts[:, :, 1]


def step_of(fn) -> int:
    """Step number from a solV_NNNNN.vts filename."""
    return int(re.search(r"solV_(\d+)\.vts", str(fn)).group(1))


def step_times(path) -> dict:
    """Map step -> time [s] from the monitor tables in outp.txt.

    `path` may be the run directory or outp.txt itself. Returns {} when the
    file is absent, so callers can fall back to SSA_evo.dat.
    """
    if os.path.isdir(path):
        path = os.path.join(path, "outp.txt")
    tmap: dict[int, float] = {}
    if not os.path.isfile(path):
        return tmap

    pat = re.compile(r"^\s+(\d+)\s+\|\s+([0-9.eE+-]+)\s+\|")
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.count("|") != 8:
                continue
            m = pat.match(line)
            if m:
                tmap[int(m.group(1))] = float(m.group(2))
    return tmap


# ---------------------------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------------------------
_TIME_UNITS = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0}


def auto_time_unit(t_max_sec: float) -> str:
    """Pick a sensible x-axis time unit for a run of this length."""
    if t_max_sec <= 600:
        return "s"
    if t_max_sec <= 7200:
        return "min"
    if t_max_sec <= 3 * 86400:
        return "h"
    return "d"


def in_time_unit(t_sec, unit: str):
    """Convert seconds to `unit` (one of s, min, h, d)."""
    return np.asarray(t_sec, dtype=float) / _TIME_UNITS[unit]


# ---------------------------------------------------------------------------
# Interface geometry from a phase field
#
# These four moved here from wedge_gt_velocity.py when contact_angle.py needed
# the same primitives. They are unchanged -- the point is that there is ONE
# implementation, so a fix to the crossing convention reaches every measurement
# that depends on it.
# ---------------------------------------------------------------------------


def crossings(x, y, level=0.5):
    """Linear-interpolation crossings of y=level in y(x).  Never a grid node."""
    s = np.flatnonzero((y[:-1] - level) * (y[1:] - level) < 0.0)
    if s.size == 0:
        return np.empty(0)
    f = (level - y[s]) / (y[s + 1] - y[s])
    return x[s] + f * (x[s + 1] - x[s])


def refine_tanh(x, phi, x0, eps, half_width=6.0, lo=0.02, hi=0.98):
    """Interface position from the whole diffuse band, not two samples.

    The equilibrium profile of this model's half-normalised well is
    phi = 0.5*(1 + tanh(x/(2*eps))), so atanh(2*phi - 1) is EXACTLY linear in x
    across the interface. Fitting that line over the band and taking its zero
    uses ~13 samples where the 0.5 crossing uses 2, which matters here: v_n is
    a time derivative of this position, so a bias that repeats with the sample
    grid shows up as a periodic ripple at the period of one cell crossing.
    Measured on rhov_eq_eq, this cuts that ripple 3x on control-net data and
    6x on true-NURBS data. Falls back to `x0` if the band is too thin to fit.
    """
    m = (np.abs(x - x0) < half_width * eps) & (phi > lo) & (phi < hi)
    if m.sum() < 5:
        return x0
    z = np.arctanh(np.clip(2.0 * phi[m] - 1.0, -1.0 + 1e-12, 1.0 - 1e-12))
    slope, intercept = np.polyfit(x[m] - x0, z, 1)
    if not np.isfinite(slope) or slope == 0.0:
        return x0
    return x0 - intercept / slope


def contour_points(x, Y, phi, level=0.5):
    """phi=level crossings row by row, as physical (x, y) pairs.

    Returns (xl, yl, xr, yr) for the rows that have exactly two crossings --
    i.e. the left and right menisci sampled along their whole arcs.
    """
    d = phi - level
    hit = (d[:, :-1] * d[:, 1:]) < 0.0
    jj, ii = np.nonzero(hit)
    if jj.size == 0:
        return (np.empty(0),) * 4
    f = (level - phi[jj, ii]) / (phi[jj, ii + 1] - phi[jj, ii])
    px = x[ii] + f * (x[ii + 1] - x[ii])
    py = Y[jj, ii] + f * (Y[jj, ii + 1] - Y[jj, ii])

    n_per_row = np.bincount(jj, minlength=phi.shape[0])
    keep = n_per_row[jj] == 2
    if not keep.any():
        return (np.empty(0),) * 4
    jj, px, py = jj[keep], px[keep], py[keep]
    order = np.lexsort((px, jj))            # per row: left crossing first
    px, py = px[order], py[order]
    return px[0::2], py[0::2], px[1::2], py[1::2]


def circle_radius(px, py):
    """Least-squares circle through (px, py).  Returns (xc, yc, r) or None.

    Linear formulation 2*a*x + 2*b*y + c = x^2 + y^2, so no iteration and no
    initial guess.  Used only to VALIDATE the analytic apex radius -- second
    differences of the contour would be the noise-amplifying alternative that
    neck_width.py:18-23 rejects for exactly this geometry.
    """
    if px.size < 5:
        return None
    A = np.column_stack((2.0 * px, 2.0 * py, np.ones(px.size)))
    b = px ** 2 + py ** 2
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    xc, yc = sol[0], sol[1]
    disc = sol[2] + xc ** 2 + yc ** 2
    if not np.isfinite(disc) or disc <= 0.0:
        return None
    return xc, yc, float(np.sqrt(disc))
