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
    data = data[np.sort(len(steps) - 1 - last)]
    # Step order, so a hand-concatenated restart that rewinds (steps 0..5,
    # then 3..4) never runs time backwards.
    return data[np.argsort(data[:, STEP], kind="stable")]


# ---------------------------------------------------------------------------
# Staged .opts files
# ---------------------------------------------------------------------------
def _parse_opts_file(path: str) -> dict:
    """{flag: value} from one .opts file. Comments (#) are stripped; a
    valueless flag maps to "". Within one file the last setting wins."""
    out: dict[str, str] = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line.startswith("-"):
                continue
            parts = line.split(None, 1)
            out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


_TS_PREFIX = re.compile(r"^\d{4}-\d\d-\d\d__\d\d\.\d\d\.\d\d_")


def opts_order(run_dir: str) -> list[str]:
    """The staged *.opts of a run, in the order the solver read them.

    run_cryosuction.sh passes solver.opts, then the geometry file, then the
    experiment file, and PETSc lets the later file win. Two folder layouts are
    in use:

        local  <scratch>/<geometry>/<YYYY-MM-DD__HH.MM.SS>_<experiment>[_tag]/
        HPC    <batch>/<geometry>__<experiment>/

    The geometry file is the one whose stem is the parent folder name (local)
    or the part before "__" (HPC); the experiment file is the one whose stem
    begins the rest of the folder name. If neither rule matches, a file that
    sets -t_final is taken as the experiment (experiment files carry the run
    length; geometry files never do). Anything still unplaced sits between
    geometry and experiment, in name order.
    """
    if not os.path.isdir(run_dir):
        return []
    files = sorted(f for f in os.listdir(run_dir) if f.endswith(".opts"))
    name = os.path.basename(os.path.normpath(run_dir))
    parent = os.path.basename(os.path.dirname(os.path.normpath(run_dir)))
    if _TS_PREFIX.match(name):
        geom, rest = parent, _TS_PREFIX.sub("", name)
    elif "__" in name:
        geom, _, rest = name.partition("__")
    else:
        geom, rest = parent, name

    def rank(fn: str) -> tuple:
        stem = fn[:-len(".opts")]
        if stem == "solver":
            return (0, fn)
        if stem == geom:
            return (1, fn)
        if rest == stem or rest.startswith(stem + "_"):
            return (3, fn)
        try:
            if "-t_final" in _parse_opts_file(os.path.join(run_dir, fn)):
                return (3, fn)
        except OSError:
            pass
        return (2, fn)

    return [os.path.join(run_dir, f) for f in sorted(files, key=rank)]


def extra_opts(run_dir: str) -> dict:
    """Flags passed after `--` to the run script, recovered from outp.txt.

    run_cryosuction.sh echoes them as "Extra opts   : -a 1 -b 2" and passes
    them after the three -options_file flags, so they beat every file.
    Returns {} when the line is absent (no extra flags, or no outp.txt).
    """
    path = os.path.join(run_dir, "outp.txt")
    if not os.path.isfile(path):
        return {}
    with open(path, errors="replace") as fh:
        for i, line in enumerate(fh):
            if i > 400:              # the line is in the header
                break
            m = re.match(r"^Extra opts\s*:\s*(.*)$", line.strip())
            if m:
                toks, out = m.group(1).split(), {}
                k = 0
                while k < len(toks):
                    if toks[k].startswith("-") and not _is_number(toks[k]):
                        has_val = k + 1 < len(toks) and (not toks[k + 1].startswith("-")
                                                         or _is_number(toks[k + 1]))
                        out[toks[k]] = toks[k + 1] if has_val else ""
                        k += 2 if has_val else 1
                    else:
                        k += 1
                return out
    return {}


def _is_number(tok: str) -> bool:
    try:
        float(tok)
        return True
    except ValueError:
        return False


def read_opts(run_dir: str) -> dict:
    """Merge the staged ``*.opts`` of a run folder into a {flag: value} dict,
    with the precedence the solver used: solver.opts < geometry < experiment
    < flags passed after `--` (see opts_order and extra_opts).

    The run scripts copy solver / geometry / experiment opts into the run
    folder, so the settings a run actually used are recoverable from the folder
    alone -- no source tree, no guessing from the folder name.

    Before 2026-11-06 this read the files in sorted name order, which put
    solver.opts LAST: a run whose experiment set -flag_BC_rhovfix 1 was
    reported as 0, and -periodic came back as solver.opts' 1. Only flags set
    in more than one file were affected.

    Keys keep their leading dash. Valueless boolean flags map to "".
    """
    opts: dict[str, str] = {}
    for path in opts_order(run_dir):
        try:
            opts.update(_parse_opts_file(path))
        except OSError:
            continue
    opts.update(extra_opts(run_dir))
    return opts


def opt_bool(opts: dict, key: str, default: bool = False) -> bool:
    """PETSc boolean flag: absent -> default; "" (bare flag), 1, true, yes, on -> True."""
    if key not in opts:
        return default
    v = str(opts[key]).strip().lower()
    if v in ("", "1", "true", "yes", "on", "t", "y"):
        return True
    if v in ("0", "false", "no", "off", "f", "n"):
        return False
    return default


def face_list(opts: dict, key: str) -> list[tuple[int, int]]:
    """Parse a face list such as -wall_faces "y0,y1" into [(axis, side), ...]."""
    raw = str(opts.get(key, "")).strip().strip('"').strip("'")
    out = []
    for tok in re.split(r"[,\s]+", raw):
        m = re.fullmatch(r"([xyz])([01])", tok.strip().lower())
        if m:
            out.append(("xyz".index(m.group(1)), int(m.group(2))))
    return out


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

    # A copy staged in the run folder first, then the path as written
    # (relative to the project root, per the run scripts), resolved against
    # the project root this file lives in and against the working directory.
    project = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates = [os.path.join(run_dir, os.path.basename(rel)), os.path.join(project, rel), rel]
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
    """Saturation vapor density over ice [kg/m^3] at T_C [deg C].

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
    """Linear-interpolation crossings of y=level in y(x), in x order.

    A sample exactly at the level is a crossing (returned once, at that node)
    when its two neighbors lie on opposite sides; a profile that only touches
    the level is not a crossing.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    d = y - level
    s = np.flatnonzero(d[:-1] * d[1:] < 0.0)
    f = (level - y[s]) / (y[s + 1] - y[s])
    out = x[s] + f * (x[s + 1] - x[s])
    if d.size >= 3:
        k = np.flatnonzero(d[1:-1] == 0.0) + 1
        k = k[d[k - 1] * d[k + 1] < 0.0]
        if k.size:
            out = np.sort(np.r_[out, x[k]])
    return out


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
    root = x0 - intercept / slope
    # A band too thin or not monotone gives a nearly flat line whose root can
    # land far outside the fitted window; then the crossing is the better guess.
    if abs(root - x0) > half_width * eps:
        return x0
    return root


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
        sol, _, rank, sv = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    # Collinear points: A = [2x, 2y, 1] is rank 2 (exactly or to round-off)
    # and there is no finite circle.
    if rank < 3 or sv[-1] <= 1e-9 * sv[0]:
        return None
    xc, yc = sol[0], sol[1]
    disc = sol[2] + xc ** 2 + yc ** 2
    if not np.isfinite(disc) or disc <= 0.0:
        return None
    return xc, yc, float(np.sqrt(disc))


# ---------------------------------------------------------------------------
# Transport laws, copied from src/material_properties.c and the defaults in
# src/cryosuction_main.c (lines 160-174). Keep them in step with the solver.
# ---------------------------------------------------------------------------
RHO_ICE = 919.0          # user.rho_ice
DIF_VAP0 = 2.178e-5      # user.dif_vap, m^2/s at 273.15 K
K_ICE, K_AIR = 2.29, 0.02  # user.thcond_ice, user.thcond_air, W/(m K)
AIR_LIM = 1.0e-6         # user.air_lim, floor on the air fraction in the vapor flux


def dif_vap(T_C, d0: float = DIF_VAP0):
    """Vapor diffusivity D_v(T) = D0 (T/273.15)^1.81 (VaporDiffus)."""
    return d0 * ((np.asarray(T_C, dtype=float) + 273.15) / 273.15) ** 1.81


def thermal_cond(phi, k_ice: float = K_ICE, k_air: float = K_AIR):
    """k(phi) = phi k_ice + (1 - phi) k_air with both fractions clipped at 0 (ThermalCond)."""
    phi = np.asarray(phi, dtype=float)
    return np.clip(phi, 0.0, None) * k_ice + np.clip(1.0 - phi, 0.0, None) * k_air


def vapor_air_fraction(phi, air_lim: float = AIR_LIM):
    """phi_aef = max(1 - clip(phi, 0, 1), air_lim), as in BoundaryFluxMonitor."""
    return np.maximum(1.0 - np.clip(np.asarray(phi, dtype=float), 0.0, 1.0), air_lim)


# ---------------------------------------------------------------------------
# Snapshot series and quadrature on a sampled (possibly curvilinear) grid
# ---------------------------------------------------------------------------
def read_pvd(run_dir: str) -> list[tuple[float, str]]:
    """[(t, absolute path to the .vts), ...] in time order, from pf.pvd.

    pf.pvd is written by plot_fields.py with the snapshot times. Falls back to
    vtkOut/solV_*.vts with times from step_times(outp.txt) when pf.pvd is
    missing; a snapshot with no known time is dropped.
    """
    pvd = os.path.join(run_dir, "pf.pvd")
    out = []
    if os.path.isfile(pvd):
        root = ET.parse(pvd).getroot()
        for ds in root.iter("DataSet"):
            f = os.path.join(run_dir, ds.get("file"))
            if os.path.isfile(f):
                out.append((float(ds.get("timestep")), f))
    else:
        tmap = step_times(run_dir)
        vdir = os.path.join(run_dir, "vtkOut")
        if os.path.isdir(vdir):
            for fn in sorted(os.listdir(vdir)):
                if re.fullmatch(r"solV_\d+\.vts", fn) and step_of(fn) in tmap:
                    out.append((tmap[step_of(fn)], os.path.join(vdir, fn)))
    out.sort(key=lambda p: p[0])
    # Keep the last file for a repeated time (a restart re-writes its first snapshot).
    dedup = {}
    for t, f in out:
        dedup[t] = f
    return sorted(dedup.items())


def node_weights(X, Y):
    """Quadrature weight of every node of a structured (ny, nx) grid [m^2].

    Each quadrilateral cell's area (shoelace formula, exact for straight-sided
    cells) is shared equally among its four corners, so sum(w * f) is the
    trapezoid rule on a uniform grid and its bilinear analogue on a mapped one.
    Integrals of a sampled field differ from the solver's Gauss quadrature by
    the sampling error, not by a definition.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    x0, x1, x2, x3 = X[:-1, :-1], X[:-1, 1:], X[1:, 1:], X[1:, :-1]
    y0, y1, y2, y3 = Y[:-1, :-1], Y[:-1, 1:], Y[1:, 1:], Y[1:, :-1]
    area = 0.5 * np.abs((x0 * y1 - x1 * y0) + (x1 * y2 - x2 * y1)
                        + (x2 * y3 - x3 * y2) + (x3 * y0 - x0 * y3))
    w = np.zeros_like(X)
    q = 0.25 * area
    w[:-1, :-1] += q
    w[:-1, 1:] += q
    w[1:, 1:] += q
    w[1:, :-1] += q
    return w


def grid_gradient(F, X, Y):
    """(dF/dx, dF/dy) at every node of a structured, possibly curvilinear grid.

    Second-order differences in index space (np.gradient), mapped to physical
    space with the inverse Jacobian of (X, Y) with respect to (i, j). Exact for
    a linear field on any grid whose mapping is itself smooth.
    """
    F = np.asarray(F, dtype=float)
    Fj, Fi = np.gradient(F)            # axis 0 = j (rows), axis 1 = i (columns)
    Xj, Xi = np.gradient(np.asarray(X, dtype=float))
    Yj, Yi = np.gradient(np.asarray(Y, dtype=float))
    det = Xi * Yj - Xj * Yi
    with np.errstate(divide="ignore", invalid="ignore"):
        dFdx = (Fi * Yj - Fj * Yi) / det
        dFdy = (Fj * Xi - Fi * Xj) / det
    return dFdx, dFdy


def is_rectilinear(X, Y, rtol: float = 1e-9) -> bool:
    """True when X depends only on the column and Y only on the row."""
    X = np.asarray(X)
    Y = np.asarray(Y)
    sx = max(np.ptp(X), 1e-300)
    sy = max(np.ptp(Y), 1e-300)
    return bool(np.ptp(X, axis=0).max() <= rtol * sx and np.ptp(Y, axis=1).max() <= rtol * sy)


class GridSampler:
    """Interpolate fields of one structured grid at arbitrary points.

    Rectilinear grids use scipy's RegularGridInterpolator (bilinear, exact for
    bilinear data). Mapped grids use a Delaunay triangulation of the nodes,
    built once (the mesh never moves) and reused for every field and snapshot.
    Points outside the domain return NaN.
    """

    def __init__(self, X, Y):
        self.X = np.asarray(X, dtype=float)
        self.Y = np.asarray(Y, dtype=float)
        self.rect = is_rectilinear(self.X, self.Y)
        if not self.rect:
            from scipy.spatial import Delaunay
            self.tri = Delaunay(np.column_stack([self.X.ravel(), self.Y.ravel()]))

    def __call__(self, F, px, py):
        px = np.asarray(px, dtype=float)
        py = np.asarray(py, dtype=float)
        F = np.asarray(F, dtype=float)
        if self.rect:
            from scipy.interpolate import RegularGridInterpolator
            f = RegularGridInterpolator((self.Y[:, 0], self.X[0, :]), F,
                                        bounds_error=False, fill_value=np.nan)
            return f(np.column_stack([py, px]))
        from scipy.interpolate import LinearNDInterpolator
        f = LinearNDInterpolator(self.tri, F.ravel(), fill_value=np.nan)
        return f(px, py)


# ---------------------------------------------------------------------------
# Vertical display (plan tool T7): the pore axis x points DOWN the page, so
# x = 0 (the face open to the surface in the cryosuction set-ups) is at the
# top and depth increases downward, as in the planning sketches.
# ---------------------------------------------------------------------------
def to_display(X, Y, vertical: bool = True):
    """Plot coordinates (h, v) for model coordinates (x, y).

    vertical=True maps x -> depth on the vertical axis (use with
    style_vertical, which inverts it) and y -> the horizontal axis.
    vertical=False returns (X, Y) unchanged.
    """
    return (Y, X) if vertical else (X, Y)


def style_vertical(ax, vertical: bool = True, unit: str = "µm"):
    """Label and orient an axes for to_display coordinates."""
    if vertical:
        if not ax.yaxis_inverted():
            ax.invert_yaxis()
        ax.set_xlabel(f"across the pore y [{unit}]")
        ax.set_ylabel(f"depth x [{unit}]")
    else:
        ax.set_xlabel(f"x [{unit}]")
        ax.set_ylabel(f"y [{unit}]")
    ax.set_aspect("equal")


# ---------------------------------------------------------------------------
# Snapshot loader for the measurement tools (ice_bodies, water_budget,
# section_flux). Two sources:
#
#   "nurbs"  sol_*.dat + igasol.dat, the spline EVALUATED on a uniform grid of
#            n_per_elem points per element (igakit). These are true field
#            values, so integrals and gradients carry only the sampling error.
#   "vts"    vtkOut/solV_*.vts. plot_fields.py writes the CONTROL-POINT values
#            there, not field values: a quadratic spline's control net is a
#            smoothed copy of the field (the 1D step 0->1 is spread over about
#            one more element) and its nodes sit at the Greville points, which
#            are unevenly spaced next to the boundary. Fine for pictures;
#            measurements from it carry an O(h) smoothing bias.
#
# "auto" takes nurbs when igakit and igasol.dat are present.
# ---------------------------------------------------------------------------
class Snapshot:
    """One saved state: step, t [s], X, Y (ny, nx) [m] and fields phi, T, rhov."""

    __slots__ = ("step", "t", "X", "Y", "phi", "T", "rhov", "source", "path")

    def __init__(self, step, t, X, Y, phi, T, rhov, source, path):
        self.step, self.t, self.X, self.Y = step, t, X, Y
        self.phi, self.T, self.rhov = phi, T, rhov
        self.source, self.path = source, path


def snapshot_times(run_dir: str) -> dict:
    """step -> t [s], from outp.txt, filled in from SSA_evo.dat when missing."""
    tmap = step_times(run_dir)
    ssa = load_ssa(run_dir)
    if ssa is not None:
        for row in ssa:
            tmap.setdefault(int(row[STEP]), float(row[TIME]))
    return tmap


def _sol_files(run_dir):
    return sorted(f for f in os.listdir(run_dir) if re.fullmatch(r"sol_\d+\.dat", f))


def iter_snapshots(run_dirs, source: str = "auto", n_per_elem: int = 2, every: int = 1):
    """Yield Snapshot objects in time order over one run or a restart chain.

    run_dirs: a folder or a list of folders (legs of a restart chain). A step
    present in two legs is taken from the later leg. Snapshots whose time is
    unknown (no outp.txt row and no SSA_evo.dat row) are skipped with a note on
    stderr, because a velocity built on a guessed time is wrong silently.
    every: keep every n-th snapshot (1 = all).
    """
    import sys
    if isinstance(run_dirs, (str, os.PathLike)):
        run_dirs = [run_dirs]
    items = {}
    for rd in run_dirs:
        rd = str(rd)
        tmap = snapshot_times(rd)
        use_nurbs = False
        if source in ("auto", "nurbs"):
            try:
                import igakit  # noqa: F401
                use_nurbs = os.path.isfile(os.path.join(rd, "igasol.dat")) and bool(_sol_files(rd))
            except ImportError:
                use_nurbs = False
            if source == "nurbs" and not use_nurbs:
                raise FileNotFoundError(f"{rd}: no igasol.dat + sol_*.dat (or igakit missing)")
        if use_nurbs:
            for f in _sol_files(rd):
                st = int(re.search(r"sol_(\d+)", f).group(1))
                items[st] = (tmap.get(st), "nurbs", os.path.join(rd, f), rd)
        else:
            vdir = os.path.join(rd, "vtkOut")
            if not os.path.isdir(vdir):
                raise FileNotFoundError(f"{rd}: no sol_*.dat/igasol.dat and no vtkOut/")
            for f in sorted(os.listdir(vdir)):
                if re.fullmatch(r"solV_\d+\.vts", f):
                    st = step_of(f)
                    items[st] = (tmap.get(st), "vts", os.path.join(vdir, f), rd)
    rows = []
    for st, (t, src, path, rd) in items.items():
        if t is None:
            print(f"pplib.iter_snapshots: no time for step {st} ({path}); skipped", file=sys.stderr)
            continue
        rows.append((t, st, src, path, rd))
    rows.sort()
    rows = rows[::max(1, int(every))] if every > 1 else rows
    cache = {}
    for t, st, src, path, rd in rows:
        if src == "nurbs":
            if rd not in cache:
                from igakit.io import PetIGA
                io = PetIGA()
                nrb = io.read(os.path.join(rd, "igasol.dat"))
                lim = [(float(nrb.knots[i][nrb.degree[i]]), float(nrb.knots[i][-nrb.degree[i] - 1]))
                       for i in range(2)]
                nel = [len(np.unique(nrb.knots[i])) - 1 for i in range(2)]
                u = np.linspace(*lim[0], n_per_elem * nel[0] + 1)
                v = np.linspace(*lim[1], n_per_elem * nel[1] + 1)
                cache[rd] = (io, nrb, u, v)
            io, nrb, u, v = cache[rd]
            C, F = nrb(u, v, fields=io.read_vec(path, nrb))
            X, Y = C[:, :, 0].T.copy(), C[:, :, 1].T.copy()
            phi, T, rv = (F[:, :, k].T.copy() for k in range(3))
        else:
            fl, X, Y = read_vts(path)
            phi, T, rv = fl["IcePhase"], fl["Temperature"], fl["VaporDensity"]
        yield Snapshot(st, t, X, Y, phi, T, rv, src, path)
