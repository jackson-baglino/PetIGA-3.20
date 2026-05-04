#!/usr/bin/env python3
"""
plot1D_profiles.py  —  Visualize 1D permafrost simulation field profiles.

Reads sol_*.dat files from a PetIGA 1D run and produces:

  Per-step phase figures (always)
    Single-panel figure per snapshot showing φ_i, φ_s, φ_a = 1−φ_i−φ_s
    plus shading where partition-of-unity is violated (air < 0).
    Output: phase_step_NNNNN.png in --out-dir/phases/.

  Per-step thermal figures (always)
    Single-panel figure per snapshot with T(x) on the left y-axis and
    ρ_v(x) on the right y-axis (twin axes).
    Output: thermal_step_NNNNN.png in --out-dir/thermal_steps/.

  First/last comparison (always)
    Single figure with φ_i, φ_s, φ_a for the first and last snapshots
    overlaid using solid (first) and dashed (last) lines.
    Output: phase_first_last.png in --out-dir.

  Thermal overlay (--thermal)
    Single-panel twin-axis figure with all snapshots overlaid:
    T(x) on the left axis, ρ_v(x) on the right axis, color-coded by time.
    Output: thermal_overlay.png

  Derived quantities (--derived)
    Time-series of ice volume fraction, interface position, slab width.

Usage
-----
  # Default: per-step phase PNGs + thermal PNGs + first/last comparison
  python plot1D_profiles.py --dir /path/to/run

  # Also produce the thermal overlay (all steps overlaid on one plot)
  python plot1D_profiles.py --dir /path/to/run --thermal

  # Derived scalar quantities only
  python plot1D_profiles.py --dir . --derived --save derived.png

  # First and last snapshot only
  python plot1D_profiles.py --dir . --first-last
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")          # safe for headless / HPC environments
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Optional cmocean — fall back gracefully if not installed
# ---------------------------------------------------------------------------
try:
    import cmocean
    _CMOCEAN = True
except ImportError:
    _CMOCEAN = False

try:
    from igakit.io import PetIGA
except ImportError:
    sys.exit(
        "ERROR: igakit is not installed.\n"
        "Install it with:  pip install igakit"
    )

# ---------------------------------------------------------------------------
# Phase field display colours (ColorBrewer-safe, high contrast)
#   ice  → steel blue  (cold / crystalline)
#   sed  → warm brown  (mineral / earth)
#   air  → muted green (pore space)
# ---------------------------------------------------------------------------
PHASE_COLORS = {
    "ice": "#2166ac",
    "sed": "#8c510a",
    "air": "#4dac26",
}

PHASE_LABELS = {
    "ice": r"$\phi_i$  (ice)",
    "sed": r"$\phi_s$  (sediment)",
    "air": r"$\phi_a$  (air)",
}

_COLOR_T    = "#d62728"   # red for temperature
_COLOR_RHOV = "#1f77b4"   # blue for vapor density


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_geometry(iga_path: str):
    if not os.path.isfile(iga_path):
        sys.exit(f"ERROR: IGA geometry file not found: {iga_path}")
    return PetIGA().read(iga_path)


def get_x_coords(nrb) -> np.ndarray:
    """Physical x-coordinates of IGA control points (1D mesh)."""
    return nrb.control[:, 0]


def step_number(path: str) -> int:
    """Extract integer step from sol_NNNNN.dat filename."""
    base   = os.path.splitext(os.path.basename(path))[0]
    digits = base.lstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")
    return int(digits) if digits else 0


def load_times(run_dir: str) -> dict:
    """Return {step: time_hours} from SSA_evo.dat. Empty dict if unavailable."""
    path   = os.path.join(run_dir, "SSA_evo.dat")
    times  = {}
    if not os.path.isfile(path):
        return times
    try:
        data = np.loadtxt(path)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        for row in data:
            if row.shape[0] >= 4:
                times[int(row[3])] = float(row[2]) / 3600.0   # s → h
    except Exception:
        pass
    return times


# ---------------------------------------------------------------------------
# Saturation vapor density (matches material_properties.c: RhoVS_I)
# ---------------------------------------------------------------------------

def rho_vs(T_C: np.ndarray) -> np.ndarray:
    """ρ_vs over ice [kg/m³] at temperature T_C [°C]."""
    return 3.25e-3 * np.exp(-6150.0 / (T_C + 273.15))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fmt_time(t_h: float) -> str:
    """Human-readable time from hours."""
    if t_h < 1.0 / 60.0:
        return f"{t_h * 3600:.1f} s"
    elif t_h < 1.0:
        return f"{t_h * 60:.1f} min"
    else:
        return f"{t_h:.2f} h"


def _find_crossings(x: np.ndarray, y: np.ndarray, level: float = 0.5) -> np.ndarray:
    """Return x positions where y crosses `level` (linear interpolation)."""
    crossings = []
    for i in range(len(y) - 1):
        y0, y1 = y[i], y[i + 1]
        if (y0 - level) * (y1 - level) < 0:
            t = (level - y0) / (y1 - y0)
            crossings.append(x[i] + t * (x[i + 1] - x[i]))
    return np.array(crossings)


def _collect_snapshots(run_dir: str, iga_file: str = "igasol.dat",
                        max_steps: int = None):
    """
    Load geometry, x-coords, all sol_*.dat files.
    Returns (x_mm, sed_list, ice_list, tem_list, rhov_list, step_list, time_h_list).
    """
    nrb       = load_geometry(os.path.join(run_dir, iga_file))
    x         = get_x_coords(nrb)
    x_mm      = x * 1e3
    times     = load_times(run_dir)

    sol_files = sorted(
        glob.glob(os.path.join(run_dir, "sol_*.dat")), key=step_number
    )
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")
    if max_steps is not None:
        sol_files = sol_files[:max_steps]

    ice_list  = []
    sed_list  = []
    tem_list  = []
    rhov_list = []
    step_list = []
    time_list = []

    for sf in sol_files:
        try:
            sol = PetIGA().read_vec(sf, nrb)
        except Exception as e:
            print(f"  WARNING: skipping {sf}: {e}")
            continue
        if sol.ndim < 2 or sol.shape[1] < 4:
            print(f"  WARNING: unexpected shape in {sf} — skipping.")
            continue

        step = step_number(sf)
        ice_list.append(sol[:, 0])
        tem_list.append(sol[:, 1])
        rhov_list.append(sol[:, 2])
        sed_list.append(sol[:, 3])
        step_list.append(step)
        time_list.append(times.get(step, None))

    if not ice_list:
        sys.exit("No valid solution files could be read.")

    print(f"Loaded {len(ice_list)} snapshots from '{run_dir}'")
    if not _CMOCEAN:
        print("  NOTE: cmocean not installed — using fallback colormaps.")

    return x_mm, sed_list, ice_list, tem_list, rhov_list, step_list, time_list


# ---------------------------------------------------------------------------
# Phase figure — all three phases on one panel
# ---------------------------------------------------------------------------

def _make_phase_fig(x_mm, phi_i, phi_s, phi_a, step, t_h):
    """Single-panel figure with φ_i, φ_s, φ_a overlaid.

    phi_a is RAW (= 1 − φ_i − φ_s, no clipping). Regions where phi_a < 0
    are shaded red to flag partition-of-unity violations.
    """
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    t_str = f"  (t = {_fmt_time(t_h)})" if t_h is not None else ""
    ax.set_title(f"Phase fields — step {step:d}{t_str}", fontsize=13)

    ax.plot(x_mm, phi_i, color=PHASE_COLORS["ice"], lw=2.0, label=PHASE_LABELS["ice"])
    ax.plot(x_mm, phi_s, color=PHASE_COLORS["sed"], lw=2.0, label=PHASE_LABELS["sed"])
    ax.plot(x_mm, phi_a, color=PHASE_COLORS["air"], lw=2.0,
            label=r"$\phi_a = 1 - \phi_i - \phi_s$  (raw)")

    violation = phi_a < 0.0
    if np.any(violation):
        ax.fill_between(x_mm, phi_a, 0.0, where=violation,
                        color="red", alpha=0.25, label="constraint violation  (air < 0)")

    ax.axhline(0.0, color="gray", lw=0.8, ls=":")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")

    ax.set_ylabel("Volume fraction", fontsize=12)
    ax.set_xlabel("x  [mm]", fontsize=12)
    y_lo = min(-0.05, np.min(phi_a) - 0.02)
    ax.set_ylim(y_lo, 1.1)
    ax.legend(fontsize=11, loc="best")
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Thermal figure — T and ρ_v on twin y-axes, single panel
# ---------------------------------------------------------------------------

def _make_thermal_step_fig(x_mm, tem, rhov, step, t_h):
    """Single-panel figure with T(x) on the left y-axis and ρ_v(x) on the right."""
    fig, ax_T = plt.subplots(1, 1, figsize=(8, 5))
    ax_rho = ax_T.twinx()

    t_str = f"  (t = {_fmt_time(t_h)})" if t_h is not None else ""
    ax_T.set_title(f"Thermal fields — step {step:d}{t_str}", fontsize=13)

    ln_T,   = ax_T.plot(x_mm, tem,  color=_COLOR_T,    lw=2.0,
                        label="Temperature  [°C]")
    ln_rho, = ax_rho.plot(x_mm, rhov, color=_COLOR_RHOV, lw=2.0,
                           label=r"Vapor density  [kg m$^{-3}$]")

    ax_T.set_xlabel("x  [mm]", fontsize=12)
    ax_T.set_ylabel("Temperature  [°C]", fontsize=12, color=_COLOR_T)
    ax_rho.set_ylabel(r"Vapor density  [kg m$^{-3}$]", fontsize=12, color=_COLOR_RHOV)
    ax_T.tick_params(axis="y", labelcolor=_COLOR_T, labelsize=10)
    ax_rho.tick_params(axis="y", labelcolor=_COLOR_RHOV, labelsize=10)
    ax_T.tick_params(axis="x", labelsize=10)

    # Combined legend
    ax_T.legend([ln_T, ln_rho],
                [ln_T.get_label(), ln_rho.get_label()],
                fontsize=11, loc="best")

    ax_T.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Per-step plot functions
# ---------------------------------------------------------------------------

def plot_phase_steps(run_dir: str, out_dir: str = None,
                     iga_file: str = "igasol.dat",
                     max_steps: int = None) -> list:
    """Save one PNG per snapshot showing φ_i, φ_s, φ_a on a single panel."""
    if out_dir is None:
        out_dir = run_dir
    phases_dir = os.path.join(out_dir, "phases")
    os.makedirs(phases_dir, exist_ok=True)

    x_mm, sed_list, ice_list, _, _, step_list, time_list = _collect_snapshots(
        run_dir, iga_file, max_steps
    )

    saved = []
    for phi_i, phi_s, step, t_h in zip(ice_list, sed_list, step_list, time_list):
        phi_a = 1.0 - phi_i - phi_s
        fig   = _make_phase_fig(x_mm, phi_i, phi_s, phi_a, step, t_h)
        path  = os.path.join(phases_dir, f"phase_step_{step:05d}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)
        print(f"  Saved: {os.path.relpath(path)}")

    return saved


def plot_thermal_steps(run_dir: str, out_dir: str = None,
                       iga_file: str = "igasol.dat",
                       max_steps: int = None) -> list:
    """Save one PNG per snapshot with T(x) and ρ_v(x) on twin y-axes."""
    if out_dir is None:
        out_dir = run_dir
    thermal_dir = os.path.join(out_dir, "thermal_steps")
    os.makedirs(thermal_dir, exist_ok=True)

    x_mm, _, _, tem_list, rhov_list, step_list, time_list = _collect_snapshots(
        run_dir, iga_file, max_steps
    )

    saved = []
    for tem, rhov, step, t_h in zip(tem_list, rhov_list, step_list, time_list):
        fig  = _make_thermal_step_fig(x_mm, tem, rhov, step, t_h)
        path = os.path.join(thermal_dir, f"thermal_step_{step:05d}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)
        print(f"  Saved: {os.path.relpath(path)}")

    return saved


# ---------------------------------------------------------------------------
# First / last snapshot comparison
# ---------------------------------------------------------------------------

def plot_first_last(run_dir: str, save_path: str = None,
                    iga_file: str = "igasol.dat"):
    """
    Overlay φ_i, φ_s, φ_a for the first (solid) and last (dashed) snapshots
    on a single panel.
    """
    x_mm, sed_list, ice_list, _, _, step_list, time_list = _collect_snapshots(
        run_dir, iga_file
    )

    if len(ice_list) < 2:
        print("WARNING: fewer than 2 snapshots found — nothing to compare.")
        return

    fig, ax = plt.subplots(1, 1, figsize=(9, 5))

    for idx, ls, label in zip([0, -1], ["-", "--"], ["first", "last"]):
        phi_i = ice_list[idx]
        phi_s = sed_list[idx]
        phi_a = 1.0 - phi_i - phi_s
        step  = step_list[idx]
        t_h   = time_list[idx]
        t_str = f", t = {_fmt_time(t_h)}" if t_h is not None else f", step {step}"

        ax.plot(x_mm, phi_i, color=PHASE_COLORS["ice"], lw=2.0, ls=ls,
                label=rf"$\phi_i$ ({label}{t_str})")
        ax.plot(x_mm, phi_s, color=PHASE_COLORS["sed"], lw=2.0, ls=ls,
                label=rf"$\phi_s$ ({label}{t_str})")
        ax.plot(x_mm, phi_a, color=PHASE_COLORS["air"], lw=2.0, ls=ls,
                label=rf"$\phi_a$ ({label}{t_str})")

    run_label = os.path.basename(run_dir.rstrip("/")) or run_dir
    ax.set_title(f"Phase fields — first vs. last snapshot\n({run_label})", fontsize=13)
    ax.set_ylabel("Volume fraction", fontsize=12)
    ax.set_xlabel("x  [mm]", fontsize=12)
    y_lo_all = min(-0.05, ax.get_ylim()[0])
    ax.set_ylim(y_lo_all, 1.1)
    ax.axhline(0.0, color="gray", lw=0.8, ls=":")
    ax.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax.legend(fontsize=10, loc="best", ncol=2)
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)

    plt.tight_layout()

    if save_path is None:
        save_path = os.path.join(run_dir, "phase_first_last.png")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"First/last comparison saved to: {save_path}")


# ---------------------------------------------------------------------------
# Thermal overlay — all snapshots, twin y-axes, time colorbar
# ---------------------------------------------------------------------------

def plot_thermal_overlay(run_dir: str, save_path: str = None,
                          iga_file: str = "igasol.dat",
                          max_steps: int = None):
    """
    Single-panel twin-axis figure with all snapshots overlaid:
    T(x) on the left y-axis (solid lines), ρ_v(x) on the right y-axis
    (dashed lines), both color-coded by simulation time.
    """
    x_mm, _, _, tem_list, rhov_list, step_list, time_list = _collect_snapshots(
        run_dir, iga_file, max_steps
    )

    n = len(tem_list)

    t_vals = np.array([
        t if t is not None else float(s)
        for s, t in zip(step_list, time_list)
    ])
    use_time = any(t is not None for t in time_list)
    cb_label = "Time  [h]" if use_time else "Snapshot index"

    t_min, t_max = t_vals[0], t_vals[-1]
    if t_min == t_max:
        t_max = t_min + 1.0

    norm_t   = plt.Normalize(vmin=t_min, vmax=t_max)
    cmap_time = cm.plasma

    fig, ax_T = plt.subplots(1, 1, figsize=(10, 6))
    ax_rho = ax_T.twinx()

    for i, (tem, rhov, tv) in enumerate(zip(tem_list, rhov_list, t_vals)):
        c  = cmap_time(norm_t(tv))
        lw = 2.0 if i in (0, n - 1) else 1.2
        ax_T.plot(x_mm, tem,  color=c, lw=lw, ls="-")
        ax_rho.plot(x_mm, rhov, color=c, lw=lw, ls="--")

    # Reference saturation vapor density at the final mean temperature
    T_mean   = float(np.mean(tem_list[-1]))
    rvs_mean = rho_vs(np.array([T_mean]))[0]
    ax_rho.axhline(rvs_mean, color="k", ls=":", lw=1.2,
                   label=rf"$\rho_{{vs}}(T={T_mean:.1f}°C) = {rvs_mean:.3e}$ kg/m³")
    ax_rho.legend(fontsize=9, loc="upper right")

    # Proxy artists for line-style legend
    from matplotlib.lines import Line2D
    proxy_T   = Line2D([0], [0], color="gray", lw=2.0, ls="-",  label="Temperature  [°C]")
    proxy_rho = Line2D([0], [0], color="gray", lw=2.0, ls="--", label=r"Vapor density  [kg m$^{-3}$]")
    ax_T.legend(handles=[proxy_T, proxy_rho], fontsize=11, loc="upper left")

    ax_T.set_xlabel("x  [mm]", fontsize=12)
    ax_T.set_ylabel("Temperature  [°C]", fontsize=12, color=_COLOR_T)
    ax_rho.set_ylabel(r"Vapor density  [kg m$^{-3}$]", fontsize=12, color=_COLOR_RHOV)
    ax_T.tick_params(axis="y", labelcolor=_COLOR_T, labelsize=10)
    ax_rho.tick_params(axis="y", labelcolor=_COLOR_RHOV, labelsize=10)
    ax_T.tick_params(axis="x", labelsize=10)

    run_label = os.path.basename(run_dir.rstrip("/")) or run_dir
    ax_T.set_title(f"1D Permafrost — thermal fields (all snapshots)\n({run_label})",
                   fontsize=13)
    ax_T.grid(True, alpha=0.3)

    # Shared time colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap_time, norm=norm_t)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=[ax_T], orientation="vertical",
                        fraction=0.025, pad=0.10)
    cbar.set_label(cb_label, fontsize=11)

    plt.tight_layout()

    if save_path is None:
        thermal_dir = os.path.join(run_dir, "thermal")
        os.makedirs(thermal_dir, exist_ok=True)
        save_path = os.path.join(thermal_dir, "thermal_overlay.png")

    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Thermal overlay saved to: {save_path}")


# ---------------------------------------------------------------------------
# Derived quantities (time-series)
# ---------------------------------------------------------------------------

def plot_derived(run_dir: str, save_path: str = None,
                 iga_file: str = "igasol.dat"):
    """Time-series of ice volume fraction, interface positions, slab width."""
    nrb = load_geometry(os.path.join(run_dir, iga_file))
    x   = get_x_coords(nrb)
    Lx  = x[-1] - x[0]

    sol_files = sorted(
        glob.glob(os.path.join(run_dir, "sol_*.dat")), key=step_number
    )
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")

    times = load_times(run_dir)

    steps       = []
    t_h_arr     = []
    ice_volumes = []
    left_edges  = []
    right_edges = []
    slab_widths = []

    for sf in sol_files:
        try:
            sol = PetIGA().read_vec(sf, nrb)
        except Exception:
            continue
        if sol.ndim < 2:
            continue

        step = step_number(sf)
        ice  = sol[:, 0]
        vol  = np.trapz(ice, x) / Lx
        crossings = _find_crossings(x, ice, level=0.5)

        steps.append(step)
        t_h_arr.append(times.get(step, None))
        ice_volumes.append(vol)

        if len(crossings) >= 2:
            left_edges.append(crossings[0]  * 1e3)
            right_edges.append(crossings[-1] * 1e3)
            slab_widths.append((crossings[-1] - crossings[0]) * 1e3)
        elif len(crossings) == 1:
            left_edges.append(crossings[0] * 1e3)
            right_edges.append(np.nan)
            slab_widths.append(np.nan)
        else:
            left_edges.append(np.nan)
            right_edges.append(np.nan)
            slab_widths.append(np.nan)

    if not steps:
        sys.exit("No valid solutions found.")

    use_time = any(t is not None for t in t_h_arr)
    if use_time:
        xdata  = np.array([t if t is not None else np.nan for t in t_h_arr])
        xlabel = "Time  [h]"
    else:
        xdata  = np.array(steps, dtype=float)
        xlabel = "Step index"

    ice_volumes = np.array(ice_volumes)
    left_edges  = np.array(left_edges)
    right_edges = np.array(right_edges)
    slab_widths = np.array(slab_widths)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    axes[0].plot(xdata, ice_volumes, "o-", color="#1f77b4", ms=4)
    axes[0].set_ylabel(r"Mean ice fraction  $\langle\phi_i\rangle$", fontsize=12)
    axes[0].axhline(ice_volumes[0], color="gray", ls="--", lw=1, label="initial")
    axes[0].legend(fontsize=10)

    axes[1].plot(xdata, left_edges,  "s-", color="#2ca02c", ms=4, label="left edge")
    axes[1].plot(xdata, right_edges, "D-", color="#d62728", ms=4, label="right edge")
    axes[1].set_ylabel("Interface position  [mm]", fontsize=12)
    axes[1].legend(fontsize=10)

    axes[2].plot(xdata, slab_widths, "^-", color="#9467bd", ms=4)
    axes[2].set_ylabel("Slab width  [mm]", fontsize=12)
    axes[2].set_xlabel(xlabel, fontsize=12)

    run_label = os.path.basename(run_dir.rstrip("/")) or run_dir
    axes[0].set_title(f"1D Permafrost — derived quantities\n({run_label})", fontsize=13)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=10)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Derived figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot 1D permafrost field profiles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--dir",        default=".",
                   help="Run directory with sol_*.dat files (default: .)")
    p.add_argument("--iga",        default="igasol.dat",
                   help="IGA geometry file (default: igasol.dat)")
    p.add_argument("--max-steps",  type=int, default=None,
                   help="Maximum number of snapshots to load")
    p.add_argument("--out-dir",    default=None,
                   help="Directory for per-step PNGs (default: same as --dir)")
    p.add_argument("--thermal",    action="store_true",
                   help="Also produce the thermal overlay (all steps on one twin-axis plot)")
    p.add_argument("--derived",    action="store_true",
                   help="Plot derived scalar quantities (ice volume, interface position, width)")
    p.add_argument("--first-last", action="store_true",
                   help="Only produce the first/last comparison figure (skip per-step PNGs)")
    p.add_argument("--save",       default=None,
                   help="Save path for --thermal or --derived figure")
    return p.parse_args()


def main():
    args    = parse_args()
    run_dir = args.dir
    out_dir = args.out_dir or run_dir

    if args.derived:
        plot_derived(run_dir, save_path=args.save, iga_file=args.iga)
        return

    if args.first_last:
        plot_first_last(run_dir, save_path=args.save, iga_file=args.iga)
        return

    # Default: per-step phase PNGs + per-step thermal PNGs + first/last comparison
    plot_phase_steps(run_dir, out_dir=out_dir,
                     iga_file=args.iga, max_steps=args.max_steps)

    plot_thermal_steps(run_dir, out_dir=out_dir,
                       iga_file=args.iga, max_steps=args.max_steps)

    plot_first_last(run_dir, iga_file=args.iga)

    # Optional: thermal overlay (all snapshots on one twin-axis figure)
    if args.thermal:
        plot_thermal_overlay(run_dir, save_path=args.save,
                             iga_file=args.iga, max_steps=args.max_steps)


if __name__ == "__main__":
    main()
