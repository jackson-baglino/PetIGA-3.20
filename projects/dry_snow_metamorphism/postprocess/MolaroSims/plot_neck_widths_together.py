"""
plot_neck_widths_together.py

Plot neck-width evolution for two simulations on the same axes.

Assumptions
-----------
- Each simulation has already been processed by `plot_neck_width.py`
  (or a similar script) to produce a CSV file in the *parent* directory
  of the `vtkOut` folder.
- That CSV file is named `neck_width_vs_time.csv` and contains at least:
    - time in seconds or minutes   (columns like: time_s, time_min, or time)
    - neck width in meters or µm   (columns like: neck_width_m, neck_width_um)
- If your column names differ, you can adjust the small helper function
  `_extract_time_and_neck` below.

How to use
----------
1. Set `PARENT_DIR_1` and `PARENT_DIR_2` to the parent directories of
   the simulations you want to compare (the same directories that contain
   `SSA_evo.dat` and `vtkOut/`).
2. Optionally adjust `CSV_NAME` if your output file has a different name.
3. Run:
       python plot_neck_widths_together.py
"""

from pathlib import Path
from typing import Tuple, List, Dict, Optional

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import re

# ============================
# USER PARAMETERS
# ============================

# Fill these in with the parent directories of your two simulations
PARENT_DIR_1 = Path(
    "/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/archived_results/BestParams/Group1/NASAv2-Molaro2D_T-05.0_hum_2024-11-26__09.41.44"
)
PARENT_DIR_2 = Path(
    "/Users/jacksonbaglino/SimulationResults/dry_snow_metamorphism/archived_results/BestParams/Group1/NASAv2-Molaro2D_T-20.0_hum_2024-11-26__09.41.38"
)

# Name of the CSV file produced by the single-simulation script
CSV_NAME = "neck_widths_vs_time.csv"

# Optional: custom labels for the legend.
# If left as None, the script will use the last directory name of each parent.
LABEL_1 = None
LABEL_2 = None

# Path to your Matplotlib style file (same one used in plot_neck_width.py)
STYLE_FILE = Path(__file__).with_name("plot_style.mplstyle")


# Output figure name
OUTPUT_FIG = "neck_widths_comparison.png"

#
# Output formats are organized into subfolders by format.
OUTPUT_FORMATS = ["pdf", "svg"]

# Output directory structure:
#   <SCRIPT_DIR>/plots/pdf/
#   <SCRIPT_DIR>/plots/svg/
PLOTS_DIR = Path(__file__).with_name("plots")
PDF_DIR = PLOTS_DIR / "pdf"
SVG_DIR = PLOTS_DIR / "svg"

# Save two variants of each plot:
#   - with title
#   - without title
SAVE_WITH_TITLE = True
SAVE_NO_TITLE = True


# Preferred font family for all plot text. Falls back to a default sans-serif if unavailable.
FONT_FAMILY = "Montserrat"

# Optional: directory containing local font files (*.ttf, *.otf).
# If you downloaded Montserrat from Google Fonts, copy the extracted folder
# (or just the .ttf files) into a folder named `fonts/` next to this script.
# Example:
#   postprocess/MolaroSims/fonts/Montserrat-Regular.ttf
FONT_DIR = Path(__file__).with_name("fonts")

# If True, attempt to register fonts from FONT_DIR before applying rcParams.
REGISTER_LOCAL_FONTS = True

PLOT_EXPERIMENTAL = True

# MATLAB-like default colors (R2014b+)
# [0, 0.4470, 0.7410] -> "#0072BD" (blue)
# [0.8500, 0.3250, 0.0980] -> "#D95319" (orange)
COLOR_5C = "#0072BD"   # MATLAB blue for ±5 °C
COLOR_20C = "#D95319"  # MATLAB orange for ±20 °C

# ------------------------------------------------------------------
# Hard-coded experimental data (absolute neck width and errors)
# ------------------------------------------------------------------
# All times are in minutes and neck widths / errors are in microns.
# Replace the example arrays below with the values from your table.
#
# -5 °C experimental data
USE_EXP_5 = True
EXP_TIME_5_MIN = np.array([
    0.0,
    26.0,
    48.0,
])
EXP_NECK_5_UM = np.array([
    32.51,
    41.57,
    44.00,
])
# Asymmetric errors for -5 °C (upper and lower, in µm)
EXP_ERR_5_UP_UM = np.array([
    3.7,
    0.7,
    1.1,
])
EXP_ERR_5_DN_UM = np.array([
    2.8,
    0.5,
    1.0,
])

# -20 °C experimental data
USE_EXP_20 = True
EXP_TIME_20_MIN = np.array([
    0.0,
    2.0,
    5.0,
    9.0,
    18.0,
    29.0,
    49.0,
    57.0,
    78.0,
])
EXP_NECK_20_UM = np.array([
    32.81,
    33.60,
    37.16,
    44.68,
    46.84,
    54.18,
    58.00,
    60.54,
    64.78,
])
# Asymmetric errors for -20 °C (upper and lower, in µm)
EXP_ERR_20_UP_UM = np.array([
    2.3,
    3.9,
    0.7,
    1.3,
    3.1,
    1.6,
    1.1,
    1.2,
    2.0,
])
EXP_ERR_20_DN_UM = np.array([
    2.4,
    2.4,
    0.6,
    1.7,
    2.5,
    3.2,
    2.1,
    1.1,
    1.7
])

# Optional time shifts (in minutes) for experimental curves.
# Positive values move the experimental curves to the left (toward t=0),
# after which times are clipped at a minimum of 0.
EXP_TIME_SHIFT_5 = 0.0   # shift for the -5 °C experimental data
EXP_TIME_SHIFT_20 = 0.0  # shift for the -20 °C experimental data




# ============================
# HELPER FUNCTIONS
# ============================

def _ensure_output_dirs() -> None:
    """Create output directories for vector formats."""
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    PDF_DIR.mkdir(parents=True, exist_ok=True)
    SVG_DIR.mkdir(parents=True, exist_ok=True)


def _save_fig_variants(fig: plt.Figure, ax: plt.Axes, stem: str) -> None:
    """Save the figure in per-format subfolders, with and without a title."""
    _ensure_output_dirs()

    title_text = ax.get_title()

    def _format_dir(fmt: str) -> Path:
        fmt = fmt.lower()
        if fmt == "pdf":
            return PDF_DIR
        if fmt == "svg":
            return SVG_DIR
        # fallback: save under plots/<fmt>
        d = PLOTS_DIR / fmt
        d.mkdir(parents=True, exist_ok=True)
        return d

    # Save WITH title
    if SAVE_WITH_TITLE:
        if title_text:
            ax.set_title(title_text)
        for fmt in OUTPUT_FORMATS:
            out_path = _format_dir(fmt) / f"{stem}.{fmt}"
            fig.savefig(out_path, bbox_inches="tight", transparent=True)

    # Save WITHOUT title
    if SAVE_NO_TITLE:
        ax.set_title("")
        for fmt in OUTPUT_FORMATS:
            out_path = _format_dir(fmt) / f"{stem}_noTitle.{fmt}"
            fig.savefig(out_path, bbox_inches="tight", transparent=True)
        # Restore
        ax.set_title(title_text)


def _register_local_fonts(font_dir: Optional[Path]) -> None:
    """Register local .ttf/.otf fonts with Matplotlib's font manager."""
    if not font_dir:
        return
    try:
        font_dir = Path(font_dir)
    except Exception:
        return

    if not font_dir.exists():
        return

    # Register any TTF/OTF fonts found in this directory (recursively).
    font_files = list(font_dir.rglob("*.ttf")) + list(font_dir.rglob("*.otf"))
    if not font_files:
        return

    added = 0
    for fp in font_files:
        try:
            font_manager.fontManager.addfont(str(fp))
            added += 1
        except Exception:
            # Skip fonts that Matplotlib cannot load.
            continue

    if added > 0:
        # Trigger font list refresh (best-effort; different Matplotlib versions behave differently).
        try:
            font_manager._get_font.cache_clear()  # type: ignore[attr-defined]
        except Exception:
            pass


def _pick_font_family(preferred: str) -> str:
    """Return preferred font family if installed (or locally registered), else fall back safely."""
    # Try to register local fonts first (so the preferred family can become available).
    if REGISTER_LOCAL_FONTS:
        _register_local_fonts(FONT_DIR)

    try:
        available = {f.name for f in font_manager.fontManager.ttflist}
        if preferred in available:
            return preferred
    except Exception:
        # If font discovery fails for any reason, fall back.
        pass
    return "DejaVu Sans"


def _apply_font_settings(preferred: str) -> None:
    """Apply global Matplotlib font settings."""
    if REGISTER_LOCAL_FONTS:
        _register_local_fonts(FONT_DIR)
    chosen = _pick_font_family(preferred)
    if chosen != preferred:
        print(f"[WARN] Font '{preferred}' not found. Falling back to {chosen}.")

    # Use a list so Matplotlib can fall back to other sans-serif fonts if needed.
    plt.rcParams["font.family"] = [chosen]
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans", "Arial", "Helvetica"]

    # Keep vector outputs text as text (not paths) where possible.
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"

def _parse_label_radius_temp(label: str) -> Tuple[Optional[float], Optional[float]]:
    """
    Parse a label like '99µm at -5ºC' into (radius_um, temperature_C).

    Returns (None, None) if parsing fails.
    """
    if not label:
        return None, None

    # Radius: number before 'µm'
    radius_um: Optional[float] = None
    m_r = re.search(r"([-+]?\d+\.?\d*)\s*µm", label)
    if m_r:
        try:
            radius_um = float(m_r.group(1))
        except ValueError:
            radius_um = None

    # Temperature: number before 'C', possibly with degree symbol
    temp_C: Optional[float] = None
    m_t = re.search(r"([-+]?\d+\.?\d*)\s*[°º]?\s*C", label)
    if m_t:
        try:
            temp_C = float(m_t.group(1))
        except ValueError:
            temp_C = None

    return radius_um, temp_C


def _load_experimental_data(csv_path: Path) -> List[Dict]:
    """
    Load experimental neck data from a Molaro-style CSV with multiple
    (X, Y) column pairs and a first row of labels.

    Returns a list of dicts with keys:
        - label
        - radius_um
        - temperature_C
        - time_min  (np.ndarray)
        - rel_neck  (np.ndarray)
    """
    if not csv_path.is_file():
        print(f"[WARN] Experimental CSV not found: {csv_path}")
        return []

    # First row: labels for each (X,Y) pair
    labels_row = pd.read_csv(csv_path, header=None, nrows=1).iloc[0]

    # Remaining rows: data, with second row as header (X, Y, X, Y, ...)
    df = pd.read_csv(csv_path, header=1)

    n_cols = df.shape[1]
    if n_cols % 2 != 0:
        raise ValueError(
            f"Expected an even number of columns in experimental CSV, got {n_cols}."
        )

    datasets: List[Dict] = []

    for i in range(0, n_cols, 2):
        label = str(labels_row[i]).strip()
        if not label or label.lower() == "nan":
            continue

        xcol = df.columns[i]
        ycol = df.columns[i + 1]

        x = df[xcol]
        y = df[ycol]

        # Drop rows where either x or y is NaN
        mask = ~(x.isna() | y.isna())
        x_vals = x[mask].to_numpy()
        y_vals = y[mask].to_numpy()

        if x_vals.size == 0 or y_vals.size == 0:
            continue

        radius_um, temp_C = _parse_label_radius_temp(label)

        datasets.append(
            dict(
                label=label,
                radius_um=radius_um,
                temperature_C=temp_C,
                rel_neck=x_vals,
                time_min=y_vals,  # experimental time is already in minutes
            )
        )

    return datasets



def _average_experimental_datasets(
    datasets: List[Dict],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Given multiple experimental datasets for a single temperature, build
    an averaged absolute neck-width curve vs time.

    Each dataset dict must contain:
        - radius_um
        - time_min  (np.ndarray)
        - rel_neck  (np.ndarray)

    Behavior:
      - If only one dataset is present, convert its relative neck size to
        absolute neck width (µm) using:
            neck_width_um = RELATIVE_TO_ABS_FACTOR * rel_neck * (2.0 * radius_um)
        and return its (time_min, neck_width_um) unchanged.
      - If two or more datasets are present:
        - Select the first dataset as the reference and use its time_min
          as the reference time grid.
        - For each dataset, compute its absolute neck width curve:
            width_um = RELATIVE_TO_ABS_FACTOR * rel_neck * (2.0 * radius_um)
        - For the reference dataset, keep its width curve aligned with t_ref.
        - For each additional dataset, linearly interpolate its width curve
          onto t_ref using np.interp; mark values outside that dataset's
          original time range as np.nan so they do not contribute to the average.
        - Stack all width curves and compute the mean at each t_ref using
          np.nanmean.
        - Drop any entries where the averaged width is NaN.

    Returns:
        (time_min, avg_neck_width_um)
    """
    if not datasets:
        return np.array([]), np.array([])

    if len(datasets) == 1:
        g = datasets[0]
        t = g["time_min"]
        rel = g["rel_neck"]
        radius_um = g["radius_um"]
        width_um = RELATIVE_TO_ABS_FACTOR * rel * (2.0 * radius_um)
        return t, width_um

    # Use the first dataset as the reference time grid
    ref = datasets[0]
    t_ref = ref["time_min"]
    rel_ref = ref["rel_neck"]
    radius_ref = ref["radius_um"]
    width_ref = RELATIVE_TO_ABS_FACTOR * rel_ref * (2.0 * radius_ref)

    width_curves = [width_ref]

    for g in datasets[1:]:
        t = g["time_min"]
        rel = g["rel_neck"]
        radius_um = g["radius_um"]
        width_um = RELATIVE_TO_ABS_FACTOR * rel * (2.0 * radius_um)

        # Interpolate onto t_ref, mark outside original time range as np.nan
        w_interp = np.interp(t_ref, t, width_um)
        valid = (t_ref >= t.min()) & (t_ref <= t.max())
        w_interp[~valid] = np.nan
        width_curves.append(w_interp)

    arr = np.vstack(width_curves)
    avg_width = np.nanmean(arr, axis=0)
    valid = ~np.isnan(avg_width)
    t_valid = t_ref[valid]
    width_valid = avg_width[valid]

    return t_valid, width_valid


def _extract_time_and_neck(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    """
    Given a DataFrame loaded from the neck-width CSV, extract
    time (in minutes) and neck width (in microns).

    This is written to be robust to a few different possible
    column-name conventions. Adjust as needed if your CSV
    uses different names.
    """
    # ----- Time -----
    if "time_min" in df.columns:
        time_min = df["time_min"]
    elif "time_s" in df.columns:
        time_min = df["time_s"] / 60.0
    elif "time" in df.columns:
        # Assume already in seconds
        time_min = df["time"] / 60.0
    else:
        raise KeyError(
            "Could not find a time column. Expected one of: "
            "'time_min', 'time_s', or 'time'."
        )

    # ----- Neck width -----
    if "neck_width_um" in df.columns:
        neck_um = df["neck_width_um"]
    elif "neck_width_m" in df.columns:
        neck_um = df["neck_width_m"] * 1e6
    elif "neck_length" in df.columns:
        # Fallback to generic name, assume meters
        neck_um = df["neck_length"] * 1e6
    else:
        raise KeyError(
            "Could not find a neck-width column. Expected one of: "
            "'neck_width_um', 'neck_width_m', or 'neck_length'."
        )

    return time_min, neck_um


def _load_neck_data(parent_dir: Path) -> Tuple[pd.Series, pd.Series]:
    """
    Load neck-width CSV for a given parent directory and return
    (time_min, neck_width_um).

    Primary behavior:
        - Look for a file named exactly `CSV_NAME` in `parent_dir`.

    Fallback behavior (if that file is not found):
        - Search for any CSV in `parent_dir` whose name starts with "neck"
          and ends with ".csv".
        - If exactly one such file is found, use that and print a diagnostic
          message.
        - If multiple candidates are found, raise an error listing them so
          the user can choose/rename the correct one.
    """
    csv_path = parent_dir / CSV_NAME

    if not csv_path.is_file():
        # Fallback: search for plausible candidates
        candidates = sorted(parent_dir.glob("neck*.csv"))
        if len(candidates) == 1:
            csv_path = candidates[0]
            print(f"[INFO] Using detected CSV file for {parent_dir.name}: {csv_path.name}")
        elif len(candidates) > 1:
            msg_lines = [
                f"Could not find CSV file named '{CSV_NAME}' in: {parent_dir}",
                "However, multiple candidate neck-width CSV files were found:",
            ]
            for c in candidates:
                msg_lines.append(f"  - {c.name}")
            msg_lines.append(
                "Please either set CSV_NAME to the desired file name, "
                "or rename the correct CSV to match CSV_NAME."
            )
            raise FileNotFoundError("\n".join(msg_lines))
        else:
            raise FileNotFoundError(
                f"Could not find CSV file: {csv_path}\n"
                f"No files matching 'neck*.csv' were found in {parent_dir}."
            )

    df = pd.read_csv(csv_path)
    return _extract_time_and_neck(df)


def _plot_simulation_curves(
    ax: plt.Axes,
    t1_min: pd.Series,
    w1_um: pd.Series,
    label1: str,
    t2_min: pd.Series,
    w2_um: pd.Series,
    label2: str,
    color1: Optional[str] = None,
    color2: Optional[str] = None,
) -> None:
    ax.plot(t1_min, w1_um, linewidth=2.0, label="_nolegend_", color=color1)
    ax.plot(t2_min, w2_um, linewidth=2.0, label="_nolegend_", color=color2)


# ============================
# MAIN
# ============================

def main() -> None:
    # Apply common plot style if available
    if STYLE_FILE.is_file():
        plt.style.use(STYLE_FILE)

    # Apply font settings after style so they take precedence
    _apply_font_settings(FONT_FAMILY)
    _ensure_output_dirs()

    # Load data for both simulations
    t1_min, w1_um = _load_neck_data(PARENT_DIR_1)
    t2_min, w2_um = _load_neck_data(PARENT_DIR_2)

    # Determine labels
    def extract_temp_from_folder(name):
        # Look for T-05.0, T+10.0, T-20.0, etc.
        m = re.search(r"T([+-]?\d+\.?\d*)", name)
        if m:
            try:
                return float(m.group(1))
            except Exception:
                return None
        return None

    label1_base = LABEL_1 or PARENT_DIR_1.name
    label2_base = LABEL_2 or PARENT_DIR_2.name
    temp1 = extract_temp_from_folder(PARENT_DIR_1.name)
    temp2 = extract_temp_from_folder(PARENT_DIR_2.name)
    # Keep the sign (negative temperatures should display as negative)
    label1 = f"Simulation, T={temp1:.0f}°C" if temp1 is not None else label1_base
    label2 = f"Simulation, T={temp2:.0f}°C" if temp2 is not None else label2_base

    # Assign colors for each simulation based on temperature
    def _color_for_temp(temp: Optional[float]) -> Optional[str]:
        if temp is None:
            return None
        t_abs = abs(temp)
        if np.isclose(t_abs, 5.0):
            return COLOR_5C
        if np.isclose(t_abs, 20.0):
            return COLOR_20C
        return None

    color1 = _color_for_temp(temp1)
    color2 = _color_for_temp(temp2)

    # Determine simulation time span
    t_sim_max = float(max(t1_min.max(), t2_min.max()))

    # Load experimental data if requested
    # exp_curves holds tuples: (time_min, neck_width_um, err_um_or_None, label)
    exp_curves: List[Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], str]] = []
    t_exp_max: Optional[float] = None

    if PLOT_EXPERIMENTAL:
        # -5 °C experimental data
        if USE_EXP_5 and EXP_TIME_5_MIN.size > 0:
            t5 = EXP_TIME_5_MIN
            w5 = EXP_NECK_5_UM
            if (
                "EXP_ERR_5_UP_UM" in globals()
                and "EXP_ERR_5_DN_UM" in globals()
                and EXP_ERR_5_UP_UM.size == EXP_TIME_5_MIN.size
                and EXP_ERR_5_DN_UM.size == EXP_TIME_5_MIN.size
            ):
                # yerr expects shape (2, N): [lower, upper]
                err5 = np.vstack([EXP_ERR_5_DN_UM, EXP_ERR_5_UP_UM])
            else:
                err5 = None

            # Apply optional time shift and clip at t=0
            t5_shift = t5 - EXP_TIME_SHIFT_5
            t5_shift = np.clip(t5_shift, a_min=0.0, a_max=None)

            if t_exp_max is None:
                t_exp_max = float(t5_shift.max())
            else:
                t_exp_max = max(t_exp_max, float(t5_shift.max()))

            exp_curves.append((t5_shift, w5, err5, "Experimental, T=-5°C"))

        # -20 °C experimental data
        if USE_EXP_20 and EXP_TIME_20_MIN.size > 0:
            t20 = EXP_TIME_20_MIN
            w20 = EXP_NECK_20_UM
            if (
                "EXP_ERR_20_UP_UM" in globals()
                and "EXP_ERR_20_DN_UM" in globals()
                and EXP_ERR_20_UP_UM.size == EXP_TIME_20_MIN.size
                and EXP_ERR_20_DN_UM.size == EXP_TIME_20_MIN.size
            ):
                # yerr expects shape (2, N): [lower, upper]
                err20 = np.vstack([EXP_ERR_20_DN_UM, EXP_ERR_20_UP_UM])
            else:
                err20 = None

            # Apply optional time shift and clip at t=0
            t20_shift = t20 - EXP_TIME_SHIFT_20
            t20_shift = np.clip(t20_shift, a_min=0.0, a_max=None)

            if t_exp_max is None:
                t_exp_max = float(t20_shift.max())
            else:
                t_exp_max = max(t_exp_max, float(t20_shift.max()))

            exp_curves.append((t20_shift, w20, err20, "Experimental, T=-20°C"))

    # ---------- Figure 1: limited to experimental duration ----------
    if exp_curves and t_exp_max is not None:
        fig1, ax1 = plt.subplots(figsize=(4.5, 3.375))

        _plot_simulation_curves(
            ax1,
            t1_min,
            w1_um,
            label1,
            t2_min,
            w2_um,
            label2,
            color1=color1,
            color2=color2,
        )

        for t_exp, w_exp_um, err_exp_um, lbl in exp_curves:
            # Match experimental color to the corresponding simulation temperature
            if "-5°C" in lbl:
                c = COLOR_5C
            elif "-20°C" in lbl:
                c = COLOR_20C
            else:
                c = None

            if err_exp_um is not None:
                ax1.errorbar(
                    t_exp,
                    w_exp_um,
                    yerr=err_exp_um,
                    linestyle="none",
                    marker="o",
                    markersize=6,
                    capsize=3,
                    color=c,
                    markerfacecolor="none",
                    markeredgewidth=1.2,
                    label="_nolegend_",
                )
            else:
                ax1.plot(
                    t_exp,
                    w_exp_um,
                    linestyle="none",
                    marker="o",
                    markersize=6,
                    color=c,
                    markerfacecolor="none",
                    markeredgewidth=1.2,
                    label="_nolegend_",
                )

        left1 = -0.02 * t_exp_max if t_exp_max > 0 else -0.1
        ax1.set_xlim(left=left1, right=t_exp_max * 1.05)
        ax1.set_xlabel("Time (min)")
        ax1.set_ylabel("Neck width (µm)")
        ax1.set_title("Neck width evolution")

        from matplotlib.lines import Line2D

        legend_elements = [
            Line2D([0], [0], color=COLOR_5C, lw=2, label="T=-5°C"),
            Line2D([0], [0], color=COLOR_20C, lw=2, label="T=-20°C"),
        ]

        ax1.legend(
            handles=legend_elements,
            frameon=True,
            loc="best",
        )
        # ax1.grid(False, which="both", linestyle="--", alpha=0.3)

        fig1.tight_layout()
        _save_fig_variants(fig1, ax1, stem="neck_widths_comparison_expWindow")
        print(f"Saved experimental-window comparison figures to: {PLOTS_DIR}")
        plt.close(fig1)

    # ---------- Figure 2: full simulation duration ----------
    fig2, ax2 = plt.subplots(figsize=(4.5, 3.375))

    _plot_simulation_curves(
        ax2,
        t1_min,
        w1_um,
        label1,
        t2_min,
        w2_um,
        label2,
        color1=color1,
        color2=color2,
    )

    for t_exp, w_exp_um, err_exp_um, lbl in exp_curves:
        # Match experimental color to the corresponding simulation temperature
        if "-5°C" in lbl:
            c = COLOR_5C
        elif "-20°C" in lbl:
            c = COLOR_20C
        else:
            c = None

        if err_exp_um is not None:
            ax2.errorbar(
                t_exp,
                w_exp_um,
                yerr=err_exp_um,
                linestyle="none",
                marker="o",
                markersize=6,
                capsize=3,
                color=c,
                markerfacecolor="none",
                markeredgewidth=1.2,
                label="_nolegend_",
            )
        else:
            ax2.plot(
                t_exp,
                w_exp_um,
                linestyle="none",
                marker="o",
                markersize=6,
                color=c,
                markerfacecolor="none",
                markeredgewidth=1.2,
                label="_nolegend_",
            )

    left2 = -0.02 * t_sim_max if t_sim_max > 0 else -0.1
    ax2.set_xlim(left=left2, right=t_sim_max * 1.05)
    ax2.set_xlabel("Time (min)")
    ax2.set_ylabel("Neck width (µm)")
    ax2.set_title("Neck width evolution")

    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], color=COLOR_5C, lw=2, label="T=-5°C"),
        Line2D([0], [0], color=COLOR_20C, lw=2, label="T=-20°C"),
    ]

    ax2.legend(
        handles=legend_elements,
        frameon=True,
        loc="best",
    )
    # ax2.grid(False, which="both", linestyle="--", alpha=0.3)

    fig2.tight_layout()
    _save_fig_variants(fig2, ax2, stem="neck_widths_comparison_simWindow")
    print(f"Saved simulation-window comparison figures to: {PLOTS_DIR}")
    plt.close(fig2)


if __name__ == "__main__":
    main()
