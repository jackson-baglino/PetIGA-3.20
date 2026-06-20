#!/usr/bin/env python3
"""
plot_timestep.py — Log-log plot of adaptive time step size vs simulation time.

Parses the monitor table written to outp.txt (pipe-delimited rows) and plots
dt vs t on log-log axes.  Useful as a quick diagnostic for solver performance:
a steadily growing dt indicates a well-behaved adaptive stepper; sudden drops
signal rejected steps or hard events (e.g. the sediment freeze transition).

Usage
-----
  python plot_timestep.py --dir /path/to/run/folder
  python plot_timestep.py --dir /path/to/run/folder --save timestep.png
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Parser for outp.txt monitor table
# ---------------------------------------------------------------------------

# outp.txt contains two kinds of pipe-delimited rows that both start with a
# number followed by '|':
#
#   Domain-integral table (8 fields after the step number):
#     STEP | TIME [s] | DT [s] | TOT_ICE | TOT_AIR | TEMP | TOT_RHOV | I-A INTERF | TOTAL_MASS
#        0 |  0.00000e+00 | 1.000e-02 |  1.125e-09 |  ...
#
#   SNES/Newton iteration table (10 fields after the iteration number):
#     it | fnorm | n0 | r0 | s0 | n1 | r1 | s1 | n2 | r2 | s2
#      0 |  6.8134e-15 |   4.818e-15 | 1.000e+00 | ...
#
# We must match only domain-integral rows.  The reliable discriminator is
# field count: domain rows have exactly _DOMAIN_NFIELDS fields after the step.
_ROW_RE = re.compile(r"^\s*(\d+)\s*\|(.+)$")
_DOMAIN_NFIELDS = 8  # TIME, DT, TOT_ICE, TOT_AIR, TEMP, TOT_RHOV, I-A INTERF, TOTAL_MASS


def _parse_float(s: str) -> float:
    try:
        return float(s.strip())
    except ValueError:
        return float("nan")


def load_monitor(outp_path: str) -> np.ndarray:
    """
    Parse outp.txt and return an (N, 3) array of [step, t, dt].
    Rows where t == 0 are kept (step 0); NaN rows are dropped.
    """
    if not os.path.isfile(outp_path):
        sys.exit(f"ERROR: monitor file not found: {outp_path}")

    rows = []
    with open(outp_path) as fh:
        for line in fh:
            m = _ROW_RE.match(line)
            if not m:
                continue
            step = int(m.group(1))
            fields = m.group(2).split("|")
            if len(fields) != _DOMAIN_NFIELDS:   # skip SNES/Newton iteration rows
                continue
            t  = _parse_float(fields[0])   # TIME [s]
            dt = _parse_float(fields[1])   # DT [s]
            if not (np.isnan(t) or np.isnan(dt)):
                rows.append((step, t, dt))

    if not rows:
        sys.exit(f"ERROR: no monitor data rows found in {outp_path}")

    return np.array(rows, dtype=float)   # shape (N, 3): step, t, dt


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_timestep(run_dir: str, save_path: str = None, title: str = None):
    outp_path = os.path.join(run_dir, "outp.txt")
    data = load_monitor(outp_path)

    step = data[:, 0].astype(int)
    t    = data[:, 1]
    dt   = data[:, 2]

    # For log-log both axes must be positive; drop t == 0 (step 0 has no
    # meaningful "elapsed" time to plot on a log scale)
    mask = t > 0
    t_plot    = t[mask]
    dt_plot   = dt[mask]
    step_plot = step[mask]

    if len(t_plot) == 0:
        sys.exit("ERROR: all time values are zero — nothing to plot on log scale.")

    fig, (ax_step, ax_time) = plt.subplots(1, 2, figsize=(13, 4.5))

    # Step-indexed view: reveals reject/grow cycling and other per-step
    # structure that a time axis compresses away when dt spans many orders
    # of magnitude (e.g. the initial ramp-up from delt_t to dtmax).
    ax_step.semilogy(step_plot, dt_plot, color="#1f77b4", lw=1.5, marker=".", markersize=3)
    ax_step.set_xlabel("Step", fontsize=12)
    ax_step.set_ylabel("Time step  $\\Delta t$  [s]", fontsize=12)
    ax_step.set_title("vs. step", fontsize=12)
    ax_step.grid(True, which="both", alpha=0.3)
    ax_step.tick_params(labelsize=11)

    # Time-indexed view: shows dt growth relative to elapsed simulation
    # time -- a steadily growing dt indicates a well-behaved adaptive
    # stepper; sudden drops signal rejected steps or hard events.
    ax_time.loglog(t_plot, dt_plot, color="#1f77b4", lw=1.5, marker=".", markersize=3)
    ax_time.set_xlabel("Simulation time  $t$  [s]", fontsize=12)
    ax_time.set_ylabel("Time step  $\\Delta t$  [s]", fontsize=12)
    ax_time.set_title("vs. simulation time", fontsize=12)
    ax_time.grid(True, which="both", alpha=0.3)
    ax_time.tick_params(labelsize=11)

    fig.suptitle(title or "Adaptive time step history", fontsize=13)

    # Annotate total step count and final dt
    ax_time.text(0.98, 0.05,
            f"steps: {step[-1]}   final $\\Delta t$: {dt[-1]:.2e} s",
            transform=ax_time.transAxes, ha="right", va="bottom",
            fontsize=9, color="gray")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Time step figure saved to: {save_path}")
    else:
        plt.show()

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="Log-log plot of dt vs t from outp.txt monitor output."
    )
    p.add_argument("--dir",   default=".",
                   help="Run folder containing outp.txt (default: current dir)")
    p.add_argument("--save",  default=None,
                   help="Save figure to this path (omit to display interactively)")
    p.add_argument("--title", default=None,
                   help="Figure title override")
    args = p.parse_args()

    plot_timestep(run_dir=args.dir, save_path=args.save, title=args.title)


if __name__ == "__main__":
    main()
