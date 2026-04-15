# plot_summary_all.py
import os, argparse, numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dsm_plotlib import (set_default_rc, collect_runs, color_for, beautify_axes)

def main():
    ap = argparse.ArgumentParser(description="Two-panel summary across all runs")
    ap.add_argument("parent_dir")
    ap.add_argument("--pattern", default=None)
    ap.add_argument("--no-normalize", action="store_true")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--outfile", default="kxx_kyy_vs_ssa_and_time_all.png")
    args = ap.parse_args()
    normalize = not args.no_normalize

    set_default_rc()
    runs = collect_runs(args.parent_dir, args.pattern or None, normalize=normalize)
    if not runs:
        print("No runs found."); return

    fig, (axL, axR) = plt.subplots(1,2, figsize=(14.5,4.6), sharey=False, gridspec_kw={'wspace':0.18})

    # Choose most common time unit across runs (for axis label)
    units = [d.time_unit for d in runs]
    unit = max(set(units), key=units.count) if units else 'seconds'

    # Left: vs time
    for i, d in enumerate(runs):
        col = color_for(i)
        t = d.time_scaled
        vx = np.isfinite(t) & np.isfinite(d.k_xx)
        if np.count_nonzero(vx): axL.plot(t[vx], d.k_xx[vx], lw=2.0, marker='o', ms=4, markevery=2, color=col)
        vy = np.isfinite(t) & np.isfinite(d.k_yy)
        if np.count_nonzero(vy): axL.plot(t[vy], d.k_yy[vy], lw=2.0, linestyle='--', marker='s', ms=4, markevery=2, color=col)

    axL.set_xlabel(f'Time ({unit})')
    axL.set_ylabel(r'$k/k_0$' if normalize else r'$k$ (W/m·K)')
    axL.set_title(r'$k_{xx}$ & $k_{yy}$ vs Time' if normalize else 'k vs Time')
    beautify_axes(axL)

    # Right: vs SSA
    for i, d in enumerate(runs):
        if d.ssa is None or d.ssa.size==0: continue
        col = color_for(i)
        vx = np.isfinite(d.ssa) & np.isfinite(d.k_xx)
        if np.count_nonzero(vx): axR.plot(d.ssa[vx], d.k_xx[vx], lw=2.0, marker='o', ms=4, markevery=2, color=col)
        vy = np.isfinite(d.ssa) & np.isfinite(d.k_yy)
        if np.count_nonzero(vy): axR.plot(d.ssa[vy], d.k_yy[vy], lw=2.0, linestyle='--', marker='s', ms=4, markevery=2, color=col)

    axR.set_xlabel(r'SSA/SSA$_0$' if normalize else 'SSA')
    axR.set_ylabel(r'$k/k_0$' if normalize else r'$k$ (W/m·K)')
    axR.set_title(r'$k_{xx}$ & $k_{yy}$ vs SSA' if normalize else 'k vs SSA')
    beautify_axes(axR)

    component_handles = [
        Line2D([0],[0], color='black', lw=2.0, linestyle='-', marker='o', ms=6, label=r'$k_{xx}$'),
        Line2D([0],[0], color='black', lw=2.0, linestyle='--', marker='s', ms=6, label=r'$k_{yy}$')
    ]
    fig.legend(component_handles, [h.get_label() for h in component_handles], loc='center left', bbox_to_anchor=(0.995, 0.5))
    fig.subplots_adjust(right=0.92)
    fig.tight_layout(rect=[0.0,0.0,0.92,1.0], pad=0.2)
    fig.savefig(os.path.join(args.parent_dir, args.outfile), dpi=args.dpi, bbox_inches='tight', transparent=True)
    plt.close(fig)

if __name__ == "__main__":
    main()