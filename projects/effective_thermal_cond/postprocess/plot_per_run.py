# plot_per_run.py
import os, argparse
import numpy as np
import matplotlib.pyplot as plt
from dsm_plotlib import (set_default_rc, collect_runs, color_for, beautify_axes,
                         place_legend_outside, sparse_markevery, add_secondary_ssa_axis)

def main():
    ap = argparse.ArgumentParser(description="Per-folder DSM plots")
    ap.add_argument("parent_dir")
    ap.add_argument("--pattern", default=None, help="Regex for folder names")
    ap.add_argument("--no-normalize", action="store_true", help="Plot absolute values")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--dualx", action="store_true", help="Add top SSA axis to time plots")
    ap.add_argument("--show", action="store_true", help="Show instead of save")
    args = ap.parse_args()

    normalize = not args.no_normalize
    set_default_rc()
    runs = collect_runs(args.parent_dir, pattern=args.pattern or None, normalize=normalize)
    if not runs:
        print("No runs found.")
        return

    for rd in runs:
        folder_path = os.path.join(args.parent_dir, rd.folder)
        print(f"Processing {rd.folder} (T={rd.temperature} K)")

        # k_yy vs time
        fig, ax = plt.subplots(figsize=(11.28, 4.2))
        ax.plot(rd.time_scaled, rd.k_yy, marker='o', lw=1.5)
        ax.set_xlabel(f'Time ({rd.time_unit})')
        ax.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize else r'$k_{yy}$ (W/m·K)')
        ax.set_title(f'T = {rd.temperature:g} K')
        beautify_axes(ax)
        fig.tight_layout()
        if args.show: plt.show()
        else: fig.savefig(os.path.join(folder_path, 'kyy_vs_time.png'), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

        # k_xx vs time
        fig, ax = plt.subplots(figsize=(10.0, 3.0))
        ax.plot(rd.time_scaled, rd.k_xx, marker='o', lw=1.5)
        ax.set_xlabel(f'Time ({rd.time_unit})')
        ax.set_ylabel(r'$k_{xx}/k_{xx,0}$' if normalize else r'$k_{xx}$ (W/m·K)')
        ax.set_title(f'T = {rd.temperature:g} K')
        beautify_axes(ax)
        fig.tight_layout()
        if args.show: plt.show()
        else: fig.savefig(os.path.join(folder_path, 'kxx_vs_time.png'), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

        # dual y (kxx & kyy) vs time
        fig, ax1 = plt.subplots(figsize=(11.28, 4.0))
        col_xx, col_yy = 'C0', 'C1'
        ln1 = ax1.plot(rd.time_scaled, rd.k_xx, marker='o', lw=1.5, label=r'$k_{xx}$', color=col_xx)
        ax1.set_xlabel(f'Time ({rd.time_unit})')
        ax1.set_ylabel(r'$k_{xx}/k_{xx,0}$' if normalize else r'$k_{xx}$ (W/m·K)')
        ax1.set_title(f'T = {rd.temperature:g} K'); beautify_axes(ax1)
        ax2 = ax1.twinx()
        ln2 = ax2.plot(rd.time_scaled, rd.k_yy, marker='s', lw=1.5, label=r'$k_{yy}$', color=col_yy)
        ax2.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize else r'$k_{yy}$ (W/m·K)')
        rect = place_legend_outside(ax1, ln1+ln2, [l.get_label() for l in ln1+ln2])
        fig.tight_layout(rect=rect)
        if args.show: plt.show()
        else: fig.savefig(os.path.join(folder_path, 'kxx_kyy_vs_time_dualy.png'), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

        # vs SSA (kyy)
        valid = np.isfinite(rd.ssa) & np.isfinite(rd.k_yy)
        if np.count_nonzero(valid):
            fig, ax = plt.subplots(figsize=(11.28, 3.8))
            ax.plot(rd.ssa[valid], rd.k_yy[valid], 'o-', lw=1.5)
            ax.set_xlabel(r'SSA/SSA$_0$' if normalize else 'Specific Surface Area')
            ax.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize else r'$k_{yy}$ (W/m·K)')
            ax.set_title(f'T = {rd.temperature:g} K'); beautify_axes(ax); fig.tight_layout()
            if args.show: plt.show()
            else: fig.savefig(os.path.join(folder_path, 'kyy_vs_ssa.png'), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

        # vs SSA (kxx)
        valid = np.isfinite(rd.ssa) & np.isfinite(rd.k_xx)
        if np.count_nonzero(valid):
            fig, ax = plt.subplots(figsize=(11.28, 3.8))
            ax.plot(rd.ssa[valid], rd.k_xx[valid], 'o-', lw=1.5)
            ax.set_xlabel(r'SSA/SSA$_0$' if normalize else 'Specific Surface Area')
            ax.set_ylabel(r'$k_{xx}/k_{xx,0}$' if normalize else r'$k_{xx}$ (W/m·K)')
            ax.set_title(f'T = {rd.temperature:g} K'); beautify_axes(ax); fig.tight_layout()
            if args.show: plt.show()
            else: fig.savefig(os.path.join(folder_path, 'kxx_vs_ssa.png'), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

        # optional dual-x (time bottom, SSA top) for kyy & kxx
        if args.dualx:
            for comp, vals, fname in (('yy', rd.k_yy, 'kyy_time_ssa_dualx.png'),
                                      ('xx', rd.k_xx, 'kxx_time_ssa_dualx.png')):
                fig, ax = plt.subplots(figsize=(11.28, 4.0))
                me = sparse_markevery(rd.time_scaled)
                ax.plot(rd.time_scaled, vals, lw=2.0, marker='o', markevery=me, ms=4)
                ax.set_xlabel(f'Time ({rd.time_unit})')
                ax.set_ylabel(rf'$k_{comp}/k_{{{comp},0}}$' if normalize else rf'$k_{comp}$ (W/m·K)')
                ax.set_title(f'T = {rd.temperature:g} K'); beautify_axes(ax)
                ok,_ = add_secondary_ssa_axis(ax, rd.time, rd.ssa, normalize)
                if not ok: print(f"[warn] SSA not strictly monotonic in {rd.folder}")
                fig.tight_layout()
                if args.show: plt.show()
                else: fig.savefig(os.path.join(folder_path, fname), dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

if __name__ == "__main__":
    main()