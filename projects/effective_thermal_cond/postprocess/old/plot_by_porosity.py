# plot_by_porosity.py
import os, argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from dsm_plotlib import (set_default_rc, collect_runs, color_for, beautify_axes)

def main():
    ap = argparse.ArgumentParser(description="Group plots by porosity (and seeds)")
    ap.add_argument("parent_dir")
    ap.add_argument("--pattern", default=None)
    ap.add_argument("--no-normalize", action="store_true")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--prefer-seed", default="7", help="Preferred seed when plotting across porosities")
    args = ap.parse_args()
    normalize = not args.no_normalize

    set_default_rc()
    runs = collect_runs(args.parent_dir, args.pattern or None, normalize=normalize)
    if not runs:
        print("No runs found."); return

    # by porosity (round to 3 decimals)
    by_por = defaultdict(list)
    for d in runs:
        por = d.porosity
        if por is not None and np.isfinite(por): by_por[round(por,3)].append(d)

    # seeds within each porosity
    for por, group in by_por.items():
        if len(group) < 2: continue
        fig, ax = plt.subplots(figsize=(11.28,4.2))
        for i, d in enumerate(group):
            col = color_for(i)
            t = d.time_scaled
            valid = np.isfinite(t) & np.isfinite(d.k_yy)
            if np.count_nonzero(valid):
                seed = d.seed if d.seed is not None else f"run{i+1}"
                ax.plot(t[valid], d.k_yy[valid], lw=2.0, marker='o', ms=4, markevery=2, color=col, label=f"seed={seed}")
        ax.set_xlabel(f"Time ({group[0].time_unit})")
        ax.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize else r'$k_{yy}$ (W/m·K)')
        ax.set_title(f"Porosity φ={por:.3f}: k_yy vs Time (seeds)")
        beautify_axes(ax); ax.legend(); fig.tight_layout()
        out = os.path.join(args.parent_dir, f'kyy_vs_time_por{por:.3f}_seeds.png')
        fig.savefig(out, dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

    # across porosities: pick preferred seed if available
    if len(by_por) > 1:
        fig, ax = plt.subplots(figsize=(11.28,4.2))
        for i, (por, group) in enumerate(sorted(by_por.items())):
            chosen = None
            for d in group:
                sd = str(d.seed).strip() if d.seed is not None else ""
                if sd == args.prefer_seed: chosen = d; break
            if chosen is None: chosen = group[0]
            t = chosen.time_scaled; col = color_for(i)
            valid = np.isfinite(t) & np.isfinite(chosen.k_yy)
            if np.count_nonzero(valid):
                lab = f"φ={por:.3f}, seed={chosen.seed}" if chosen.seed is not None else f"φ={por:.3f}"
                ax.plot(t[valid], chosen.k_yy[valid], lw=2.0, marker='o', ms=4, markevery=2, color=col, label=lab)
        ax.set_xlabel(f"Time ({runs[0].time_unit})")
        ax.set_ylabel(r'$k_{yy}/k_{yy,0}$' if normalize else r'$k_{yy}$ (W/m·K)')
        ax.set_title("k_yy vs Time across porosities")
        beautify_axes(ax); ax.legend(); fig.tight_layout()
        out = os.path.join(args.parent_dir, 'kyy_vs_time_multiple_porosities.png')
        fig.savefig(out, dpi=args.dpi, bbox_inches='tight', transparent=True); plt.close(fig)

if __name__ == "__main__":
    main()