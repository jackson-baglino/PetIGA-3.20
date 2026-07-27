#!/usr/bin/env python3
"""
make_regolith_texture.py — generate a tileable lunar-regolith-colored grain
texture for use as the sediment-region background in make_movie.py.

Run once with this project's regular Python (matplotlib + numpy, no
ParaView needed) -- NOT pvpython:

    python3 postprocess/make_regolith_texture.py

Writes postprocess/textures/lunar_regolith.png, checked into git so
make_movie.py's default --sediment-texture doesn't require regenerating
it. Re-run with --seed/--n-grains/--size to tweak the look.
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Warm gray-tan base tone, in the range commonly used to depict lunar
# regolith (Apollo soil samples: medium-to-dark gray with a brownish cast).
DEFAULT_BASE_RGB = (120, 110, 96)


def make_texture(out_path, base_rgb=DEFAULT_BASE_RGB, n_grains=2200,
                  size_px=1200, dpi=200, seed=42):
    rng = np.random.default_rng(seed)
    base = np.array(base_rgb) / 255.0

    fig, ax = plt.subplots(figsize=(size_px / dpi, size_px / dpi))
    ax.set_facecolor(base)
    fig.patch.set_facecolor(base)

    xs = rng.uniform(0, 1, n_grains)
    ys = rng.uniform(0, 1, n_grains)
    sizes = rng.uniform(0.0008, 0.004, n_grains)
    shade = rng.normal(0, 0.18, n_grains)  # multiplicative brightness jitter

    for x, y, s, sh in zip(xs, ys, sizes, shade):
        c = np.clip(base * (1 + sh), 0, 1)
        ax.add_patch(Circle((x, y), s, color=c, alpha=0.55, linewidth=0))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout(pad=0)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, facecolor=base)
    plt.close(fig)
    print(f"wrote {out_path}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="postprocess/textures/lunar_regolith.png")
    p.add_argument("--base-rgb", type=int, nargs=3, default=DEFAULT_BASE_RGB,
                    help="base color as R G B in 0-255 (default: warm gray-tan)")
    p.add_argument("--n-grains", type=int, default=2200)
    p.add_argument("--size-px", type=int, default=1200)
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    make_texture(Path(args.out), tuple(args.base_rgb), args.n_grains,
                 args.size_px, seed=args.seed)
