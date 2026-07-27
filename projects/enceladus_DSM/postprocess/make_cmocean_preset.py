#!/usr/bin/env python3
"""
make_cmocean_preset.py — export a cmocean colormap as a ParaView JSON color
preset, so pvpython scripts (which run under ParaView's own embedded Python,
not this project's venv) can apply it without importing cmocean directly.

Run once with this project's regular Python (the one with cmocean installed,
e.g. venv_pf311) -- NOT pvpython:

    python3 postprocess/make_cmocean_preset.py ice
    python3 postprocess/make_cmocean_preset.py --out postprocess/colormaps/cmocean_dense.json dense

The output JSON goes in postprocess/colormaps/ and is checked into git --
make_movie.py (run via pvpython) loads it with servermanager.LoadPreset()
without ever needing cmocean itself.
"""

import argparse
import json
from pathlib import Path

import cmocean
import matplotlib.cm as mcm


def export_preset(cmap_name: str, n_samples: int, out_path: Path) -> None:
    try:
        cmap = getattr(cmocean.cm, cmap_name)
    except AttributeError:
        cmap = mcm.get_cmap(cmap_name)

    rgb_points = []
    for i in range(n_samples):
        t = i / (n_samples - 1)
        r, g, b, _a = cmap(t)
        rgb_points += [t, r, g, b]

    preset = [{
        "ColorSpace": "RGB",
        "Name": f"cmocean_{cmap_name}",
        "NanColor": [0.0, 0.0, 0.0],
        "RGBPoints": rgb_points,
    }]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(preset, fh, indent=2)
    print(f"wrote {out_path} (preset name: cmocean_{cmap_name})")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cmap", help="cmocean colormap name, e.g. ice, dense, haline, thermal")
    p.add_argument("--n-samples", type=int, default=64,
                   help="number of RGB control points sampled from the colormap (default: 64)")
    p.add_argument("--out", default=None,
                   help="output JSON path (default: postprocess/colormaps/cmocean_<name>.json)")
    args = p.parse_args()

    out_path = Path(args.out) if args.out else \
        Path(__file__).parent / "colormaps" / f"cmocean_{args.cmap}.json"
    export_preset(args.cmap, args.n_samples, out_path)


if __name__ == "__main__":
    main()
