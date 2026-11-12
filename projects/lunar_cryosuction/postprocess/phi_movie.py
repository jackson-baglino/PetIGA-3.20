"""Ice-field movie of one 2D run (or a restart chain): one frame per saved
snapshot, phi on a fixed 0-1 color scale, phi = 1/2 contour, the domain
boundary rows of the mesh drawn as walls, and the time in days from the
snapshot time (pplib.snapshot_times: outp.txt, then SSA_evo.dat), not from the
frame index.

Plain matplotlib + ffmpeg; no ParaView needed (make_movie.py is the ParaView
version with the vapor field). Snapshots are not evenly spaced in time (steps
0 and 1 are both at t ~ 0; later ones follow -n_out), so every frame is held
for the same playback time and the label gives the true time.

Repeated frames are reported, not dropped: a snapshot whose phi equals the
previous one to machine precision is printed as IDENTICAL.

    venv_cryosuction/bin/python postprocess/phi_movie.py <run> [<run leg 2> ...] --out movie.mp4 \
        [--fps 4] [--frames-dir DIR] [--title TEXT]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib  # noqa: E402

UM = 1e6
WALL = "#8c7350"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run folder(s); later legs of a restart chain last")
    ap.add_argument("--out", required=True, help=".mp4 (ffmpeg) or .gif")
    ap.add_argument("--fps", type=float, default=4.0)
    ap.add_argument("--frames-dir", default=None, help="keep the PNG frames here (default: temporary)")
    ap.add_argument("--title", default=None)
    ap.add_argument("--width-in", type=float, default=8.0)
    ap.add_argument("--dpi", type=int, default=200)
    a = ap.parse_args()

    out = Path(a.out)
    fdir = Path(a.frames_dir) if a.frames_dir else Path(tempfile.mkdtemp(prefix="phi_movie_"))
    fdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "Arial", "font.size": 9})
    title = a.title or Path(a.runs[-1]).name
    prev, n, rows = None, 0, []
    for s in pplib.iter_snapshots(a.runs):
        X, Y, phi = s.X * UM, s.Y * UM, s.phi
        if np.ptp(Y[0, :]) > np.ptp(Y[:, 0]):    # along-wall index second
            X, Y, phi = X.T, Y.T, phi.T
        same = prev is not None and prev.shape == phi.shape and np.array_equal(prev, phi)
        dmax = float(np.abs(phi - prev).max()) if prev is not None and prev.shape == phi.shape else float("nan")
        rows.append((n, s.step, s.t / 86400, dmax, same))
        prev = phi
        lx, ly = X.max() - X.min(), Y.max() - Y.min()
        h = a.width_in * ly / lx + 0.9
        fig, ax = plt.subplots(figsize=(a.width_in, h))
        fig.subplots_adjust(left=0.06, right=0.90, top=1 - 0.42 / h, bottom=0.36 / h)
        pc = ax.pcolormesh(X, Y, phi, cmap="Blues", vmin=0, vmax=1, shading="gouraud", rasterized=True)
        ax.contour(X, Y, phi, levels=[0.5], colors=["#1f4e8c"], linewidths=0.7)
        for j in (0, -1):
            ax.plot(X[j], Y[j], color=WALL, lw=1.6)
        ax.set_aspect("equal")
        ax.set_xlabel("x [µm]")
        ax.set_yticks([])
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.set_title(f"{title}\n t = {s.t / 86400:6.2f} d   (step {s.step})", fontsize=8, loc="left")
        cax = fig.add_axes((0.92, 0.36 / h + 0.05, 0.015, 1 - 0.78 / h - 0.1))
        fig.colorbar(pc, cax=cax, ticks=[0, 0.5, 1]).set_label(r"$\phi$")
        fig.savefig(fdir / f"frame_{n:05d}.png", dpi=a.dpi)
        plt.close(fig)
        n += 1
    print(f"{n} frames in {fdir}")
    for k, st, td, dmax, same in rows:
        flag = "  IDENTICAL to previous" if same else ""
        print(f"  frame {k:3d}  step {st:6d}  t = {td:8.3f} d  max|dphi| = {dmax:.3e}{flag}")
    if out.suffix == ".gif" or not shutil.which("ffmpeg"):
        from PIL import Image
        ims = [Image.open(f) for f in sorted(fdir.glob("frame_*.png"))]
        gif = out.with_suffix(".gif")
        ims[0].save(gif, save_all=True, append_images=ims[1:], duration=int(1000 / a.fps), loop=0)
        print("wrote", gif)
    else:
        cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", str(fdir / "frame_%05d.png"),
               "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white", "-c:v", "libx264", "-pix_fmt", "yuv420p",
               "-crf", "18", str(out)]
        subprocess.run(cmd, check=True)
        print("wrote", out)
    if not a.frames_dir:
        shutil.rmtree(fdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
