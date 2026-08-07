# preprocess/scratch — quarantine, not an archive

Everything here was in `preprocess/` and is no longer used. It is kept only so
a function can be pulled back out if it turns out to be needed.

**This directory will be cleared out wholesale.** Do not add anything here that
you want to keep, and do not treat its contents as reference material.

To bring a file back: `git mv preprocess/scratch/<file> preprocess/<file>`.

What stayed in `preprocess/`:

- `comp_eps.py` — sizing `eps` from Kaempfer & Plapp sharp-interface bounds.
  Required by `CLAUDE.md` for every geometry.
- `build_geometry_*.py` and `igakit_complexgeocodes/` — the igakit complex-mesh
  builders. `build_geometry_ripening.py` and `build_geometry_two_throat.py` do
  not import igakit directly but write their meshes through
  `build_geometry_multi_grain`; `src/lunar_main.c` names the latter as the
  authoritative generator for the two-throat opts.

This project does **not** get a grain-packing generator — packed-grain
geometries are an `enceladus_DSM` concern only. `generateCircularGrains.py`,
the one packing script that was here, is in this folder and is not coming back.

The loose `*.png` files are geometry previews written next to the source by the
`build_geometry_*.py` scripts. They are build products, not inputs.
