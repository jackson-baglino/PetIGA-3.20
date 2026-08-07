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
  builders. `build_geometry_ripening.py` does not import igakit itself but
  writes its mesh through `build_geometry_multi_grain`.
- `packing_lib.py` / `generate_packing_gravity.py` — the current grain-packing
  generator. It replaces `generate_packing.py`, `generateCircularGrains.py`,
  `FCC_grain_packing.py`, `circle_packing_touching.py` and
  `select_packing_seeds.py`, all of which are in here.
