#!/usr/bin/env python3
"""make_run_previews.py: one preview PDF per planned run, sorted into
preprocess/previews/type<N>/.

Reads the batch tests files (studies/cryosuction/type<N>_tests.txt: lines
"geometry:experiment   # ID") and renders each run with plot_geometry.render:
the domain and its initial ice, every element line of the mesh, and a header
with the run ID, the geometry and the experiment's conditions (thermal set,
vapor set, wall angle, length). Files are named "<ID>_<geometry>.pdf".

    venv_cryosuction/bin/python preprocess/make_run_previews.py [--types 1,2,3,4,5,6,7] [--jobs 6]

The PDFs are generated files (git-ignored); rerun after changing inputs.
"""
from __future__ import annotations

import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "preprocess"))
STUDY = ROOT / "studies/cryosuction"
OUT = ROOT / "preprocess/previews"


def find_opts(base: Path, name: str) -> Path:
    hits = list(base.rglob(f"{name}.opts"))
    hits = [h for h in hits if "scratch" not in h.parts]
    if len(hits) != 1:
        raise SystemExit(f"{name}: {len(hits)} matches under {base}")
    return hits[0]


def runs(t: str):
    f = STUDY / f"type{t}_tests.txt"
    for line in f.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        body, _, rid = line.partition("#")
        geom, exp = body.split()[0].split(":")
        rid = rid.strip().split()[0] if rid.strip() else f"{geom}__{exp}"
        yield t, re.sub(r"[^A-Za-z0-9_.-]", "_", rid), geom, exp


def one(job):
    t, rid, geom, exp = job
    from plot_geometry import render
    g = find_opts(ROOT / "inputs/geometry", geom)
    e = find_opts(ROOT / "inputs/experiment", exp)
    out = OUT / f"type{t}" / f"{rid}_{geom}.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    bad = render(g, out, title=f"run {rid}: {geom}", exp_path=e)
    return out, bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--types", default="1,2,3,4,5,6,7")
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args(argv)
    jobs = [j for t in a.types.split(",") if (STUDY / f"type{t}_tests.txt").exists() for j in runs(t)]
    missing = [t for t in a.types.split(",") if not (STUDY / f"type{t}_tests.txt").exists()]
    if missing:
        print(f"no tests file for type(s) {', '.join(missing)}; skipped")
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        res = list(ex.map(one, jobs))
    bad = [(o, b) for o, b in res if b]
    print(f"{len(res)} previews in {OUT.relative_to(ROOT)}; figcheck problems: {len(bad)}")
    for o, b in bad:
        print(f"  {o.relative_to(ROOT)}: {b}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
