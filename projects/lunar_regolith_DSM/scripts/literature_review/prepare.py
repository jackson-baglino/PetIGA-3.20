#!/usr/bin/env python3
"""P1, first part: page images, tiles, parsed text, equations, sections and page map (WORKFLOW.md §3, §6).

For one copy of a paper (a PDF), writes into the run folder

    <run>/untrusted/content/<paper_id>/<copy_hash>/
      images/p001.png ...    300 dpi page images
      layout.json            the parser's blocks per page (type and box), used to cut the tiles
      tiles/p001_t1.png ...  tiles: whole layout blocks, about 4.5 inches of a column, at most
                             2,000 px on the long edge (Claude's cap when a request holds more
                             than 20 images); tiles.json gives each tile's box on the page
      pages/p001.md ...      the parser's text of each page (PDF page numbers)
      equations.jsonl        displayed equations the parser found: LaTeX, printed number, page
      sections.json          the parser's headings with their pages
      page_map.json          PDF page to printed page (from the printed page numbers)
      manifest.json          every file's sha256; tool versions, flags and model files

and a summary with counts only, outside untrusted/, in <run>/summary_prepare_<paper_id>.json.
Everything under untrusted/ holds paper text and is denied to sessions; this script never
prints paper text, only counts, page numbers and hashes.

Not built yet (P1 later parts): text from arXiv LaTeX or JATS source, tables (two
extractions compared), figure crops, disagreements between sources, and publishing into
the library through curator.py.

Tiles. With the parser's layout, tiles are made of whole layout blocks: a
table, or a figure with its axis labels and caption, is one tile, scaled down if it is over
2,000 px; text and equations are packed in reading order into tiles of about 4.5 inches of
a column, so cuts fall between blocks. Without a parser (--parser none), or on a page with no
blocks, the older ink-only method below is used.

Ink-only tiles. A page is split into columns at a gutter: a band of blank pixel columns near the
middle that the ink crosses on few rows. Rows where ink crosses the gutter (titles,
wide figures and equations) are full-width bands; the rest is cut per column. Each band is
cut into tiles of about half a column high at blank rows, so no line is cut in two; where
no blank row lies near the target, the next tile starts 1 inch above the cut. A tile wider
than 2,000 px is scaled down to fit, and its dpi is recorded.

Parser. MinerU2.5-Pro (docparse_run.py, in venv_docparse) reads the 300 dpi page images
through MLX on this machine: a layout pass, then each block read as LaTeX, HTML or text.
Chosen over Marker on 2026-10-29: it read all 10 equations of Kaempfer & Plapp p. 7 exactly
as the user transcribed them (Marker 7) and finds the printed equation numbers as blocks.
Marker (marker_run.py, `--parser marker`) is kept as a fallback, never with --use_llm. Both
run offline from local weights; nothing leaves the machine.

    venv_workflows/bin/python scripts/literature_review/prepare.py --run <run folder> \\
        --paper-id <id> --pdf <file> [--parser mineru|marker|none] [--allow-dirty]
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_calibration as cc  # noqa: E402

ROOT = Path(cc.ROOT)
DPI = 300
MAX_EDGE = 2000                 # px; Claude scales larger images down when a request holds > 20
OVERLAP = DPI                   # 1 inch, where a tile must be cut without a gap
INK = 200                       # gray level below which a pixel is ink
MIN_GAP = 6                     # px of blank rows that count as a gap between lines (0.5 pt)
GUTTER = int(0.12 * DPI)        # a column gutter is at least this wide
MARKER_FLAGS = ["--paginate_output", "--output_format", "markdown", "--disable_image_extraction",
                "--disable_tqdm"]
PAGE_SEP = re.compile(r"^\{(\d+)\}-{48}\s*$", flags=re.M)
MODEL_CACHE = Path.home() / "Library" / "Caches" / "datalab" / "models"     # Marker's weights
DOCPARSE_PY = Path(__file__).resolve().parents[2] / "venv_docparse" / "bin" / "python"
DEFINITIONS = ["scripts/literature_review", "workflows/literature-review", ".claude/permissions"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def uncommitted_definitions():
    out = subprocess.run(["git", "status", "--porcelain", "--", *DEFINITIONS], cwd=ROOT,
                         capture_output=True, text=True).stdout
    return [l[3:] for l in out.splitlines() if l.strip()]


# ---- page images ----

def render(pdf, outdir):
    """300 dpi PNG per page, named p001.png; returns the number of pages."""
    outdir.mkdir(parents=True, exist_ok=True)
    subprocess.run(["mutool", "draw", "-q", "-r", str(DPI), "-c", "gray", "-o",
                    str(outdir / "p%03d.png"), str(pdf)], check=True)
    return len(list(outdir.glob("p*.png")))


# ---- tiles ----

def runs(flags):
    """[(start, end)) of consecutive True values."""
    out, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            out.append((start, i))
            start = None
    if start is not None:
        out.append((start, len(flags)))
    return out


def find_gutter(ink, x0, x1):
    """(g0, g1) of a column gutter near the middle of [x0, x1), or None.

    A gutter is a band at least GUTTER wide that contains the center (within 10 % of the
    width), that ink crosses on at most 5 % of the inked rows, and that has ink on both
    sides on at least 30 % of the inked rows each (text in two columns).
    """
    rows = ink[:, x0:x1].any(axis=1)
    if rows.sum() == 0:
        return None
    sub = ink[rows][:, x0:x1]
    crossed = sub.sum(axis=0) / sub.shape[0]                # share of inked rows that cross each column
    w = x1 - x0
    c = w // 2
    best = None
    for a, b in runs(list(crossed <= 0.05)):
        if b - a < GUTTER or not (a - w // 10 <= c <= b + w // 10):
            continue
        left = sub[:, :a].any(axis=1).mean()
        right = sub[:, b:].any(axis=1).mean()
        if left >= 0.3 and right >= 0.3 and (best is None or b - a > best[1] - best[0]):
            best = (a + x0, b + x0)
    return best


def cut_band(ink, y0, y1, x0, x1, target):
    """Vertical cuts of the band [y0, y1) x [x0, x1) into pieces of about target rows.

    Returns [(top, bottom)]. Cuts fall in blank rows (no ink in the band's columns); where
    there is none between half and the full maximum height, the piece ends at the maximum
    height and the next starts OVERLAP rows above it.
    """
    blank = ~ink[y0:y1, x0:x1].any(axis=1)
    gaps = [((a + b) // 2 + y0) for a, b in runs(list(blank)) if b - a >= MIN_GAP]
    limit = MAX_EDGE
    pieces, top = [], y0
    while y1 - top > limit or (y1 - top > target * 1.5):
        want = top + target
        cands = [g for g in gaps if top + target // 2 <= g <= top + limit]
        if cands:
            cut = min(cands, key=lambda g: abs(g - want))
            pieces.append((top, cut))
            top = cut
        elif y1 - top > limit:
            cut = top + limit
            pieces.append((top, cut))
            top = cut - OVERLAP
        else:
            break                                         # no gap, but what is left fits in one tile
    pieces.append((top, y1))
    return pieces


def plan_tiles(gray):
    """Tile boxes (x0, y0, x1, y1) in pixels covering all the ink of a page image."""
    ink = gray < INK
    ys = np.where(ink.any(axis=1))[0]
    xs = np.where(ink.any(axis=0))[0]
    if len(ys) == 0:
        return []
    pad = 4
    Y0, Y1 = max(0, ys[0] - pad), min(gray.shape[0], ys[-1] + 1 + pad)
    X0, X1 = max(0, xs[0] - pad), min(gray.shape[1], xs[-1] + 1 + pad)
    g = find_gutter(ink, X0, X1)
    if g is None:
        target = max(1, (Y1 - Y0) // 2)                     # single column: half a page
        return [(X0, a, X1, b) for a, b in cut_band(ink, Y0, Y1, X0, X1, target)]
    mid = (g[0] + g[1]) // 2
    # Blocks: runs of inked rows across the page, split at blank rows of at least MIN_GAP.
    # A block is wide if ink crosses the gutter in it; neighboring blocks of the same kind
    # form one band, and bands meet in the middle of the blank rows between them.
    rows = ink[Y0:Y1, X0:X1].any(axis=1)
    blocks = [(a + Y0, b + Y0) for a, b in runs(list(rows))]
    joined = []
    for a, b in blocks:
        if joined and a - joined[-1][1] < MIN_GAP:
            joined[-1] = (joined[-1][0], b)
        else:
            joined.append((a, b))
    third = max(1, (g[1] - g[0]) // 3)                      # the gutter's middle third, so a line
    wide = ink[:, g[0] + third:g[1] - third].any(axis=1)    # ending at its edge is not wide
    bands = []
    for a, b in joined:
        kind = "wide" if wide[a:b].any() else "cols"
        if bands and bands[-1][0] == kind:
            bands[-1][2] = b
        else:
            bands.append([kind, a, b])
    for i in range(len(bands) - 1):                         # meet in the middle of the gap
        m = (bands[i][2] + bands[i + 1][1]) // 2
        bands[i][2], bands[i + 1][1] = m, m
    bands[0][1], bands[-1][2] = Y0, Y1
    col_h = sum(b - a for k, a, b in bands if k == "cols")
    target = max(1, col_h // 2) if col_h else max(1, (Y1 - Y0) // 2)
    boxes = []
    for kind, a, b in bands:
        if kind == "wide":
            boxes += [(X0, p, X1, q) for p, q in cut_band(ink, a, b, X0, X1, target)]
        else:
            for c0, c1 in ((X0, mid), (mid, X1)):
                if ink[a:b, c0:c1].any():
                    boxes += [(c0, p, c1, q) for p, q in cut_band(ink, a, b, c0, c1, target)]
    return boxes


KEEP_WHOLE = {"Table", "TableGroup", "Figure", "FigureGroup", "Picture", "PictureGroup", "Form"}
FLOATS = KEEP_WHOLE - {"Form"}
TARGET = DPI * 9 // 2           # px; a tile holds about 4.5 inches of a column, as before
NEAR = DPI // 4                 # px; ink outside every block joins a block this close
PAD = 6                         # px around each block


def _union(a, b):
    return (min(a[0], b[0]), min(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3]))


def _gap(a, b):
    """Distance in px between two boxes (0 if they touch or overlap)."""
    dx = max(0, a[0] - b[2], b[0] - a[2])
    dy = max(0, a[1] - b[3], b[1] - a[3])
    return max(dx, dy)


def _area(a):
    return (a[2] - a[0]) * (a[3] - a[1])


def _overlap(a, b):
    return max(0, min(a[2], b[2]) - max(a[0], b[0])) * max(0, min(a[3], b[3]) - max(a[1], b[1]))


def layout_atoms(gray, layout):
    """The parser's blocks on one page as [type, box] in pixels, in reading order.

    Ink outside every block (axis labels, stray numbers) joins the nearest block within
    NEAR, or becomes its own atom; a caption outside a figure group joins the figure or
    table next to it in reading order.
    """
    from scipy import ndimage
    H, W = gray.shape
    sx, sy = W / layout["bbox"][2], H / layout["bbox"][3]
    atoms = []
    for b in layout["blocks"]:
        x0, y0, x1, y1 = b["bbox"]
        box = (max(0, int(x0 * sx) - PAD), max(0, int(y0 * sy) - PAD),
               min(W, int(np.ceil(x1 * sx)) + PAD), min(H, int(np.ceil(y1 * sy)) + PAD))
        if box[2] > box[0] and box[3] > box[1]:
            atoms.append([b["type"], box])
    ink = gray < INK
    inside = np.zeros_like(ink)
    for _, (x0, y0, x1, y1) in atoms:
        inside[y0:y1, x0:x1] = True
    stray = ink & ~inside
    if stray.any():
        lab, n = ndimage.label(ndimage.binary_dilation(stray, iterations=MIN_GAP))
        for sl in ndimage.find_objects(lab):
            ys, xs = sl
            sub = stray[sl] & (lab[sl] > 0)
            if not sub.any():
                continue
            box = (xs.start, ys.start, xs.stop, ys.stop)
            near = [(i, _gap(box, a[1])) for i, a in enumerate(atoms)]
            i, d = min(near, key=lambda t: t[1]) if near else (None, None)
            # Join only what barely grows the block (axis labels under a plot); a margin
            # stamp running down the page would otherwise swallow the page.
            if i is not None and d <= NEAR and _area(_union(atoms[i][1], box)) <= 1.5 * _area(atoms[i][1]):
                atoms[i][1] = _union(atoms[i][1], box)
            else:
                atoms.append(["Ink", box])
    inside = [i for i, a in enumerate(atoms)                  # blocks nested in another block
              if any(j != i and _overlap(a[1], b[1]) >= 0.9 * _area(a[1])
                     and (_area(b[1]) > _area(a[1]) or j < i) for j, b in enumerate(atoms))]
    atoms = [a for i, a in enumerate(atoms) if i not in inside]
    # Page furniture (running heads, page numbers, rules) in the top or bottom 15 % of the
    # page becomes one header and one footer atom, unless that would overlap the body.
    for kind, band in (("PageHeader", lambda b: b[3] <= 0.15 * H), ("PageFooter", lambda b: b[1] >= 0.85 * H)):
        group = [i for i, (t, b) in enumerate(atoms) if t in ("PageHeader", "PageFooter", "Ink") and band(b)]
        if len(group) < 2:
            continue
        box = atoms[group[0]][1]
        for i in group[1:]:
            box = _union(box, atoms[i][1])
        if all(_overlap(box, b) == 0 for i, (t, b) in enumerate(atoms) if i not in group):
            rest = [a for i, a in enumerate(atoms) if i not in group]
            atoms = [[kind, box]] + rest if kind == "PageHeader" else rest + [[kind, box]]
    merged = []
    for t, box in atoms:
        if t == "Caption" and merged and merged[-1][0] in FLOATS:
            merged[-1][1] = _union(merged[-1][1], box)
        elif merged and merged[-1][0] == "Caption" and t in FLOATS:
            merged[-1] = [t, _union(merged[-1][1], box)]
        else:
            merged.append([t, box])
    return merged


def plan_layout_tiles(gray, layout):
    """Tile boxes from the parser's layout: whole blocks, packed in reading order.

    Consecutive blocks share a tile when they share columns, the tile stays under TARGET
    rows (or the blocks overlap, as parts of one figure do) and it overlaps no other block,
    which keeps a tile inside its column. Tables and figures (with their captions) are
    never cut; a text block taller than MAX_EDGE is cut at blank rows as in plan_tiles.
    """
    atoms = layout_atoms(gray, layout)
    ink = gray < INK
    tiles = []                                              # [box, [atom indices]]

    def fits(box, members):
        return all(_overlap(box, a[1]) <= 0.02 * (a[1][2] - a[1][0]) * (a[1][3] - a[1][1])
                   for j, a in enumerate(atoms) if j not in members)

    for i, (_t, box) in enumerate(atoms):
        if tiles:
            cur, members = tiles[-1]
            side = min(cur[2], box[2]) - max(cur[0], box[0]) > 0   # share some columns
            u = _union(cur, box)
            touching = _overlap(cur, box) > 0                      # parts of one figure overlap
            if side and (u[3] - u[1] <= TARGET or touching) and fits(u, members + [i]):
                tiles[-1] = [u, members + [i]]
                continue
        tiles.append([box, [i]])
    boxes = []
    for box, members in tiles:
        x0, y0, x1, y1 = box
        thin = (x1 - x0) * 20 < (y1 - y0)                    # a margin stamp: keep it in one piece
        if y1 - y0 > MAX_EDGE and not thin and not any(atoms[j][0] in KEEP_WHOLE for j in members):
            boxes += [(x0, a, x1, b) for a, b in cut_band(ink, y0, y1, x0, x1, TARGET)]
        else:
            boxes.append(box)
    return boxes


def write_tiles(img_path, outdir, page, layout=None):
    gray = np.asarray(Image.open(img_path).convert("L"))
    boxes = plan_layout_tiles(gray, layout) if layout and layout["blocks"] else plan_tiles(gray)
    tiles = []
    boxes = [b for b in boxes if b[2] > b[0] and b[3] > b[1]]
    for k, (x0, y0, x1, y1) in enumerate(boxes, 1):
        tile = Image.fromarray(gray[y0:y1, x0:x1])
        scale = min(1.0, MAX_EDGE / max(tile.size))
        if scale < 1.0:
            tile = tile.resize((max(1, round(tile.size[0] * scale)), max(1, round(tile.size[1] * scale))),
                               Image.LANCZOS)
        name = f"p{page:03d}_t{k}.png"
        tile.save(outdir / name, optimize=False)
        tiles.append({"file": name, "page": page, "box_px": [int(x0), int(y0), int(x1), int(y1)],
                      "size_px": list(tile.size), "dpi": round(DPI * scale, 1)})
    return tiles, gray


def coverage(gray, tiles):
    """Fraction of ink pixels inside some tile box."""
    ink = gray < INK
    covered = np.zeros_like(ink)
    for t in tiles:
        x0, y0, x1, y1 = t["box_px"]
        covered[y0:y1, x0:x1] = True
    n = ink.sum()
    return 1.0 if n == 0 else float((ink & covered).sum() / n)


# ---- parsers ----

def split_pages(md):
    """{page (1-based): markdown} from Marker's paginated output ({0}----... separators, 0-based)."""
    parts = PAGE_SEP.split(md)
    return {int(parts[i]) + 1: parts[i + 1].strip("\n") + "\n" for i in range(1, len(parts) - 1, 2)}


def run_marker(pdf):
    """Marker's pages {page: Markdown} and layout [{bbox, blocks}] (marker_run.py, one pass)."""
    env = dict(os.environ, HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1")
    with tempfile.TemporaryDirectory() as tmp:
        env["MPLCONFIGDIR"] = tmp
        p = subprocess.run([sys.executable, str(Path(__file__).with_name("marker_run.py")), str(pdf), tmp],
                           capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise RuntimeError(f"Marker failed with exit code {p.returncode}")
        return split_pages((Path(tmp) / "marker.md").read_text()), json.loads((Path(tmp) / "layout.json").read_text())


def run_mineru(images):
    """MinerU2.5-Pro's pages {page: Markdown}, layout [{bbox, blocks}] and tool record (docparse_run.py)."""
    env = dict(os.environ, HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1")
    with tempfile.TemporaryDirectory() as tmp:
        p = subprocess.run([str(DOCPARSE_PY), str(Path(__file__).with_name("docparse_run.py")), str(images), tmp, str(DPI)],
                           capture_output=True, text=True, env=env)
        if p.returncode != 0:
            raise RuntimeError(f"MinerU failed with exit code {p.returncode}")
        pages = {int(k): v for k, v in json.loads((Path(tmp) / "pages.json").read_text()).items()}
        return (pages, json.loads((Path(tmp) / "layout.json").read_text()),
                json.loads((Path(tmp) / "tool.json").read_text()))


DISPLAY = re.compile(r"\$\$(.+?)\$\$", flags=re.S)
EQ_NUMBER = r"\d+(?:\.\d+)*[a-z]?|[A-Z]\.?\d+(?:\.\d+)*[a-z]?"   # 12, 3b, 2.14, A1, A.3
TAG = re.compile(r"\\tag\{([^}]*)\}|\((" + EQ_NUMBER + r")\)\s*$")
HEADING = re.compile(r"^(#{1,6})\s+(.+?)\s*$", flags=re.M)


def equations(pages):
    out = []
    for p in sorted(pages):
        for m in DISPLAY.finditer(pages[p]):
            latex = m.group(1).strip()
            t = TAG.search(latex)
            # The number may follow the equation directly, or end a short line after it
            # that carries a qualifier, as in "in $\Omega$ , (2.21)".
            after = pages[p][m.end():m.end() + 80]
            t2 = (re.match(r"\s*\((" + EQ_NUMBER + r")\)", after)
                  or re.match(r"[ \t]*\n?(?:(?!\$\$)[^\n]){0,60}?\((" + EQ_NUMBER + r")\)[ \t]*(?:\n|$)", after))
            number = (t.group(1) or t.group(2)) if t else (t2.group(1) if t2 else None)
            out.append({"page_pdf": p, "number": number, "latex": latex})
    return out


def sections(pages):
    return [{"page_pdf": p, "level": len(m.group(1)), "heading": m.group(2)}
            for p in sorted(pages) for m in HEADING.finditer(pages[p])]


def marker_versions():
    v = {p: importlib.metadata.version(p) for p in ("marker-pdf", "surya-ocr", "torch")}
    v["model_files"] = {str(f.relative_to(MODEL_CACHE)): f.stat().st_size
                        for f in sorted(MODEL_CACHE.rglob("*")) if f.is_file()} if MODEL_CACHE.exists() else None
    return v


# ---- driver ----

def prepare(run, paper_id, pdf, parser="mineru"):
    copy_hash = sha256(pdf)
    final = Path(run) / "untrusted" / "content" / paper_id / copy_hash[:16]
    if final.exists():
        raise SystemExit(f"{final.relative_to(run)} exists; content is written once")
    out = final.with_name(final.name + ".partial")       # written here, renamed when complete, so a P1 stopped
    if out.exists():                                     # midway leaves no half-written content: the next run
        shutil.rmtree(out)                               # clears the partial folder and starts again
    t0 = time.time()
    n = render(pdf, out / "images")
    pages, layout, tool = {}, [], None
    if parser == "mineru":
        pages, layout, tool = run_mineru(out / "images")
    elif parser == "marker":
        pages, layout = run_marker(pdf)
        tool = {"flags": MARKER_FLAGS, **marker_versions()}
    if layout:
        (out / "layout.json").write_text(json.dumps(layout, indent=1))
    (out / "tiles").mkdir()
    all_tiles, failed, low_cover = [], [], []
    for p in range(1, n + 1):
        img = out / "images" / f"p{p:03d}.png"
        tiles, gray = write_tiles(img, out / "tiles", p, layout[p - 1] if p <= len(layout) else None)
        c = coverage(gray, tiles)
        if c < 0.999:
            low_cover.append({"page_pdf": p, "ink_covered": round(c, 4)})
        if any(max(t["size_px"]) > MAX_EDGE for t in tiles):
            failed.append({"page_pdf": p, "why": "tile over the size limit"})
        all_tiles += tiles
    (out / "tiles.json").write_text(json.dumps(all_tiles, indent=1))
    pmap = cc.PageMap(str(pdf))
    pm = {"offset": pmap.offset, "agreed_by": pmap.agree,
          "pages": {str(p): (p + pmap.offset if pmap.offset is not None else None) for p in range(1, n + 1)}}
    (out / "page_map.json").write_text(json.dumps(pm, indent=1))
    record = {"paper_id": paper_id, "copy_sha256": copy_hash, "copy_file": os.path.basename(pdf),
              "made": datetime.datetime.now().isoformat(timespec="seconds"), "dpi": DPI,
              "max_edge_px": MAX_EDGE, "pages": n,
              "tools": {"mutool": subprocess.run(["mutool", "-v"], capture_output=True, text=True)
                        .stderr.strip().splitlines()[0]}}
    eqs, secs, no_text = [], [], []
    if parser != "none":
        (out / "pages").mkdir()
        for p in range(1, n + 1):
            text = pages.get(p, "")
            if not text.strip():
                no_text.append(p)
            (out / "pages" / f"p{p:03d}.md").write_text(text)
        eqs, secs = equations(pages), sections(pages)
        with open(out / "equations.jsonl", "w") as f:
            for e in eqs:
                f.write(json.dumps(e) + "\n")
        (out / "sections.json").write_text(json.dumps(secs, indent=1))
        record["tools"]["parser"] = {"name": parser, **tool}
    failed += [{"page_pdf": p, "why": f"no text from {parser}"} for p in no_text]
    record["failed_pages"] = failed
    record["files"] = {str(f.relative_to(out)): sha256(f) for f in sorted(out.rglob("*"))
                       if f.is_file() and f.name != "manifest.json"}
    (out / "manifest.json").write_text(json.dumps(record, indent=1))
    out.rename(final)
    out = final
    summary = {"stage": "P1", "paper_id": paper_id, "copy_sha256": copy_hash, "pages": n,
               "tiles": len(all_tiles), "tiles_per_page_max": max(
                   (sum(1 for t in all_tiles if t["page"] == p) for p in range(1, n + 1)), default=0),
               "parser": parser, "tiles_from_layout": bool(layout),
               "tiles_scaled_down": sum(1 for t in all_tiles if t["dpi"] < DPI),
               "min_tile_dpi": min((t["dpi"] for t in all_tiles), default=None),
               "pages_ink_not_fully_covered": low_cover, "failed_pages": failed,
               "equations": len(eqs), "numbered_equations": sum(1 for e in eqs if e["number"]),
               "headings": len(secs), "page_map_offset": pmap.offset,
               "manifest_sha256": sha256(out / "manifest.json"), "seconds": round(time.time() - t0, 1)}
    (Path(run) / f"summary_prepare_{paper_id}.json").write_text(json.dumps(summary, indent=1))
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True, help="run folder, runs/literature-review/<review>/<run_id>")
    ap.add_argument("--paper-id", required=True)
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--parser", choices=["mineru", "marker", "none"], default="mineru",
                    help="mineru (default); marker (fallback); none: images and tiles only")
    ap.add_argument("--allow-dirty", action="store_true", help="development only")
    a = ap.parse_args(argv)
    print(f"P1 prepare: run={a.run} paper={a.paper_id} pdf={os.path.basename(a.pdf)} dpi={DPI} "
          f"max_edge={MAX_EDGE} parser={a.parser}")
    dirty = uncommitted_definitions()
    if dirty and not a.allow_dirty:
        raise SystemExit("workflow definitions have uncommitted changes; commit first:\n  " + "\n  ".join(dirty))
    s = prepare(a.run, a.paper_id, a.pdf, a.parser)
    print(json.dumps(s, indent=1))
    return 1 if s["failed_pages"] or s["pages_ink_not_fully_covered"] else 0


if __name__ == "__main__":
    sys.exit(main())
