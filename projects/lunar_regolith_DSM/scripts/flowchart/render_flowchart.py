#!/usr/bin/env python3
"""Render a workflow flowchart (HTML with inline SVG) from a JSON spec, then lint it.

    python scripts/flowchart/render_flowchart.py workflows/<name>/flowchart.json [--png]
    python scripts/flowchart/render_flowchart.py --all [--png]      # every workflows/*/flowchart.json

Layout is deterministic: each node sits on a grid (its lane is the column, its
"row" the row); Graphviz (neato -n2, orthogonal splines) routes the edges. Edge
labels are placed by this script, not by Graphviz: each label is measured with
the Helvetica font metrics and put at the first candidate spot beside its edge
that clears every node, edge, arrowhead, lane title and earlier label by a
margin. The agent writing the spec never places coordinates.

The lint then checks the rendered SVG: label/node, label/edge (any edge,
including the label's own), label/label and label/lane-title overlaps, and
edges that pass through a node that is not their end. Segments are tested
exactly (segment-rectangle intersection), not by sampling. Exit status 1 if
anything overlaps. --png also writes a screenshot (headless Chrome) to the
system temporary folder for visual review; it is never saved in the repo.

Colors, type and spacing come from scripts/flowchart/theme.css, which is
inlined into every page, so all workflow pages share one theme.

Spec (see scripts/flowchart/README.md for the full schema):
  {"title", "summary", "lanes": [{"id","label"}],
   "nodes": [{"id","lane","row","kind","title","job","model","tools","details"}],
   "edges": [{"from","to","label","kind"}]}
kind (node): user | script | agent | record.  kind (edge): flow | stop | loop | data.
"""
import html
import json
import re
import subprocess
import sys
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
THEME_CSS = HERE / "theme.css"

NODE_W = 290.0          # node width [pt]
LANE_W = 330.0          # lane width [pt]
ROW_GAP = 64.0          # vertical gap between rows [pt]
EDGE_PAD = 18.0         # gap above the first and below the last row [pt]
TOP = 50.0              # space for lane titles [pt]
WRAP = 46               # characters per line inside a node
LINE_H = 14.0           # line height estimate [pt]
LABEL_FS = 10.0         # edge label font size [pt]
LABEL_GAP = 5.0         # gap between a label and its own edge [pt]
CLEAR = 4.0             # clearance a label keeps from other things [pt]

KIND_TEXT = {"user": "user decision point", "script": "deterministic script: same input, same output",
             "agent": "agent (model and tools listed in the box)", "record": "record written by a stage"}
EDGE_TEXT = {"flow": "main flow", "data": "writes a record", "stop": "stop and ask the user",
             "loop": "loop back for a recheck"}
EDGE_W = {"flow": 1.5, "data": 1.0, "stop": 1.5, "loop": 1.4}


def load_theme(path=THEME_CSS):
    """CSS custom properties of the theme (first :root block) as a dict."""
    css = path.read_text()
    root = re.search(r":root\s*\{(.*?)\}", css, re.S).group(1)
    return css, {k: v.strip() for k, v in re.findall(r"--([\w-]+)\s*:\s*([^;]+);", root)}


THEME_TEXT, T = load_theme()
KINDS = ("user", "script", "agent", "record")
EDGES = ("flow", "data", "stop", "loop")


# ---------- text measurement ----------
_FONTS = {}


def _font(bold):
    if bold not in _FONTS:
        try:
            from PIL import ImageFont
            ttc = "/System/Library/Fonts/Helvetica.ttc"
            for i in range(6):
                f = ImageFont.truetype(ttc, 100, index=i)
                name = " ".join(f.getname()).lower()
                if ("bold" in name) == bold and "oblique" not in name and "light" not in name:
                    _FONTS[bold] = f
                    break
        except Exception:
            _FONTS[bold] = None
    return _FONTS.get(bold)


def text_w(s, fs, bold=False):
    """Rendered width of s in Helvetica at font size fs [pt]; a safe estimate if the font is unavailable."""
    f = _font(bold)
    if f is None:
        return 0.62 * fs * len(s)
    return f.getlength(s) * fs / 100.0


def esc(s):
    return html.escape(str(s), quote=False)


def wrap(s, width=WRAP):
    return textwrap.wrap(str(s), width) or [""]


def node_label(n):
    """Graphviz HTML-like label: title, model/tools line, one-sentence job, then details."""
    accent = T[f"kind-{n['kind']}"]
    rows = [f'<TR><TD ALIGN="LEFT"><FONT COLOR="{T["text"]}"><B>{esc(n["title"])}</B></FONT></TD></TR>']
    lines = 1
    if n["kind"] == "agent":
        meta = f'{n["model"]} · tools: {", ".join(n["tools"]) or "none"}'
        for l in wrap(meta):
            rows.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="10" COLOR="{accent}">{esc(l)}</FONT></TD></TR>')
            lines += 1
    for l in wrap(n.get("job", "")):
        rows.append(f'<TR><TD ALIGN="LEFT"><FONT COLOR="{T["text-2"]}"><I>{esc(l)}</I></FONT></TD></TR>')
        lines += 1
    for d in n.get("details", []):
        for i, l in enumerate(wrap(d, WRAP - 2)):
            rows.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="10" COLOR="{T["text-muted"]}">'
                        f'{"• " if i == 0 else "  "}{esc(l)}</FONT></TD></TR>')
            lines += 1
    label = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="1">' + "".join(rows) + "</TABLE>>"
    return label, lines * LINE_H + 22.0


def validate(spec):
    errs = []
    lanes = {l["id"] for l in spec["lanes"]}
    ids = set()
    for n in spec["nodes"]:
        for k in ("id", "lane", "row", "kind", "title", "job"):
            if k not in n:
                errs.append(f'node {n.get("id", "?")}: missing "{k}"')
        if n.get("lane") not in lanes:
            errs.append(f'node {n.get("id")}: unknown lane {n.get("lane")}')
        if n.get("kind") not in KINDS:
            errs.append(f'node {n.get("id")}: kind must be one of {list(KINDS)}')
        if n.get("kind") == "agent" and (not n.get("model") or "tools" not in n):
            errs.append(f'node {n.get("id")}: agent needs "model" and "tools"')
        if n.get("id") in ids:
            errs.append(f'duplicate node id {n.get("id")}')
        ids.add(n.get("id"))
        if len(n.get("job", "")) > 160:
            errs.append(f'node {n.get("id")}: "job" should be one short sentence (<= 160 characters)')
    cells = {}
    for n in spec["nodes"]:
        key = (n.get("lane"), n.get("row"))
        if key in cells:
            errs.append(f'nodes {cells[key]} and {n["id"]} share lane {key[0]}, row {key[1]}')
        cells[key] = n.get("id")
    row_of = {n["id"]: n.get("row") for n in spec["nodes"]}
    pairs = {(e["from"], e["to"]) for e in spec["edges"]}
    for a, b in pairs:
        if (b, a) in pairs and a < b and row_of.get(a) == row_of.get(b):
            errs.append(f'edges {a}->{b} and {b}->{a} join two nodes in the same row and would overlap; '
                        'state the return path in a node detail instead (e.g. "missing: stop and ask the user")')
    for e in spec["edges"]:
        if e["from"] not in ids or e["to"] not in ids:
            errs.append(f'edge {e["from"]}->{e["to"]}: unknown node')
        if e.get("kind", "flow") not in EDGES:
            errs.append(f'edge {e["from"]}->{e["to"]}: kind must be one of {list(EDGES)}')
    return errs


def layout(spec):
    lane_x = {l["id"]: i * LANE_W + LANE_W / 2 for i, l in enumerate(spec["lanes"])}
    labels, heights = {}, {}
    for n in spec["nodes"]:
        labels[n["id"]], heights[n["id"]] = node_label(n)
    rows = sorted({n["row"] for n in spec["nodes"]})
    row_h = {r: max(heights[n["id"]] for n in spec["nodes"] if n["row"] == r) for r in rows}
    total = TOP + 2 * EDGE_PAD + sum(row_h.values()) + ROW_GAP * (len(rows) - 1)
    y, row_y = total - TOP - EDGE_PAD, {}
    for r in rows:                               # Graphviz y grows upward
        row_y[r] = y - row_h[r] / 2
        y -= row_h[r] + ROW_GAP
    pos = {n["id"]: (lane_x[n["lane"]], row_y[n["row"]]) for n in spec["nodes"]}
    return pos, labels, total


def to_dot(spec, pos, labels, total):
    width = LANE_W * len(spec["lanes"])
    out = ['digraph G {', 'graph [splines=ortho, outputorder=edgesfirst, pad=0, bgcolor="transparent", '
           'fontname="Helvetica"];', 'node [fontname="Helvetica", fontsize=11];',
           'edge [fontname="Helvetica", fontsize=10, arrowsize=0.8];',
           '__c1 [shape=point, style=invis, width=0.01, pos="0,0!"];',
           f'__c2 [shape=point, style=invis, width=0.01, pos="{width:.1f},{total:.1f}!"];']
    for n in spec["nodes"]:
        x, y = pos[n["id"]]
        shape = "note" if n["kind"] == "record" else "box"
        style = "filled" if n["kind"] == "record" else "rounded,filled"
        out.append(f'"{n["id"]}" [shape={shape}, style="{style}", fillcolor="{T["node-fill"]}", '
                   f'color="{T["kind-" + n["kind"]]}", penwidth=1.5, width={NODE_W / 72:.3f}, margin="0.12,0.06", '
                   f'pos="{x:.1f},{y:.1f}!", label={labels[n["id"]]}];')
    for e in spec["edges"]:
        k = e.get("kind", "flow")
        out.append(f'"{e["from"]}" -> "{e["to"]}" [color="{T["edge-" + k]}", penwidth={EDGE_W[k]}];')
    out.append("}")
    return "\n".join(out)


# ---------- geometry ----------
def _nums(s):
    return [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", s)]


def _pad(b, p):
    return (b[0] - p, b[1] - p, b[2] + p, b[3] + p)


def _hit(a, b):
    return a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]


def _seg_hits(p, q, r):
    """True if segment p-q meets rectangle r = (x0, y0, x1, y1) (Liang-Barsky clipping)."""
    (x1, y1), (x2, y2) = p, q
    dx, dy = x2 - x1, y2 - y1
    u0, u1 = 0.0, 1.0
    for pi, qi in ((-dx, x1 - r[0]), (dx, r[2] - x1), (-dy, y1 - r[1]), (dy, r[3] - y1)):
        if pi == 0:
            if qi < 0:
                return False
        else:
            t = qi / pi
            if pi < 0:
                u0 = max(u0, t)
            else:
                u1 = min(u1, t)
            if u0 > u1:
                return False
    return True


def parse_svg(svg):
    """Node boxes, edge polylines and arrowhead boxes from the Graphviz SVG (graph coordinates)."""
    nodes, edges, heads = {}, {}, []
    for g in re.finditer(r'<g id="(node\d+|edge\d+)" class="(node|edge)">(.*?)</g>', svg, re.S):
        gid, cls, body = g.groups()
        title = html.unescape(re.search(r"<title>(.*?)</title>", body).group(1))
        if cls == "node":
            if title.startswith("__"):
                continue
            pts = []
            for p in re.findall(r'(?:points|d)="([^"]+)"', body):
                v = _nums(p)
                pts += list(zip(v[0::2], v[1::2]))
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            nodes[title] = (min(xs), min(ys), max(xs), max(ys))
        else:
            lines = []
            for p in re.findall(r'<path[^>]* d="([^"]+)"', body):
                v = _nums(p)
                lines.append(list(zip(v[0::2], v[1::2])))
            edges[title] = (gid, lines)
            for p in re.findall(r'<polygon[^>]* points="([^"]+)"', body):
                v = _nums(p)
                xs, ys = v[0::2], v[1::2]
                heads.append((title, (min(xs), min(ys), max(xs), max(ys))))
    return nodes, edges, heads


def place_labels(spec, svg, lane_titles):
    """Put each edge label beside its edge at the first spot that clears everything; insert it in the SVG."""
    nodes, edges, heads = parse_svg(svg)
    segs = [(t, a, b) for t, (_, lines) in edges.items() for pl in lines for a, b in zip(pl, pl[1:])]
    placed = []
    h = LABEL_FS * 1.05
    for e in spec["edges"]:
        if not e.get("label"):
            continue
        key = f'{e["from"]}->{e["to"]}'
        gid, lines = edges[key]
        txt = e["label"]
        w = text_w(txt, LABEL_FS, bold=True)
        own = sorted(((a, b) for pl in lines for a, b in zip(pl, pl[1:])),
                     key=lambda s: -abs(s[1][0] - s[0][0]) - abs(s[1][1] - s[0][1]))
        best, best_bad = None, None
        for a, b in own:
            vertical = abs(b[0] - a[0]) < abs(b[1] - a[1])
            for t in (0.5, 0.35, 0.65, 0.2, 0.8):
                px, py = a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t
                if vertical:
                    cands = [(px + LABEL_GAP, py - h / 2, px + LABEL_GAP + w, py + h / 2),
                             (px - LABEL_GAP - w, py - h / 2, px - LABEL_GAP, py + h / 2)]
                else:
                    cands = [(px - w / 2, py - LABEL_GAP - h, px + w / 2, py - LABEL_GAP),
                             (px - w / 2, py + LABEL_GAP, px + w / 2, py + LABEL_GAP + h)]
                for box in cands:
                    pb = _pad(box, CLEAR)
                    bad = sum(_hit(pb, nb) for nb in nodes.values())
                    bad += sum(_seg_hits(p, q, _pad(box, CLEAR if tt != key else 1.0)) for tt, p, q in segs)
                    bad += sum(_hit(pb, hb) for _, hb in heads)
                    bad += sum(_hit(pb, lb) for lb in placed)
                    bad += sum(_hit(pb, tb) for tb in lane_titles)
                    if best_bad is None or bad < best_bad:
                        best, best_bad = box, bad
                    if bad == 0:
                        break
                if best_bad == 0:
                    break
            if best_bad == 0:
                break
        placed.append(best)
        col = T["edge-" + e.get("kind", "flow")]
        tag = (f'<text class="elabel" x="{best[0]:.1f}" y="{best[3] - 0.22 * LABEL_FS:.1f}" '
               f'font-size="{LABEL_FS:g}" font-weight="600" text-anchor="start" fill="{col}">{esc(txt)}</text>')
        svg = re.sub(rf'(<g id="{gid}" class="edge">.*?)(</g>)', lambda m: m.group(1) + tag + "\n" + m.group(2),
                     svg, count=1, flags=re.S)
    return svg


def add_lanes(svg, spec, total):
    """Draw lane backgrounds and titles under the graph (same coordinates as the nodes)."""
    from collections import Counter
    kinds = {}
    for n in spec["nodes"]:
        kinds.setdefault(n["lane"], Counter())[n["kind"]] += 1
    accent = {lane: T["kind-" + c.most_common(1)[0][0]] for lane, c in kinds.items()}
    parts, titles = [], []
    for i, l in enumerate(spec["lanes"]):
        x0, w = i * LANE_W + 6, LANE_W - 12
        col = accent.get(l["id"], T["text-muted"])
        parts.append(f'<rect x="{x0:.1f}" y="{-total:.1f}" width="{w:.1f}" height="{total:.1f}" rx="10" '
                     f'fill="{T["lane-bg"]}" stroke="{T["lane-border"]}"/>')
        parts.append(f'<rect x="{x0 + 14:.1f}" y="{-total + 12:.1f}" width="{w - 28:.1f}" height="2" rx="1" '
                     f'fill="{col}" opacity="0.85"/>')
        label = l["label"].upper()
        tw = text_w(label, 12) * 1.25
        cx, ty = x0 + w / 2, -total + 34
        parts.append(f'<text class="lane-title" x="{cx:.1f}" y="{ty:.1f}" text-anchor="middle" '
                     f'font-size="12" fill="{col}">{esc(label)}</text>')
        titles.append((cx - tw / 2, ty - 12, cx + tw / 2, ty + 4))
    m = re.search(r'(<g id="graph0"[^>]*>\s*<title>.*?</title>\s*(<polygon[^>]*/>)?)', svg, re.S)
    return svg[:m.end()] + "\n" + "\n".join(parts) + svg[m.end():], titles


# ---------- lint ----------
def lint(svg, spec, lane_titles):
    """Label/node, label/edge, label/label, label/lane-title overlaps and edges through nodes."""
    nodes, edges, _ = parse_svg(svg)
    labels = []
    for title, (gid, _) in edges.items():
        body = re.search(rf'<g id="{gid}" class="edge">(.*?)</g>', svg, re.S).group(1)
        for t in re.finditer(r'<text[^>]*x="([-\d.]+)" y="([-\d.]+)"[^>]*font-size="([\d.]+)"[^>]*>(.*?)</text>', body):
            x, y, fs, txt = float(t.group(1)), float(t.group(2)), float(t.group(3)), html.unescape(t.group(4))
            w = text_w(txt, fs, bold='font-weight="600"' in t.group(0))
            labels.append((title, txt, (x, y - 0.83 * fs, x + w, y + 0.22 * fs)))
    issues = []
    for e, txt, bb in labels:
        for n, nb in nodes.items():
            if _hit(_pad(bb, 2), nb):
                issues.append({"type": "label_over_node", "edge": e, "label": txt, "node": n})
        for e2, (_, lines) in edges.items():
            pad = 0.5 if e2 == e else 1.5
            if any(_seg_hits(p, q, _pad(bb, pad)) for pl in lines for p, q in zip(pl, pl[1:])):
                issues.append({"type": "label_over_edge", "edge": e, "label": txt, "other_edge": e2})
        for tb in lane_titles:
            if _hit(_pad(bb, 2), tb):
                issues.append({"type": "label_over_lane_title", "edge": e, "label": txt})
    for i, (e, txt, bb) in enumerate(labels):
        for (e2, txt2, bb2) in labels[i + 1:]:
            if _hit(_pad(bb, 2), bb2):
                issues.append({"type": "label_over_label", "labels": [txt, txt2]})
    for e, (_, lines) in edges.items():
        a, b = e.split("->")
        for n, nb in nodes.items():
            if n in (a, b):
                continue
            if any(_seg_hits(p, q, _pad(nb, -3)) for pl in lines for p, q in zip(pl, pl[1:])):
                issues.append({"type": "edge_through_node", "edge": e, "node": n})
    return {"nodes": len(nodes), "edge_labels": len(labels), "issues": issues}


def page(spec, svg, spec_path):
    nodes = spec["nodes"]
    legend = "".join(f'<span class="chip"><span class="sw" style="border-color:{T["kind-" + k]}"></span>{KIND_TEXT[k]}</span>'
                     for k in KINDS if any(n["kind"] == k for n in nodes))
    legend += "".join(f'<span class="chip"><span class="ln" style="background:{T["edge-" + k]}"></span>{EDGE_TEXT[k]}</span>'
                      for k in EDGES if any(e.get("kind", "flow") == k for e in spec["edges"]))
    rows = "".join(
        f'<tr><td class="stage">{esc(n["title"])}</td>'
        f'<td class="kind" style="color:{T["kind-" + n["kind"]]}">{esc(n["kind"])}</td>'
        f'<td class="mono">{esc(n.get("model", "—"))}</td>'
        f'<td class="mono">{esc(", ".join(n.get("tools", [])) or "—")}</td><td>{esc(n["job"])}</td></tr>'
        for n in nodes if n["kind"] in ("agent", "script"))
    summary = "".join(f"<p>{esc(p)}</p>" for p in (spec["summary"] if isinstance(spec["summary"], list)
                                                    else [spec["summary"]]))
    try:
        src = spec_path.resolve().relative_to(ROOT)
    except ValueError:
        src = spec_path
    page_w = int(LANE_W * len(spec["lanes"]) * 4 / 3) + 104
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(spec["title"])}</title>
<style>
{THEME_TEXT}
:root {{ --page-width: {page_w}px; }}
</style></head><body><main class="page">
<div class="eyebrow">Workflow design</div>
<h1>{esc(spec["title"])}</h1>
{summary}
<div class="badge">{esc(spec.get("status", ""))}</div>
<div class="card">
{svg}
</div>
<h2>Legend</h2>
<div class="legend">{legend}</div>
<h2>Agents and scripts</h2>
<table><tr><th>Stage</th><th>Kind</th><th>Model</th><th>Tools</th><th>Job</th></tr>{rows}</table>
<footer>Rendered from {esc(src)} by scripts/flowchart/render_flowchart.py with scripts/flowchart/theme.css</footer>
</main></body></html>
"""


def render(spec_path, png=False):
    spec = json.loads(spec_path.read_text())
    errs = validate(spec)
    if errs:
        return 2, {"spec": str(spec_path), "spec_errors": errs}
    pos, labels, total = layout(spec)
    dot = to_dot(spec, pos, labels, total)
    r = subprocess.run(["neato", "-n2", "-Tsvg"], input=dot, capture_output=True, text=True)
    if r.returncode != 0:
        return 3, {"spec": str(spec_path), "graphviz_error": r.stderr}
    svg = re.sub(r"^.*?(<svg)", r"\1", r.stdout, flags=re.S)
    svg, lane_titles = add_lanes(svg, spec, total)
    svg = place_labels(spec, svg, lane_titles)
    out = spec_path.with_suffix(".html")
    out.write_text(page(spec, svg, spec_path))
    report = {"spec": str(spec_path), **lint(svg, spec, lane_titles)}
    warnings = [l for l in r.stderr.splitlines() if l.strip()]
    if warnings:
        report["graphviz_warnings"] = warnings
    if png:
        import tempfile
        shot = Path(tempfile.gettempdir()) / f"flowchart_{spec_path.parent.name}.png"
        chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        h = int(total * 4 / 3) + 1100
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                        f"--window-size={int(LANE_W * len(spec['lanes']) * 4 / 3) + 180},{h}",
                        f"--screenshot={shot}", out.resolve().as_uri()], capture_output=True)
        report["png"] = str(shot)
    report["html"] = str(out)
    return (1 if report["issues"] else 0), report


def main(argv):
    png = "--png" in argv
    if "--all" in argv:
        specs = sorted((ROOT / "workflows").glob("*/flowchart.json"))
    else:
        specs = [Path(a) for a in argv[1:] if not a.startswith("--")]
    status = 0
    for s in specs:
        code, report = render(s, png)
        print(json.dumps(report, indent=1))
        status = max(status, code)
    return status


if __name__ == "__main__":
    sys.exit(main(sys.argv))
