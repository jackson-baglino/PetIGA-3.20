#!/usr/bin/env python3
"""Render a workflow flowchart (HTML with inline SVG) from a JSON spec, then lint it.

    python scripts/flowchart/render_flowchart.py workflows/<name>/flowchart.json [--png]

Layout is deterministic: each node sits on a grid (its lane is the column, its
"row" the row); Graphviz (neato -n2, orthogonal splines) routes the edges and
places edge labels away from nodes. The agent writing the spec never places
coordinates. After rendering, the script lints the SVG for label/node,
label/edge and label/label overlaps and prints a JSON report; exit status 1
if anything overlaps. --png also writes a screenshot (headless Chrome) for
visual review.

Spec (see scripts/flowchart/README.md for the full schema):
  {"title", "summary", "lanes": [{"id","label"}],
   "nodes": [{"id","lane","row","kind","title","job","model","tools","details"}],
   "edges": [{"from","to","label","kind"}]}
kind (node): user | script | agent | record.  kind (edge): flow | stop | loop | data.
"""
import html
import json
import math
import re
import subprocess
import sys
import textwrap
from pathlib import Path

NODE_W = 290.0          # node width [pt]
LANE_W = 330.0          # lane width [pt]
ROW_GAP = 64.0          # vertical gap between rows [pt]
EDGE_PAD = 18.0         # gap above the first and below the last row [pt]
TOP = 46.0              # space for lane titles [pt]
WRAP = 46               # characters per line inside a node
LINE_H = 14.0           # line height estimate [pt]

STYLE = {  # node kind -> (border, lane tint)
    "user":   ("#c98a2b", "#fdf6ec"),
    "script": ("#2b6cb0", "#eef4fb"),
    "agent":  ("#2f855a", "#f1f7ee"),
    "record": ("#7a7a7a", "#f5f5f5"),
}
EDGE = {  # edge kind -> (colour, width, style)
    "flow": ("#4a4a4a", 1.4, "solid"),
    "data": ("#8a8a8a", 1.0, "solid"),
    "stop": ("#b3261e", 1.4, "solid"),
    "loop": ("#b07a00", 1.3, "solid"),
}
KIND_TEXT = {"user": "user decision point", "script": "deterministic script: same input, same output",
             "agent": "agent (model and tools listed in the box)", "record": "record written by a stage"}
EDGE_TEXT = {"flow": "main flow", "data": "writes a record", "stop": "stop and ask the user",
             "loop": "loop back for a re-check"}


def esc(s):
    return html.escape(str(s), quote=False)


def wrap(s, width=WRAP):
    return textwrap.wrap(str(s), width) or [""]


def node_label(n):
    """Graphviz HTML-like label: title, model/tools line, one-sentence job, then details."""
    rows = [f'<TR><TD ALIGN="LEFT"><B>{esc(n["title"])}</B></TD></TR>']
    lines = 1
    if n["kind"] == "agent":
        meta = f'{n["model"]} · tools: {", ".join(n["tools"]) or "none"}'
        for l in wrap(meta):
            rows.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="10" COLOR="#2f5d3a">{esc(l)}</FONT></TD></TR>')
            lines += 1
    job = n.get("job", "")
    for l in wrap(job):
        rows.append(f'<TR><TD ALIGN="LEFT"><I>{esc(l)}</I></TD></TR>')
        lines += 1
    for d in n.get("details", []):
        for i, l in enumerate(wrap(d, WRAP - 2)):
            rows.append(f'<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="10" COLOR="#444444">{"• " if i == 0 else "  "}{esc(l)}</FONT></TD></TR>')
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
        if n.get("kind") not in STYLE:
            errs.append(f'node {n.get("id")}: kind must be one of {sorted(STYLE)}')
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
        if e.get("kind", "flow") not in EDGE:
            errs.append(f'edge {e["from"]}->{e["to"]}: kind must be one of {sorted(EDGE)}')
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
    out = ['digraph G {', 'graph [splines=ortho, forcelabels=true, outputorder=edgesfirst, pad=0, '
           'fontname="Helvetica"];', f'node [fontname="Helvetica", fontsize=11];',
           'edge [fontname="Helvetica", fontsize=10, arrowsize=0.8];',
           f'__c1 [shape=point, style=invis, width=0.01, pos="0,0!"];',
           f'__c2 [shape=point, style=invis, width=0.01, pos="{width:.1f},{total:.1f}!"];']
    for n in spec["nodes"]:
        x, y = pos[n["id"]]
        border = STYLE[n["kind"]][0]
        shape = "note" if n["kind"] == "record" else "box"
        style = "filled" if n["kind"] == "record" else "rounded,filled"
        out.append(f'"{n["id"]}" [shape={shape}, style="{style}", fillcolor="white", color="{border}", '
                   f'penwidth=1.5, width={NODE_W / 72:.3f}, margin="0.12,0.06", '
                   f'pos="{x:.1f},{y:.1f}!", label={labels[n["id"]]}];')
    for e in spec["edges"]:
        col, w, st = EDGE[e.get("kind", "flow")]
        lab = f', xlabel="{esc(e["label"])}", fontcolor="{col}"' if e.get("label") else ""
        out.append(f'"{e["from"]}" -> "{e["to"]}" [color="{col}", penwidth={w}, style={st}{lab}];')
    out.append("}")
    return "\n".join(out)


def add_lanes(svg, spec, total):
    """Draw lane backgrounds and titles under the graph (same coordinates as the nodes)."""
    from collections import Counter
    kinds = {}
    for n in spec["nodes"]:
        kinds.setdefault(n["lane"], Counter())[n["kind"]] += 1
    tint = {lane: STYLE[c.most_common(1)[0][0]][1] for lane, c in kinds.items()}  # colour of the lane's main kind
    parts = []
    for i, l in enumerate(spec["lanes"]):
        x0 = i * LANE_W + 6
        parts.append(f'<rect x="{x0:.1f}" y="{-total:.1f}" width="{LANE_W - 12:.1f}" height="{total:.1f}" '
                     f'fill="{tint.get(l["id"], "#f7f7f7")}" stroke="#dddddd"/>')
        parts.append(f'<text x="{x0 + (LANE_W - 12) / 2:.1f}" y="{-total + 26:.1f}" text-anchor="middle" '
                     f'font-family="Helvetica" font-size="13" font-weight="bold" fill="#444">{esc(l["label"])}</text>')
    # insert right after the graph group's opening tag and its background polygon
    m = re.search(r'(<g id="graph0"[^>]*>\s*<title>.*?</title>\s*(<polygon[^>]*/>)?)', svg, re.S)
    return svg[:m.end()] + "\n" + "\n".join(parts) + svg[m.end():]


# ---------- lint ----------
def _nums(s):
    return [float(v) for v in re.findall(r"-?\d+(?:\.\d+)?", s)]


def lint(svg):
    """Label/node, label/edge and label/label overlaps, from the rendered SVG."""
    boxes, labels, paths = [], [], []
    for g in re.finditer(r'<g id="(node\d+|edge\d+)" class="(node|edge)">(.*?)</g>', svg, re.S):
        gid, cls, body = g.groups()
        title = re.search(r"<title>(.*?)</title>", body).group(1)
        if cls == "node":
            if title.startswith("__"):
                continue
            pts = []
            for p in re.findall(r'(?:points|d)="([^"]+)"', body):
                v = _nums(p)
                pts += list(zip(v[0::2], v[1::2]))
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            boxes.append((html.unescape(title), min(xs), min(ys), max(xs), max(ys)))
        else:
            for p in re.findall(r'<path[^>]* d="([^"]+)"', body):
                v = _nums(p)
                paths.append((html.unescape(title), list(zip(v[0::2], v[1::2]))))
            for t in re.finditer(r'<text[^>]*x="([-\d.]+)" y="([-\d.]+)"[^>]*font-size="([\d.]+)"[^>]*>(.*?)</text>', body):
                x, y, fs, txt = float(t.group(1)), float(t.group(2)), float(t.group(3)), html.unescape(t.group(4))
                w = 0.52 * fs * len(txt)
                anchor = re.search(r'text-anchor="(\w+)"', t.group(0))
                x0 = x - w / 2 if anchor and anchor.group(1) == "middle" else (x - w if anchor and anchor.group(1) == "end" else x)
                labels.append((html.unescape(title), txt, x0, y - 0.8 * fs, x0 + w, y + 0.25 * fs))
    issues = []
    hit = lambda a, b: a[0] < b[2] and b[0] < a[2] and a[1] < b[3] and b[1] < a[3]
    for (e, txt, *bb) in labels:
        for (n, *nb) in boxes:
            if hit(bb, nb):
                issues.append({"type": "label_over_node", "edge": e, "label": txt, "node": n})
        for (e2, pts) in paths:
            if e2 == e:
                continue
            for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
                for k in range(11):
                    px, py = x1 + (x2 - x1) * k / 10, y1 + (y2 - y1) * k / 10
                    if bb[0] < px < bb[2] and bb[1] < py < bb[3]:
                        issues.append({"type": "label_over_edge", "edge": e, "label": txt, "other_edge": e2})
                        break
                else:
                    continue
                break
    for i, (e, txt, *bb) in enumerate(labels):
        for (e2, txt2, *bb2) in labels[i + 1:]:
            if hit(bb, bb2):
                issues.append({"type": "label_over_label", "labels": [txt, txt2]})
    return {"nodes": len(boxes), "edge_labels": len(labels), "issues": issues}


def page(spec, svg):
    nodes = {n["id"]: n for n in spec["nodes"]}
    legend_nodes = "".join(f'<tr><td><span class="sw" style="border-color:{STYLE[k][0]}"></span></td><td>{KIND_TEXT[k]}</td></tr>'
                           for k in STYLE if any(n["kind"] == k for n in spec["nodes"]))
    legend_edges = "".join(f'<tr><td><span class="ln" style="background:{EDGE[k][0]}"></span></td><td>{EDGE_TEXT[k]}</td></tr>'
                           for k in EDGE if any(e.get("kind", "flow") == k for e in spec["edges"]))
    rows = "".join(
        f'<tr><td>{esc(n["title"])}</td><td>{esc(n["kind"])}</td><td>{esc(n.get("model", "—"))}</td>'
        f'<td>{esc(", ".join(n.get("tools", [])) or "—")}</td><td>{esc(n["job"])}</td></tr>'
        for n in nodes.values() if n["kind"] in ("agent", "script"))
    summary = "".join(f"<p>{esc(p)}</p>" for p in (spec["summary"] if isinstance(spec["summary"], list) else [spec["summary"]]))
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><title>{esc(spec["title"])}</title>
<style>
body {{ font-family: -apple-system, "Helvetica Neue", Arial, sans-serif; margin: 24px; color: #1d1d1f; max-width: {LANE_W * len(spec["lanes"]) + 60:.0f}px; }}
h1 {{ font-size: 20px; margin: 0 0 8px; }} h2 {{ font-size: 15px; margin: 22px 0 6px; }}
p {{ font-size: 14px; line-height: 1.45; margin: 0 0 8px; max-width: 900px; }}
table {{ border-collapse: collapse; font-size: 12.5px; }} td, th {{ border: 1px solid #ddd; padding: 4px 8px; text-align: left; vertical-align: top; }}
th {{ background: #f4f4f4; }} .sw {{ display: inline-block; width: 22px; height: 12px; border: 2px solid; border-radius: 3px; }}
.ln {{ display: inline-block; width: 26px; height: 3px; vertical-align: middle; }}
.meta {{ color: #666; font-size: 12px; }}
</style></head><body>
<h1>{esc(spec["title"])}</h1>
{summary}
<p class="meta">{esc(spec.get("status", ""))}</p>
{svg}
<h2>Legend</h2>
<table>{legend_nodes}{legend_edges}</table>
<h2>Agents and scripts</h2>
<table><tr><th>Stage</th><th>Kind</th><th>Model</th><th>Tools</th><th>Job</th></tr>{rows}</table>
</body></html>
"""


def main(argv):
    spec_path = Path(argv[1])
    spec = json.loads(spec_path.read_text())
    errs = validate(spec)
    if errs:
        print(json.dumps({"spec_errors": errs}, indent=1))
        return 2
    pos, labels, total = layout(spec)
    dot = to_dot(spec, pos, labels, total)
    r = subprocess.run(["neato", "-n2", "-Tsvg"], input=dot, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr)
        return 3
    svg = re.sub(r"^.*?(<svg)", r"\1", r.stdout, flags=re.S)
    svg = add_lanes(svg, spec, total)
    out = spec_path.with_suffix(".html")
    out.write_text(page(spec, svg))
    report = lint(svg)
    warnings = [l for l in r.stderr.splitlines() if l.strip()]
    if warnings:
        report["graphviz_warnings"] = warnings
    if "--png" in argv:
        png = spec_path.with_suffix(".png")
        chrome = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
        h = int(total * 1.34) + 900
        subprocess.run([chrome, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                        f"--window-size={int(LANE_W * len(spec['lanes']) * 1.34) + 120},{h}",
                        f"--screenshot={png}", out.resolve().as_uri()], capture_output=True)
        report["png"] = str(png)
    report["html"] = str(out)
    print(json.dumps(report, indent=1))
    return 1 if report["issues"] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
