# Workflow flowcharts

`render_flowchart.py` turns `workflows/<name>/flowchart.json` into
`flowchart.html` (and, with `--png`, a screenshot in the system temporary folder for review; screenshots are not saved in the repo). Layout and
routing are done by the script and Graphviz, not by hand, so every chart has
the same structure and arrows always end on their target.

```bash
venv_lunar/bin/python scripts/flowchart/render_flowchart.py workflows/<name>/flowchart.json --png
venv_lunar/bin/python scripts/flowchart/render_flowchart.py --all          # every workflows/*/flowchart.json
venv_lunar/bin/python -m unittest scripts/flowchart/test_render_flowchart.py
```

**Theme.** `theme.css` is the one place that sets colors, type and spacing for
every workflow page: the openclaw.ai layout conventions on the Monokai Pro
palette. The renderer inlines it into each page (pages stay self-contained)
and takes the chart colors from its `--kind-*`, `--edge-*`, `--node-*`,
`--lane-*` and `--text-*` values. To change the look of all pages, edit
`theme.css` and run `--all`.

Exit 0: clean. Exit 1: the lint found overlapping text (listed in the JSON
report). Exit 2: the spec is invalid (errors listed). Requires Graphviz
(`neato`) and, for `--png`, Google Chrome.

## Spec

```json
{
  "title": "workflow name",
  "summary": ["2-4 sentences: what the workflow does and how it checks itself."],
  "status": "Design for review, <date>.",
  "lanes": [{"id": "user", "label": "User"}, {"id": "script", "label": "Scripts"},
            {"id": "agent", "label": "Agents"}, {"id": "record", "label": "Records"}],
  "nodes": [
    {"id": "S2", "lane": "agent", "row": 3, "kind": "agent", "title": "S2 morty-unit-reader",
     "model": "claude-sonnet-5-5", "tools": ["Read"],
     "job": "One sentence: what this agent does.",
     "details": ["Further points, one short line each"]}
  ],
  "edges": [{"from": "S2", "to": "S3", "kind": "flow", "label": "optional, short"}]
}
```

- `kind` (node): `user`, `script`, `agent`, `record`. Agents must have `model`
  (full ID) and `tools`; every node has a one-sentence `job`.
- `kind` (edge): `flow` (main path), `data` (writes a record), `stop` (stop and
  ask the user), `loop` (re-check).
- `row` orders the stages top to bottom; nodes in one row line up. One node per
  lane per row.

## Rules the renderer enforces

- Missing fields, unknown lanes or kinds, two nodes in one cell, a `job`
  longer than 160 characters.
- No pair of opposite edges between two nodes in the same row (they would draw
  on top of each other). Put the return path in a node detail instead.
- Edge labels are placed by the renderer, not by Graphviz: each is measured
  with the Helvetica metrics and put at the first spot beside its edge that
  clears every node, edge, arrowhead, lane title and earlier label.
- Lint after rendering: no edge label over a node, over any edge (its own
  included), over another label or over a lane title, and no edge through a
  node that is not one of its ends. Segments are tested exactly, not sampled.

## Conventions

- Box styling encodes only the node kind (user, script, agent, record);
  run-once or optional stages say so in their title, not by line style.
- The page shows the title, the summary, the chart, a legend and a table of
  every agent and script with its model, tools and job.
