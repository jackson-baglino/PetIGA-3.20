#!/usr/bin/env python3
"""M3 evaluation view: the agent's derivation next to doc 18 (the 2026-10-18 F2 derivation), for the user.

Doc 18 is hidden from the agents, so this comparison measures the derivation independently (WORKFLOW.md §1,
baseline). It is a view, not a check: notations differ (doc 18 writes w for the indicator and psi for the ice
fraction; the line writes psi for the indicator and phi for the ice field), so the user compares term by term.
Writes <run>/gates/doc18_comparison.md.

    venv_workflows/bin/python scripts/model_development/compare_doc18.py --run <run folder>
"""
import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC18 = ROOT / "docs" / "planning" / "18_f2_implementation_plan.md"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    a = ap.parse_args(argv)
    rdir = Path(a.run)
    led = json.loads((rdir / "derivation" / "ledger.json").read_text())
    text = DOC18.read_text()
    model = re.search(r"## 1\. The model\n(.*?)\n## 2\.", text, re.S)
    rev7 = re.search(r"## 7\..*?\n(.*?)\n## 8\.", text, re.S)
    rev8 = re.search(r"## 8\..*?\n(.*)$", text, re.S)
    out = ["# Derivation vs doc 18 (for the user; not a check)", "",
           "Notation: doc 18 uses w = pore indicator and psi = ice fraction of the pore; this run uses psi = "
           "grain indicator and phi = ice phase field. Compare term by term: storage, gradient, well, source, "
           "wall term, vapor row, heat row, floors.", "", "## Doc 18, section 1 (the model)", "",
           model.group(1).strip() if model else "(not found)", "", "## Doc 18, sections 7 and 8 (revisions)", "",
           (rev7.group(1).strip() if rev7 else "") + "\n\n" + (rev8.group(1).strip() if rev8 else ""), "",
           "## This run: derived steps (in ledger order)", ""]
    for s in led["steps"]:
        if s.get("type") in ("derived", "choice", "assumption"):
            out.append(f"- **{s['id']}** ({s['type']}{', ' + s.get('operation', '') if s.get('operation') else ''})"
                       f": ${s['statement']}$")
    (rdir / "gates" / "doc18_comparison.md").write_text("\n".join(out) + "\n")
    print(f"wrote {rdir / 'gates' / 'doc18_comparison.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
