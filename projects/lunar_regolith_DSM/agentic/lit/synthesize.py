#!/usr/bin/env python3
"""Gate memo: one `topic-synthesizer` (Opus) call over VERIFIED RECORDS ONLY (never the PDFs). Each record is given as its markdown
(claims with page and quote, image-checked equations) followed by the independent checker's flags, so the synthesizer can discount
claims that said more than their quote. Output: $LIT_REVIEW_DIR/memos/<name>.md

  python -m agentic.lit.synthesize --slugs-file .agent_state/wave2_slugs.txt --brief brief.md --name gate_R8_diffuse_domain
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from agentic.lit import common


def record_block(slug: str) -> str:
    rec = common.REVIEW / "records"
    md = (rec / f"{slug}.md").read_text()
    v = json.loads((rec / f"{slug}.verification.json").read_text())
    ck = v.get("checker") or {}
    flags = [f"- {cid}: {d['verdict']}. {d.get('problem', '')}" for cid, d in (ck.get("verdicts") or {}).items() if d["verdict"] != "supported"]
    head = f"\n\n######## RECORD {slug} ########\n"
    tail = ("\nINDEPENDENT CHECKER (Opus, saw only claim + quote + page): supported rate "
            f"{ck.get('supported_rate')}; claims it flagged (use the quote's own words for these, not the claim's):\n" + ("\n".join(flags) if flags else "- none"))
    return head + md + tail


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slugs-file", required=True)
    ap.add_argument("--brief", required=True, help="markdown file: the gate's question and what the project needs from the memo")
    ap.add_argument("--name", required=True)
    ap.add_argument("--model", default="opus")
    a = ap.parse_args()
    slugs = [s.strip() for s in Path(a.slugs_file).read_text().split() if s.strip()]
    body = Path(a.brief).read_text() + "\n\nVERIFIED RECORDS FOLLOW.\n" + "".join(record_block(s) for s in slugs)
    memo, usage = common.call_claude("topic-synthesizer", body, model=a.model, label=f"synth:{a.name}", timeout=3600)
    out = common.REVIEW / "memos" / f"{a.name}.md"
    out.write_text(memo if isinstance(memo, str) else json.dumps(memo, indent=1))
    print(f"{len(slugs)} records, {len(body)} chars in -> {out}; ${usage['cost_usd_list']}")


if __name__ == "__main__":
    main()
