#!/usr/bin/env python3
"""Run one review wave from an approved manifest: extract -> verify -> cross-tier claim check -> replicate sample -> report.

  python -m agentic.lit.run_wave --manifest research/Review/waves/wave1.json --dry-run          # cost estimate only
  python -m agentic.lit.run_wave --manifest research/Review/waves/wave1.json --workers 4        # the real run
  python -m agentic.lit.run_wave --manifest ... --only slug1,slug2                               # subset / resume

Stages (course Lecture 5: a loop that acts needs an external check; Lecture 4: vote only over independent samples):
  1. extract   one fresh `paper-extractor` (Sonnet) per paper, full text in its context, schema-enforced JSON
  2. verify    deterministic: every quote string-matched on its page, coverage, DOI, image-check queue (verify_record)
  3. check     `claim-checker` on a DIFFERENT tier (Opus): all mandatory claims + 25 % of the rest, isolated context
  4. replicate a second independent Sonnet extraction of 20 % of the papers; agreement and dissent are logged (consistency)
  5. report    wave_report.md: metrics, dissent, image-check queue, token ledger. Stops and says so if a threshold fails.
Pauses for a human: after the manifest (approve the reading list) and after the report (read the flagged items).
"""
from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from agentic.lit import common
from agentic.lit.consistency import compare
from agentic.lit import run_extract

# Gates: a wave is reported as FAILED (not silently accepted) if any of these is violated.
MIN_QUOTE_MATCH = 0.90
MIN_COVERAGE = 1.0
MIN_CHECKER_SUPPORTED = 0.90
EST_USD_PER_PAGE_SONNET = 0.016     # measured in the Ji and Karma pilot (list prices); replace after the first wave
EST_USD_CHECK_PER_PAPER = 0.30


def pages_of(pdf: str) -> int:
    out = subprocess.run(["mutool", "info", pdf], capture_output=True, text=True).stdout
    for tok in out.split("\n"):
        if tok.strip().startswith("Pages:"):
            return int(tok.split(":")[1])
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--only", default="")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--replicate-fraction", type=float, default=0.20)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    papers = json.loads(a.manifest.read_text())
    if a.only:
        keep = set(a.only.split(","))
        papers = [p for p in papers if p["slug"] in keep]
    for p in papers:
        p["n_pages"] = pages_of(p["pdf"])
    est = sum(p["n_pages"] * EST_USD_PER_PAGE_SONNET + EST_USD_CHECK_PER_PAPER for p in papers) * (1 + a.replicate_fraction * 0.6)
    print(f"{len(papers)} papers, {sum(p['n_pages'] for p in papers)} pages; estimated list-price cost ${est:.0f} "
          f"(extract + Opus check + {a.replicate_fraction:.0%} replicate). Nothing is billed to a personal account.")
    if a.dry_run:
        return

    specs = [f"{p['id']}={p['slug']}:{','.join(p['topics'])}" for p in papers]
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(lambda s: run_extract.run_one(s, "sonnet", "", True), specs))

    rd = common.REVIEW / "records"
    sub = [p["slug"] for p in papers]
    for slug in sub:                                   # stage 3 (sequential per record; the checker itself is parallel inside)
        subprocess.run([sys.executable, "-m", "agentic.lit.check_claims", "--record", slug, "--model", "opus", "--workers", "4"], cwd=common.ROOT, check=False)

    rng = random.Random(a.seed)
    rep = rng.sample(papers, max(1, int(round(a.replicate_fraction * len(papers)))))
    rep_specs = [f"{p['id']}={p['slug']}:{','.join(p['topics'])}" for p in rep]
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(lambda s: run_extract.run_one(s, "sonnet", "_rep2", False), rep_specs))

    L = [f"# Wave report: {a.manifest.stem}", "", f"{len(papers)} papers.", "",
         "| record | claims | quote match | exact page | coverage | checker supported | image queue | status |", "|---|---|---|---|---|---|---|---|"]
    failed = []
    for slug in sub:
        v = json.loads((rd / f"{slug}.verification.json").read_text())
        m = v["metrics"]
        ck = v.get("checker", {}).get("supported_rate")
        bad = (m["quote_match_rate"] or 0) < MIN_QUOTE_MATCH or m["coverage"] < MIN_COVERAGE or (ck is not None and ck < MIN_CHECKER_SUPPORTED) or v.get("problems")
        if bad:
            failed.append(slug)
        L.append(f"| {slug} | {m['n_claims']} | {m['quote_match_rate']} | {m['exact_page_rate']} | {m['coverage']} | {ck} | {m['image_queue_items']} | {'REVIEW' if bad else 'ok'} |")
    L += ["", f"**{len(failed)} records need your review before use:** {', '.join(failed) if failed else 'none'}", "", "## Replicate agreement (independent second Sonnet extraction)"]
    for p in rep:
        a1 = json.loads((rd / f"{p['slug']}.json").read_text())
        a2 = json.loads((rd / f"{p['slug']}_rep2.json").read_text())
        c = compare(a1, a2)
        L.append(f"- {p['slug']}: probe agreement {c['probe_status_agreement']} (n={c['n_probes_compared']}), number overlap {c['number_value_jaccard']}, "
                 f"equation overlap {c['equation_number_jaccard']}; dissent: " + ("; ".join(f"{d['qid']} {d['a']} vs {d['b']}" for d in c["probe_dissent"]) or "none"))
    tot = {}
    for line in open(common.LEDGER):
        d = json.loads(line)
        t = tot.setdefault(d["agent"] + "/" + d["model"], [0, 0, 0, 0.0])
        t[0] += 1; t[1] += d["input"] + d["cache_read"] + d["cache_write"]; t[2] += d["output"]; t[3] += d["cost_usd_list"] or 0
    L += ["", "## Token ledger (list prices)", "| agent/model | calls | input tokens | output tokens | $ |", "|---|---|---|---|---|"]
    L += [f"| {k} | {v[0]} | {v[1]:,} | {v[2]:,} | {v[3]:.2f} |" for k, v in tot.items()]
    (common.REVIEW / "memos").mkdir(parents=True, exist_ok=True)
    (common.REVIEW / "memos" / f"{a.manifest.stem}_report.md").write_text("\n".join(L) + "\n")
    print("report -> memos/" + f"{a.manifest.stem}_report.md;", "FAILED gates:" if failed else "all gates passed", ", ".join(failed))


if __name__ == "__main__":
    main()
