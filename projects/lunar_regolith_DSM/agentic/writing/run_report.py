#!/usr/bin/env python3
"""Driver for the report-writing workflow (a short course report).  Steps:
  1 facts (evidence.py, no model)  2 outline (Sonnet)  3 draft with {F:id} placeholders (Sonnet)  4 deterministic checks (no model)
  5 two independent critics (Opus): clarity/voice and rubric  6 revise (Sonnet) up to --rounds times  7 render + checklist of what a person must verify.
Everything is saved under --out-dir; every model call is logged by agentic.lit.common.call_claude in .agent_state/ledger.jsonl.

  python -m agentic.writing.run_report --brief B --guideline G --facts F --profile P --samples S --out-dir D [--budget 1000] [--rounds 2]
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from agentic.lit.common import call_claude
from agentic.writing import checks, render

STR = {"type": "string"}; STRS = {"type": "array", "items": STR}
OUTLINE = {"type": "object", "required": ["reasoning", "sections", "questions_for_author"], "properties": {
    "reasoning": STR, "questions_for_author": STRS,
    "sections": {"type": "array", "items": {"type": "object", "required": ["title", "purpose", "key_points", "fact_ids", "target_words"],
                 "properties": {"title": STR, "purpose": STR, "key_points": STRS, "fact_ids": STRS, "target_words": {"type": "integer"}}}}}}
CRITIC = {"type": "object", "required": ["reasoning", "issues", "reads_like_ai", "strongest_passage"], "properties": {
    "reasoning": STR, "strongest_passage": STR,
    "issues": {"type": "array", "items": {"type": "object", "required": ["id", "severity", "quote", "problem", "fix"], "properties": {
        "id": STR, "severity": {"enum": ["high", "medium", "low"]}, "quote": STR, "problem": STR, "fix": STR}}},
    "reads_like_ai": {"type": "array", "items": {"type": "object", "required": ["quote", "why"], "properties": {"quote": STR, "why": STR}}}}}
AUDIT = {"type": "object", "required": ["reasoning", "rubric", "required_elements", "claims_without_support", "human_checks"], "properties": {
    "reasoning": STR,
    "rubric": {"type": "array", "items": {"type": "object", "required": ["category", "weight_percent", "coverage", "evidence_in_draft", "gap"], "properties": {
        "category": STR, "weight_percent": {"type": "integer"}, "coverage": {"enum": ["strong", "adequate", "weak"]}, "evidence_in_draft": STR, "gap": STR}}},
    "required_elements": {"type": "array", "items": {"type": "object", "required": ["element", "present", "where"], "properties": {"element": STR, "present": {"type": "boolean"}, "where": STR}}},
    "claims_without_support": {"type": "array", "items": {"type": "object", "required": ["quote", "why"], "properties": {"quote": STR, "why": STR}}},
    "human_checks": {"type": "array", "items": {"type": "object", "required": ["what", "why", "where_to_look"], "properties": {"what": STR, "why": STR, "where_to_look": STR}}}}}


def ledger_view(facts: dict) -> list[dict]:
    """kind 'placeholder': a number, date or link, written as {F:id}.  kind 'paraphrase': a sentence-like fact, to be put in the author's own words and tagged <!-- F:id --> after it."""
    out = []
    for k, f in facts.items():
        v = f["value"]; ph = isinstance(v, (int, float)) or (isinstance(v, str) and (len(v) <= 40 or k.startswith("ref_")))
        out.append({"id": k, "kind": "placeholder" if ph else "paraphrase", "text": f["text"], "value": v, "source": f["source"], "author_must_confirm": f["needs_human"]})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    for k in ("brief", "guideline", "facts", "profile", "samples", "out-dir"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--budget", type=int, default=1000); ap.add_argument("--rounds", type=int, default=2)
    a = ap.parse_args()
    out = Path(getattr(a, "out_dir")); out.mkdir(parents=True, exist_ok=True)
    facts = checks.load_facts(a.facts); profile = json.loads(Path(a.profile).read_text())
    brief = Path(a.brief).read_text(); guide = Path(a.guideline).read_text()
    samples = re.sub(r"\bsubequations\b", "", Path(a.samples).read_text())
    ctx = {"ledger": ledger_view(facts), "style_profile": profile, "word_budget": a.budget}

    outline, _ = call_claude("report-outliner", json.dumps({"guideline": guide, "brief": brief, **ctx}, indent=1), schema=OUTLINE, label="report:outline")
    (out / "outline.json").write_text(json.dumps(outline, indent=1))
    draft, _ = call_claude("report-drafter", json.dumps({"mode": "draft", "outline": outline, "voice_samples": samples, **ctx}, indent=1), label="report:draft")
    draft = draft if isinstance(draft, str) else json.dumps(draft)
    history = []; versions = []
    for rnd in range(a.rounds + 1):
        (out / f"draft_v{rnd + 1}.md").write_text(draft)
        res = checks.check(draft, facts, profile, max_words=a.budget, repo_root=Path(__file__).resolve().parents[2])
        text = checks.render(draft, {**facts, **{b["text"]: {"value": "X"} for b in res["problems"] if b["kind"] == "unknown_fact"}})
        critic, _ = call_claude("report-critic", json.dumps({"draft": text, "style_profile": profile}, indent=1), schema=CRITIC, label=f"report:critic{rnd}")
        audit, _ = call_claude("report-rubric-auditor", json.dumps({"guideline": guide, "brief": brief, "draft": text, "facts_used": [facts[i] for i in checks.used_ids(draft) if i in facts]}, indent=1), schema=AUDIT, label=f"report:audit{rnd}")
        hi = [i for i in critic["issues"] if i["severity"] == "high"]
        history.append({"round": rnd + 1, "program_problems": len(res["problems"]), "critic_high": len(hi), "critic_total": len(critic["issues"]),
                        "weak_rubric": [r["category"] for r in audit["rubric"] if r["coverage"] == "weak"], "metrics": res["metrics"]})
        (out / f"review_v{rnd + 1}.json").write_text(json.dumps({"checks": res, "critic": critic, "audit": audit}, indent=1))
        score = 3 * len(res["problems"]) + 2 * len(hi) + 2 * len([r for r in audit["rubric"] if r["coverage"] == "weak"]) + (6 if res["metrics"]["words"] > a.budget else 0)
        versions.append((score, rnd + 1, draft, res, critic, audit)); history[-1]["score"] = score
        print(json.dumps(history[-1]))
        if rnd == a.rounds or (not res["problems"] and not hi and not history[-1]["weak_rubric"]):
            break
        feedback = {"mode": "revise", "words_now": res["metrics"]["words"], "word_budget_hard": a.budget, "previous_draft": draft, "failed_checks": res["problems"], "critic_issues": critic["issues"], "reads_like_ai": critic["reads_like_ai"],
                    "rubric_gaps": [r for r in audit["rubric"] if r["coverage"] != "strong"], "missing_required": [e for e in audit["required_elements"] if not e["present"]],
                    "unsupported": audit["claims_without_support"], "outline": outline, "voice_samples": samples, **ctx}
        draft, _ = call_claude("report-drafter", json.dumps(feedback, indent=1), label=f"report:revise{rnd + 1}")
        draft = draft if isinstance(draft, str) else json.dumps(draft)

    best = min(versions, key=lambda v: (v[0], -v[1]))
    _, best_round, draft, res, critic, audit = best
    print(f"best version: v{best_round} (score {best[0]}; lower is better)")
    files = render.write_all(draft, a.facts, out / "report")
    used = sorted(checks.used_ids(draft))
    human = ["# What you should check before submitting", "", "Each item says what to check, why, and where to find it.", ""]
    for f in (facts[i] for i in used if i in facts and facts[i]["needs_human"]):
        human.append(f"- **Confirm `{f['id']}` = {f['value']}**: {f['text']}. Source: {f['source']}.")
    for h in audit["human_checks"]:
        human.append(f"- **{h['what']}**: {h['why']} Look at: {h['where_to_look']}.")
    for m in re.finditer(r"\[(AUTHOR TO COMPLETE|NEED FACT)[^\]]*\]", draft):
        human.append(f"- **Fill in or resolve** {m.group(0)}")
    for p in res["problems"]:
        human.append(f"- **Unresolved automatic check ({p['kind']})**: {p['text']} ({p['where'][:100]})")
    for r in critic["reads_like_ai"]:
        human.append(f"- **Sounds like AI; reword in your own voice**: \"{r['quote']}\" ({r['why']})")
    human += ["", "## Every number in the report and where it comes from", "", "| fact | value | meaning | source |", "|---|---|---|---|"]
    human += [f"| {i} | {facts[i]['value']} | {facts[i]['text']} | {facts[i]['source']} |" for i in used if i in facts]
    (out / "HUMAN_CHECKLIST.md").write_text("\n".join(human) + "\n")
    (out / "history.json").write_text(json.dumps(history, indent=1))
    print(json.dumps(files)); print("history:", [(h["round"], h["program_problems"], h["critic_high"]) for h in history])


if __name__ == "__main__":
    main()
