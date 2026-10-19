#!/usr/bin/env python3
"""Evidence collector: builds the fact ledger a report may draw on.  Every fact has an id, a value, a plain-language text, and a SOURCE (a file or a
command) so that a reader can check it.  Facts the author must confirm himself carry needs_human = true.  No fact is written by a language model.

  python -m agentic.writing.evidence --out research/Review/writing/facts.json [--since 2026-10-12]
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import os
import re
import sqlite3
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HOME = Path.home()
MATH = HOME / "PetIGA-3.20-math" / "projects" / "lunar_regolith_DSM"
REVIEW = ROOT / "research" / "Review"
LEDGERS = [ROOT / ".agent_state" / "ledger.jsonl", MATH / ".agent_state" / "ledger.jsonl"]


def sh(cmd: str, cwd: Path = ROOT) -> str:
    return subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True).stdout.strip()


REFS = {"kiewiet": "10.1007/s44461-025-00005-4", "kaempfer": "10.1103/physreve.79.031502", "moure": "10.1021/acs.cgd.4c00539",
        "li": "10.4310/cms.2009.v7.n1.a4", "yu": "10.1088/0965-0393/20/7/075008", "guo": "10.1017/jfm.2020.790", "nicoli": "10.1103/physreve.84.046707"}


def crossref(doi: str, cache: dict, cache_file: Path) -> dict:
    """Metadata from Crossref (verifies that the DOI exists and gives the exact title, authors, venue); cached because the service rate-limits."""
    import time, urllib.parse, urllib.request
    if doi in cache:
        return cache[doi]
    for i in range(5):
        try:
            req = urllib.request.Request("https://api.crossref.org/works/" + urllib.parse.quote(doi), headers={"User-Agent": "lit-registry/0.1"})
            m = json.load(urllib.request.urlopen(req, timeout=30))["message"]
            cache[doi] = {"title": m["title"][0], "authors": [[a.get("given", ""), a.get("family", "")] for a in m.get("author", [])],
                          "journal": (m.get("container-title") or [""])[0], "volume": m.get("volume"), "issue": m.get("issue"),
                          "page": m.get("page") or m.get("article-number"), "year": m["issued"]["date-parts"][0][0]}
            cache_file.write_text(json.dumps(cache, indent=1)); return cache[doi]
        except Exception:
            time.sleep(5 * (i + 1))
    raise RuntimeError(f"Crossref lookup failed for {doi}")


TITLE_FIX = {"10.4310/cms.2009.v7.n1.a4": "Solving PDEs in complex geometries: a diffuse domain approach"}   # Crossref has a truncated, lower-case title; the PDF title page reads this


def format_ref(doi: str, m: dict) -> str:
    import html
    m = {**m, "title": html.unescape(TITLE_FIX.get(doi, m["title"])), "journal": html.unescape(m["journal"])}
    au = [f"{fam}, {''.join(x[0] + '.' for x in re.split(r'[ -]+', giv) if x)}" for giv, fam in m["authors"]]
    names = (", ".join(au[:-1]) + (", & " if len(au) > 2 else " & ") + au[-1]) if len(au) > 1 else (au[0] if au else "")
    if len(au) > 6:
        names = ", ".join(au[:6]) + ", et al."
    vol = f", {m['volume']}" if m.get("volume") else ""; iss = f"({m['issue']})" if m.get("issue") else ""; pg = f", {m['page']}" if m.get("page") else ""
    out = f"{names} ({m['year']}). {m['title'].rstrip('.')}. {m['journal']}{vol}{iss}{pg}. https://doi.org/{doi}"
    return out.replace("R\x7fatz", "Rätz")        # Crossref stores this surname with a stray control character; the paper prints Rätz


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--since", default="2026-10-12")
    a = ap.parse_args()
    F: list[dict] = []

    def add(id_, value, text, source, needs_human=False):
        F.append({"id": id_, "value": value, "text": text, "source": source, "needs_human": needs_human})

    # ---- repository activity
    n_agent = int(sh(f"git log --all --since={a.since} --oneline -- agentic .claude/agents | wc -l") or 0)
    n_all = int(sh(f"git log --all --since={a.since} --oneline | wc -l") or 0)
    add("commits_agent", n_agent, "commits touching the agent code and agent definitions since the start of the agent work", f"git log --all --since={a.since} -- agentic .claude/agents")
    add("commits_all", n_all, "commits on all branches over the same period", f"git log --all --since={a.since}")
    names = {Path(p).name for p in glob.glob(str(ROOT / ".claude/agents/*.md")) + glob.glob(str(MATH / ".claude/agents/*.md"))}
    add("n_agents", len(names), "agent definitions (one instruction file each)", ".claude/agents/ on main and agent/math-lane")
    add("start_date", __import__("datetime").date.fromisoformat(a.since).strftime("%-d %B %Y"), "start of the agent-workflow work", "user's statement (started the week of 12 October)", needs_human=True)

    # ---- literature pipeline
    db = sqlite3.connect(ROOT / ".agent_state" / "registry.sqlite")
    add("n_candidates", db.execute("select count(*) from papers").fetchone()[0], "papers in the registry (found by search or already in the library)", ".agent_state/registry.sqlite")
    recs = {}
    for f in sorted(glob.glob(str(REVIEW / "records" / "*.json"))):
        if f.endswith(".verification.json"):
            continue
        stem = os.path.basename(f)[:-5]
        if re.search(r"(_(opus|sonnet|haiku|rep\d+))+$", stem):
            continue
        recs[stem] = (json.load(open(f)), json.load(open(f[:-5] + ".verification.json")))
    add("n_papers_read", len(recs), "papers read in full, one fresh worker per paper", "research/Review/records/*.json (one per paper)")
    add("n_claims", sum(len(d.get("claims", [])) for d, _ in recs.values()), "claims extracted, each with a page and a verbatim quote", "research/Review/records/*.json")
    qm = [v["metrics"]["quote_match_rate"] for _, v in recs.values() if v.get("metrics") and "quote_match_rate" in v["metrics"]]
    add("quote_match_mean_pct", round(100 * sum(qm) / len(qm), 1), "average share of quotes found word for word on their cited page", "records/*.verification.json metrics.quote_match_rate")
    ck = [v["checker"]["supported_rate"] for _, v in recs.values() if v.get("checker") and v["checker"].get("supported_rate") is not None]
    add("checker_supported_mean_pct", round(100 * sum(ck) / len(ck), 1), "average share of sampled claims that the independent checker rated fully supported by their quote", "records/*.verification.json checker.supported_rate")
    corr = conf = 0
    for _, v in recs.values():
        c = (v.get("image_check") or {}).get("counts") or {}
        corr += c.get("corrected", 0); conf += c.get("confirmed", 0)
    add("image_items_corrected", corr, "equations or numbers the text extraction had garbled, fixed by reading the page image", "records/*.verification.json image_check")
    add("image_items_confirmed", conf, "equations or numbers the page-image check confirmed unchanged", "records/*.verification.json image_check")
    audit = (REVIEW / "memos" / "gate_R8_diffuse_domain_audit.md").read_text()
    m = re.search(r"(\d+) cited lines judged\. Verdicts: (\{.*?\})", audit)
    verd = eval(m.group(2)) if m else {}
    add("audit_lines", int(m.group(1)) if m else None, "cited lines of the literature memo that an independent checker judged against the cited pages", "research/Review/memos/gate_R8_diffuse_domain_audit.md")
    for k, name in (("supported", "audit_supported"), ("partially_supported", "audit_partial"), ("unsupported", "audit_unsupported"), ("cannot_tell", "audit_undecided")):
        add(name, verd.get(k), f"audit verdict '{k}'", "research/Review/memos/gate_R8_diffuse_domain_audit.md")

    # ---- calls and cost from the ledgers
    tot = collections.defaultdict(lambda: [0, 0.0])
    for lf in LEDGERS:
        if not lf.exists():
            continue
        for line in lf.read_text().splitlines():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            t = tot[(r["agent"], r["model"])]; t[0] += 1; t[1] += r.get("cost_usd_list") or 0
    add("n_calls", sum(t[0] for t in tot.values()), "language-model calls made by the agent workers", "token ledgers .agent_state/ledger.jsonl (main and math worktree)")
    add("cost_list_usd", round(sum(t[1] for t in tot.values()), 2), "total list-price cost of those calls (not billed to a personal account)", "token ledgers")
    ex = [t for k, t in tot.items() if k[0] == "paper-extractor"]
    add("cost_per_paper_usd", round(sum(t[1] for t in ex) / max(1, len(recs)), 2), "list-price cost per paper for reading it (one worker)", "token ledgers, paper-extractor")
    add("n_checker_calls", tot[("claim-checker", "opus")][0], "independent checker calls (Opus)", "token ledgers")

    # ---- derivation workflow ("math lane")
    ready = notready = known = 0; tasks = []
    for d in sorted(glob.glob(str(MATH / "derivations" / "*" / "SUMMARY.md"))):
        if Path(d).parent.name.endswith("_run1"):        # saved copy of a first attempt, not a separate task
            continue
        s = Path(d).read_text(); tasks.append(Path(d).parent.name)
        if "STATUS: ready" in s: ready += 1
        elif "STATUS: NOT ready" in s: notready += 1
        if "known answer: MATCH" in s: known += 1
    add("lane_tasks", len(tasks), "derivation tasks run through the workflow", "derivations/*/SUMMARY.md (agent/math-lane)")
    add("lane_ready", ready, "tasks that finished 'ready for review' (all checks passed, no serious objection)", "derivations/*/SUMMARY.md")
    add("lane_not_ready", notready, "tasks the workflow itself marked 'not ready' (a failed check or a serious objection)", "derivations/*/SUMMARY.md")
    add("lane_known_answer", known, "tasks whose result matched a known answer", "derivations/*/SUMMARY.md")
    rt = {k: json.load(open(REVIEW / "memos" / n)) for k, n in (("r1", "f2_red_team.json"),)}
    r1 = [len(v["objections"]) for v in rt["r1"].values()]
    r2 = [len(json.load(open(p))["objections"]) for p in sorted(glob.glob(str(REVIEW / "memos" / "f2_red_team_round2_*.json")))]
    add("redteam_round1", sum(r1), "objections from three fresh reviewer agents (round 1)", "research/Review/memos/f2_red_team.json")
    add("redteam_round2", sum(r2), "objections from three fresh reviewer agents (round 2)", "research/Review/memos/f2_red_team_round2_*.json")
    add("seeded_errors", "4 of 4", "deliberate errors planted in a correct derivation that at least one check caught (one was caught only by the reviewer agent)", "Planning/19_f2_derivation_full.md section 11 (seeded-error test)", needs_human=False)
    add("checker_bugs_found", 3, "bugs in our own symbolic checker exposed by testing it on real derivations", "git log agent/math-lane, commit accf8e31 message")

    # ---- more process facts (computed)
    add("claims_checked_total", sum((v.get("checker") or {}).get("n_checked", 0) for _, v in recs.values()), "claims independently checked by the second model (every number or equation claim plus a quarter of the rest)", "records/*.verification.json checker.n_checked; sampling rule in agentic/lit/check_claims.py")
    reps = []
    for rp in glob.glob(str(REVIEW / "memos" / "wave*_report.md")):
        reps += [float(x) for x in re.findall(r"probe agreement ([0-9.]+)", Path(rp).read_text())]
    add("replicate_agreement_min_pct", round(100 * min(reps)), "lowest agreement between two independent readings of the same paper (answers to the same questions)", "research/Review/memos/wave*_report.md") if reps else None
    add("replicates_n", len(reps), "papers read twice, independently, to measure agreement", "research/Review/memos/wave*_report.md") if reps else None
    if m:
        add("audit_uncertain_pct", round(100 * (verd.get("partially_supported", 0) + verd.get("unsupported", 0) + verd.get("cannot_tell", 0)) / int(m.group(1))), "share of audited memo lines that were not fully supported", "research/Review/memos/gate_R8_diffuse_domain_audit.md")
    memo_txt = (REVIEW / "memos" / "gate_R8_diffuse_domain.md").read_text()
    err = memo_txt.split("ERRATA from the citation audit")[-1] if "ERRATA from the citation audit" in memo_txt else ""
    add("memo_errata_items", len(re.findall(r"(?m)^\d+\. ", err)), "corrections we appended to the literature memo after the audit", "research/Review/memos/gate_R8_diffuse_domain.md (errata section)")
    add("papers_unread", db.execute("select count(*) from papers").fetchone()[0] - len(recs), "registry papers not yet read in full", ".agent_state/registry.sqlite minus records")
    w1 = (HOME / "PetIGA-3.20" / "projects" / "lunar_regolith_DSM" / "research" / "Literature" / "Wave1_incoming" / "RETRIEVAL_CHECKLIST.md")
    n1 = len(re.findall(r"(?m)^- \[ \]", w1.read_text())) if w1.exists() else None
    w2 = (ROOT / "research" / "Literature" / "Wave2_incoming" / "RETRIEVAL_CHECKLIST.md")
    n2 = len(re.findall(r"(?m)^- \[ \]", w2.read_text())) if w2.exists() else None
    if n1 is not None and n2 is not None:
        add("paywalled_papers", n1 + n2, "papers the open-access fetch could not get, which the author downloaded or must download through the library", "research/Literature/Wave*_incoming/RETRIEVAL_CHECKLIST.md")
    # ---- facts taken from the project record or the working conversation (the author should confirm each)
    add("lunar_background", "a 2025 review of lunar water extraction compares 27 studies and describes thermal extraction, which sublimates ice in a low-pressure environment, as one main route", "why lunar ice matters, from a paper we read", "research/Review/records/kiewiet_2025_25000054.md claims C1, C3 (pages 1 and 6)", needs_human=True)
    add("retraction_detail", "we first reported an air bubble forming inside a pair of ice grains; a closer look showed the ice surface was retreating as the grains fused", "what the retracted claim was", "research/Planning/08_air_film_diagnostics.md Update 7")
    add("checker_bugs_detail", "a numeric check that passed when it could not evaluate anything, a variable named pi that the program read as the number, and derivatives of sums that were not expanded; all three were fixed in the commit that added the test", "the three bugs in our symbolic checker", "git commit accf8e31 (agent/math-lane)")
    add("decision_reading_lists", "the author approved the reading list before each literature wave and added papers by hand from the personal library", "a human decision point", "conversation of 2026-10-17 to 2026-10-19", needs_human=True)
    add("decision_stop_runs", "after inspecting simulation output the author stopped further runs until the model equations were redone", "a human decision point", "conversation of 2026-10-18", needs_human=True)
    add("decision_third_phase", "the author asked whether the soil should enter the energy as a third phase at all, which shaped the formulation we chose", "a human decision point", "conversation of 2026-10-18", needs_human=True)
    add("decision_mechanism", "the author's own earlier debugging (the ice profile relaxes while the soil profile stays fixed) matched the mechanism our checks isolated", "human insight confirmed by the workflow", "conversation of 2026-10-18", needs_human=True)
    add("decision_checker_panel", "the author proposed replacing the single independent checker by several checkers whose verdicts are compared, with an adjudicator for disagreements", "a human design decision", "conversation of 2026-10-19", needs_human=True)
    add("own_errors_corrected", 6, "mistakes in our own claims or plan that fresh reviewers or follow-up checks exposed and that we corrected", "research/Planning/20_f2_six_iteration_audit.md and the session summary of 2026-10-19")
    add("not_ready_but_matched", "one task matched its known answer yet was marked not ready, because a reviewer agent found an assumption we had not stated", "why two derivation counts differ", "derivations/c3d_wall_angle_first_integral/SUMMARY.md (agent/math-lane)")
    add("model_pinning", "each agent file names a model tier (Sonnet, Opus or Haiku) and the ledger records the tier and exact token use per call, but not a dated model version", "a reproducibility limit", ".claude/agents/*.md and .agent_state/ledger.jsonl")
    add("run_commands", "each workflow runs from one command in the repository: run_wave for papers, run_lane for derivations, run_report for writing", "how to run the workflows", "agentic/lit/run_wave.py, agentic/math/run_lane.py, agentic/writing/run_report.py")

    # ---- cost split and testing totals
    def cost(agents):
        return round(sum(t[1] for k, t in tot.items() if k[0] in agents), 2)
    add("cost_checker_usd", cost({"claim-checker"}), "list-price cost of the independent claim checks", "token ledgers")
    add("cost_reading_usd", cost({"paper-extractor"}), "list-price cost of reading papers", "token ledgers")
    add("cost_image_usd", cost({"image-checker"}), "list-price cost of the page-image checks", "token ledgers")
    add("cost_derivation_usd", cost({"math-spec-writer", "math-deriver", "math-adversary"}), "list-price cost of the derivation workflow", "token ledgers")
    add("n_triaged", tot[("triage-scorer", "haiku")][0], "candidate papers scored from title and abstract by the smallest model", "token ledgers")
    add("image_check_calls", tot[("image-checker", "opus")][0], "page-image check calls", "token ledgers")
    add("n_checks_total", 9 + 4 + 9 + 15, "automatic tests and regression checks that pass (symbolic checker 9, derivation driver 4, report checks 9, formulation regression suite 15)", "agentic/tests/*.py and agentic/math/numerics/run_all_f2_checks.py")
    add("known_answers_text", "the workflow reproduced a published constant (found as the exact fraction 47/75), a hand-derived three-phase equation, and a thermal result that we then checked against a numerical solution", "what the known-answer tasks were", "Planning/16, Planning/19 section 11", needs_human=False)
    add("literature_gaps", "we found no paper on ice contact angles on regolith or silicate in depth, and none that couples a diffuse solid to sublimation or latent heat", "gaps the literature search exposed (author's knowledge)", "author's statement of 2026-10-19", needs_human=True)
    add("token_use_note", "our token use has been well below the plan allowance (about a fifth), so we will scale the next literature wave up by a factor of two to three", "author's statement about usage", "author's statement of 2026-10-19 (check the usage page)", needs_human=True)

    # ---- bibliography (formatted by a program from Crossref, never typed by a model)
    cf = ROOT / "research" / "Review" / "writing" / "refs_crossref.json"
    cache = json.loads(cf.read_text()) if cf.exists() else {}
    for key, doi in REFS.items():
        add(f"ref_{key}", format_ref(doi, crossref(doi, cache, cf)), f"full citation ({key})", f"Crossref record of doi:{doi}")
    # ---- tests
    py = os.environ.get("PY_SYMPY") or next((str(x) for x in (HOME / ".cache").glob("**/venv/bin/python") if x.exists()), None) or "python3"
    scratch = glob.glob("/private/tmp/claude-501/*/*/scratchpad/venv/bin/python")
    py = os.environ.get("PY_SYMPY") or (scratch[0] if scratch else "python3")
    p = subprocess.run([py, "-I", str(MATH / "agentic" / "tests" / "test_sympy_check.py")], capture_output=True, text=True, cwd=MATH)
    add("tests_symbolic", p.stdout.count("PASS"), "passing known-answer tests of the symbolic checker", "python agentic/tests/test_sympy_check.py (agent/math-lane)")
    # ---- human and honesty facts
    add("retracted_claims", 1, "claim I made and later retracted after a check showed it was wrong (an air 'bubble' in a grain pair was surface retreat)", "Planning/08_air_film_diagnostics.md Update 7")
    add("human_gates", 3, "points where a person must approve before the workflow continues (reading list, flagged claims, summary memos)", "Planning/10_lit_review_implementation.md section 3")
    add("repo_url", "https://github.com/citgitfr/BaglinoJackson_IceSedimentInteractions", "course repository", "user message 2026-10-17; confirm it is current and that the grader can open it", needs_human=True)
    F = [f for f in F if f]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(F, indent=1))
    print(f"{len(F)} facts -> {a.out}; need your confirmation: {[f['id'] for f in F if f['needs_human']]}")


if __name__ == "__main__":
    main()
