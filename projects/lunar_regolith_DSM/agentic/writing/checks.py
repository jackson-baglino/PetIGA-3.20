"""Deterministic checks for a report draft.  No language model is involved, so an agent cannot talk its way past them.
A draft states every number as a placeholder {F:fact_id}; `render` fills in the value from the fact ledger, so a number cannot drift from its source.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
PLACE = re.compile(r"\{F:([a-z0-9_]+)\}")
COMMENT = re.compile(r"<!--\s*F:([a-z0-9_, ]+?)\s*-->")     # traceability tag after a paraphrased sentence: <!-- F:id1, id2 -->


def used_ids(draft: str) -> set:
    ids = set(PLACE.findall(draft))
    for g in COMMENT.findall(draft):
        ids |= {x.strip() for x in g.split(",") if x.strip()}
    return ids
WORDS_PER_PAGE = 500


def load_facts(path) -> dict:
    return {f["id"]: f for f in json.loads(Path(path).read_text())}


def render(draft: str, facts: dict) -> str:
    def sub(m):
        fid = m.group(1)
        if fid not in facts:
            raise KeyError(f"unknown fact id {fid!r}")
        v = facts[fid]["value"]
        return f"{v:,}" if isinstance(v, int) and v >= 10000 else str(v)
    return COMMENT.sub("", PLACE.sub(sub, draft))


def sentences(text: str):
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", re.sub(r"\s+", " ", text)) if len(s.split()) > 3]


def check(draft: str, facts: dict, profile: dict | None = None, *, max_words: int = 1000, repo_root: Path | None = None) -> dict:
    """Return {'problems': [...], 'metrics': {...}}.  Each problem: {'kind','text','where'}."""
    P = []
    lex = json.loads((HERE / "jargon.json").read_text())
    bad_ids = sorted({m for m in used_ids(draft) if m not in facts})
    for b in bad_ids:
        P.append({"kind": "unknown_fact", "text": b, "where": "placeholder"})
    try:
        text = render(draft, {**facts, **{b: {"value": "X"} for b in bad_ids}})
    except KeyError:
        text = draft
    text_main = text.split("## References")[0]
    body = re.sub(r"(?m)^(#+ .*|\|.*|<<<PAGEBREAK>>>)$", "", text_main)
    words = len(body.split())
    if words > max_words:
        P.append({"kind": "too_long", "text": f"{words} words; budget {max_words} (about {max_words / WORDS_PER_PAGE:.1f} pages)", "where": "document"})
    # raw numbers in the draft (not from the ledger): flag everything except small counts, years, dates, list numbers and the AI-declaration appendix
    stripped = COMMENT.sub("", PLACE.sub("", draft))
    stripped = re.sub(r"(?m)^\s*\d+[.)]\s", "", stripped)
    for m in re.finditer(r"(?<![\w./`-])(\d[\d,]*\.?\d*)(?![\w`/-])", stripped):
        tok = m.group(1).rstrip(".,")
        if re.fullmatch(r"20\d\d|[0-9]|1[0-9]|2[0-9]|3[01]", tok):   # years, day numbers and small counts are allowed
            continue
        ctx = stripped[max(0, m.start() - 40):m.end() + 20].replace("\n", " ")
        P.append({"kind": "raw_number", "text": tok, "where": ctx})
    # sentence length
    cap = 35
    for s in sentences(body):
        n = len(s.split())
        if n > cap:
            P.append({"kind": "long_sentence", "text": f"{n} words: {s[:110]}...", "where": "sentence"})
    # first person
    n_I = len(re.findall(r"\bI\b", body)); n_we = len(re.findall(r"\bwe\b", body, re.I))
    if n_I > 0:
        P.append({"kind": "voice", "text": f"'I' used {n_I} time(s); the author writes 'we'", "where": "document"})
    if n_we == 0:
        P.append({"kind": "voice", "text": "no 'we' found", "where": "document"})
    # jargon: first use must be glossed (parenthesis, 'that is', dash, 'i.e.') within 90 characters
    for term, plain in lex["physics_jargon"].items():
        m = re.search(re.escape(term), body, re.I)
        if m and re.search(r"\(\s*" + re.escape(term) + r"\s*\)", body, re.I):      # the gloss comes first and the term is named in brackets
            m = None
        if m:
            after = body[m.end():m.end() + 90]
            if not re.search(r"\(|\s[—-]{1,2}\s|\bthat is\b|\bi\.e\.|\bmeaning\b|\bwhich means\b|:", after[:70]):
                P.append({"kind": "jargon_unglossed", "text": term, "where": body[max(0, m.start() - 40):m.end() + 40].replace("\n", " "), "plain": plain})
    for w in lex["ai_sounding_words"]:
        if re.search(rf"\b{re.escape(w)}", body, re.I):
            P.append({"kind": "ai_sounding", "text": w, "where": "word list"})
    # structure
    for need, rx in (("AI Use Declaration", r"AI Use Declaration"), ("progress: done", r"(?i)\bdone\b"), ("progress: remaining", r"(?i)remaining"), ("progress: bottlenecks", r"(?i)bottleneck"),
                     ("next steps", r"(?i)next (steps|two weeks)"), ("repository link", r"github\.com|\{F:repo_url\}")):
        if not re.search(rx, draft):
            P.append({"kind": "missing_section", "text": need, "where": "document"})
    # paths in backticks must exist (relative to the repo root) when a root is given
    if repo_root:
        for m in re.finditer(r"`([\w./-]+/[\w./-]+)`", draft):
            if not (repo_root / m.group(1)).exists():
                P.append({"kind": "missing_path", "text": m.group(1), "where": "backticked path"})
    # facts the author must confirm
    used = used_ids(draft)
    human = [f for fid, f in facts.items() if f.get("needs_human") and fid in used]
    sl = sentences(body)
    metrics = {"words": words, "approx_pages": round(words / WORDS_PER_PAGE, 2), "sentences": len(sl),
               "mean_sentence_words": round(sum(len(s.split()) for s in sl) / max(1, len(sl)), 1), "we": n_we, "I": n_I,
               "facts_used": len(used), "facts_needing_confirmation": [h["id"] for h in human]}
    if profile:
        metrics["target_mean_sentence_words"] = profile.get("mean_sentence_words")
    return {"problems": P, "metrics": metrics}
