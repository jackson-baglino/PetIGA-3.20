#!/usr/bin/env python3
"""Style profile of the author's own writing, measured (not guessed).  Reads a LaTeX manuscript, strips markup, and writes numbers only
(sentence length, pronoun use, passive rate, paragraph length).  Optionally writes a few sample paragraphs for the drafter to imitate; keep that
file outside the repository (it is unpublished writing).

  python -m agentic.writing.style --tex manuscript.tex --out agentic/writing/style_profile.json [--samples research/Review/writing/voice_samples.md]
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path


def prose(tex: str) -> str:
    b = tex[tex.index("\\begin{document}"):] if "\\begin{document}" in tex else tex
    b = re.sub(r"(?m)(?<!\\)%.*$", "", b)
    b = re.sub(r"\\begin\{(equation|align|figure|table)\*?\}.*?\\end\{\1\*?\}", "", b, flags=re.S)
    b = re.sub(r"\\(cite[a-z]*|ref|eqref|label)\{[^}]*\}", "", b)
    b = re.sub(r"\$[^$]*\$", "X", b)
    b = re.sub(r"\\(section|subsection|subsubsection)\*?\{([^}]*)\}", r"\n\n## \2\n", b)
    b = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?\{([^}]*)\}", r"\2", b)
    b = re.sub(r"\\[a-zA-Z]+", "", b); b = re.sub(r"[{}]", "", b)
    return re.sub(r"\n{3,}", "\n\n", b)


def sentences(text: str) -> list[str]:
    return [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if len(s.split()) > 3]


def measure(text: str) -> dict:
    ss = sentences(text); L = sorted(len(s.split()) for s in ss)
    paras = [p for p in re.split(r"\n\s*\n", text) if len(p.split()) > 25 and not p.startswith("##")]
    n = max(1, len(ss))
    return {"sentences": len(ss), "words": len(text.split()), "mean_sentence_words": round(statistics.mean(L), 1), "median_sentence_words": statistics.median(L),
            "p90_sentence_words": L[int(0.9 * (len(L) - 1))], "share_over_35_words": round(sum(x > 35 for x in L) / n, 3),
            "we_per_100_sentences": round(100 * len(re.findall(r"\bwe\b", text, re.I)) / n, 1), "I_per_100_sentences": round(100 * len(re.findall(r"\bI\b", text)) / n, 1),
            "passive_per_100_sentences": round(100 * len(re.findall(r"\b(is|are|was|were|been|be)\s+\w+ed\b", text)) / n, 1),
            "mean_paragraph_words": round(statistics.mean(len(p.split()) for p in paras), 1) if paras else None}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tex", required=True); ap.add_argument("--out", required=True); ap.add_argument("--samples")
    a = ap.parse_args()
    text = prose(Path(a.tex).read_text())
    prof = measure(text)
    Path(a.out).write_text(json.dumps(prof, indent=1)); print(prof)
    if a.samples:
        paras = [p.strip() for p in re.split(r"\n\s*\n", text) if 60 < len(p.split()) < 170 and not p.startswith("##") and "X" not in p[:200]]
        pick = paras[:2] + paras[len(paras) // 2:len(paras) // 2 + 2]
        Path(a.samples).write_text("# Paragraphs in the author's own voice (for imitation only; do not copy sentences)\n\n" + "\n\n---\n\n".join(pick) + "\n")


if __name__ == "__main__":
    main()
