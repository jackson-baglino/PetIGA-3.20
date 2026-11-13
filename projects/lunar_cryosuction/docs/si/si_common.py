"""Shared helpers for the supplementary-methods documents (docs/si/).

Reuses the run-plan Word helpers (docs/run_plan/docx_helpers.py: styles, LaTeX ->
Word equations). Import this module, then build with doc, H, P, T and EQ.
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "run_plan"))
os.environ.setdefault("RUNPLAN_DATE", "2026-11-13")          # footer date
os.environ.setdefault("DOC_FOOTER", "Supplementary Methods, paper 2 (draft)")

from docx_helpers import (  # noqa: E402,F401
    GRAYT, H, Inches, P, Pt, T, WD_ALIGN_PARAGRAPH, add_math, doc, rich)

OUTDIR = Path("/Users/jacksonbaglino/Library/Mobile Documents/com~apple~CloudDocs/Documents/"
              "Research/Projects/lunar_regolith_DSM/SupplementalInfo_paper2")

_n = [0]


def EQ(latex, prefix="S"):
    """Centered display equation with a right-hand number, (S1), (S2), ..."""
    _n[0] += 1
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(3)
    p.paragraph_format.space_after = Pt(3)
    add_math(p, latex, size=10.5)
    r = p.add_run(f"\t({prefix}{_n[0]})")
    r.font.size = Pt(10)
    return f"{prefix}{_n[0]}"


def TBL(header, rows, widths, size=8.5):
    """Table with a light header, so Word equations in it stay legible."""
    return T(header, rows, widths, size=size, head_fill="DCE3EE", head_color=(0x1F, 0x3A, 0x5F))


def title(text, sub):
    doc.add_paragraph(text, style="Title")
    P(sub, size=9, color=GRAYT, after=6)


def save(name):
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / name
    doc.save(out)
    print("wrote", out)
    return out
