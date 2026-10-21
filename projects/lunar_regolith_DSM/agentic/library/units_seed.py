"""One-off seed of the units standard from the image-checked candidates of Kaempfer & Plapp 2009 and Moure & Fu 2024.
Only items that the image check confirmed AND whose quote was found verbatim are added; everything else is reported, not added.
Existing symbols are never touched (first source wins; Kaempfer & Plapp first)."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from agentic.library.units import Curator, Frozen

GREEK = {"ρ": "rho", "κ": "kappa", "β": "beta", "γ": "gamma", "Δ": "Delta", "ε": "eps"}
SKIP = {"∇T", "∇ρv"}                       # gradients are derived by calculus, not conventions of their own
MOURE_SKIP = {"T", "T0", "TB", "Tnucl", "Tmelt", "Dv0", "Pvs_J", "water film thickness"}   # degC in the paper / quote not verbatim / not on page: see the report
def ascii_sym(s):
    for g, a in GREEK.items(): s = s.replace(g, a)
    return s
def expand(sym):
    return [x.strip() for x in sym.split(",")]
def sub(sym):                              # rhoi -> rho_i, kappaa -> kappa_a, Dv -> D_v
    for base in ("rho", "kappa"):
        if sym.startswith(base) and len(sym) == len(base) + 1: return f"{base}_{sym[-1]}"
    if sym in ("Dv", "Lsg", "Rda", "Rv", "DSTP", "TSTP", "xs", "vn", "tD", "tn", "jT", "jv"):
        return {"Dv": "D_v", "Lsg": "L_sg", "Rda": "R_da", "Rv": "R_v", "DSTP": "D_STP", "TSTP": "T_STP", "xs": "x_s", "vn": "v_n", "tD": "t_D", "tn": "t_n", "jT": "j_T", "jv": "j_v"}[sym]
    return sym
c = Curator(); added, skipped = [], []
for paper in ("kaempfer_plapp_2009", "moure_fu_2024"):
    d = json.loads((ROOT / f"research/Review/units/candidates_{paper}.json").read_text())
    for it in d["items"]:
        if it["symbol"] in SKIP or (paper == "moure_fu_2024" and it["symbol"] in MOURE_SKIP):
            skipped.append((paper, it["symbol"], "not added: see report")); continue
        if it["image"]["status"] != "confirmed" or not it["quote_found"]:
            skipped.append((paper, it["symbol"], f"{it['image']['status']}, quote_found={it['quote_found']}")); continue
        for s in expand(it["symbol"]):
            name = sub(ascii_sym(s))
            try:
                c.add(name, f"{it['meaning']} (printed as {s})", it["unit"], it["kind"], paper_id=paper, page=it["page"], quote=it["quote"], verified="image-checked (Opus) against the page image")
                added.append(name)
            except Frozen:
                skipped.append((paper, name, "already in the standard"))
print("added", len(added), added); print("skipped", len(skipped))
for s in skipped: print("  ", s)
