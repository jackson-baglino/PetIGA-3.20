"""Second seed: the Moure & Fu 2024 Supplemental Information (Table S1, S3: units of every model parameter), image-checked per paper.
Project notation is used for the symbol; the printed symbol stays in the meaning. Temperatures printed in degC are NOT added (see the report to the user).
Existing symbols are never touched."""
import json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(ROOT))
from agentic.library.units import Curator, Frozen
NAME = {"βsol": "beta_sol", "dj": "d_j", "uj": "u_j", "vn,j": "v_n_j", "Mj": "M_j", "αj": "alpha_j", "Λ": "Lambda", "σiw": "gamma_iw", "σia": "gamma_ia", "σwa": "gamma_wa",
        "αnucl": "alpha_nucl", "ρw": "rho_w", "cp,a": "c_p_a", "cp,i": "c_p_i", "cp,w": "c_p_w", "Ka": "K_a", "Ki": "K_i", "Kw": "K_w", "Lsol": "L_sol", "Lsub": "L_sub",
        "Dv0": "D_v0", "ξv": "xi_v", "ξT": "xi_T", "χ": "chi"}
SKIP = {"Tmelt", "Tnucl", "Pa", "ε", "ρa", "ρi", "L", "D*iw"}      # degC-printed; or already in the standard / not on the page
c = Curator(); added, skipped = [], []
d = json.loads((ROOT / "research/Review/units/candidates_moure_fu_2024_SI.json").read_text())
for it in d["items"]:
    s = it["symbol"]
    if s in SKIP or it["image"]["status"] != "confirmed" or not it["quote_found"] or s not in NAME:
        skipped.append((s, it["image"]["status"], it["quote_found"])); continue
    try:
        c.add(NAME[s], f"{it['meaning']} (printed as {s})", it["unit"], it["kind"], paper_id="moure_fu_2024_SI", page=it["page"], quote=it["quote"], verified="image-checked (Opus) against the page image")
        added.append(NAME[s])
    except Frozen:
        skipped.append((s, "already in the standard", True))
# user-origin entry: the user's own instruction
c.add("alpha", "sublimation attachment rate in the ice source term, alpha = lam/tau_sub (Moure and Fu alpha_j; SI Table S1)", "s^-1", "parameter", origin="user",
      user_instruction="user, 2026-10-21: 'The unit of \\alpha should be 1/s.' (see Moure & Fu 2024 SI p.7, alpha_j)")
added.append("alpha")
c.add("M0", "mobility of the phase-field equation (mob_sub); same role and unit as M_j of Moure and Fu", "m s^-1", "parameter", origin="dimensional-analysis",
      basis="Moure & Fu 2024 SI p.7 gives M_j in m/s; w psi_t [1/s] = M0 [m/s] * G [1/m] (doc 19 units check)")
added.append("M0")
print("added", len(added), added); print("skipped", skipped)
