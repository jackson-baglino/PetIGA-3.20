# flinComputationsSpecificSurface2014 — Flin, Lesaffre, Dufour, Gillibert, Hasan, Rolland du Roscoat, Cabanes & Pugliese, "On the computations of specific surface area and specific grain contact area from snow 3D images", Physics and Chemistry of Ice 2010 (ed. Furukawa), Hokkaido University Press, 2011, pp. 321–328

**Bib problem:** the key and bib year say 2014, but the PDF header reads "Physics and Chemistry of Ice 2010 ... Hokkaido University Press, 2011, Sapporo" (every page); other papers in the library cite it as Flin et al. (2011) (e.g., Calonne 2014 ref list; Bouvet 2022). No DOI. Fix year to 2011.

**Summary.** Compares three SSA estimators on binary μCT images of 10 snow samples (PP to MF): stereology (ST, line intercepts), marching cubes (MC), voxel projection (VP); tests resolution (voxel doubling) and REV; defines Specific Grain Contact Area (SGCA = SSA_tot − SSA) using curvature-driven grain segmentation.

**Findings relevant to the manuscript**
- Methods are for binary voxel images (p.3): ST, MC, VP. "There is no ideal numerical method to estimate the SSA" (p.7); MC "always overestimates"; VP sensitive to resolution.
- SSA REV p.5: "the REV seems to be attained with a sufficient precision before the side of the image cube reaches 2.5 mm".
- SGCA p.6–7: bonding area between grains, "data that is not available from the sole classical SSA. This could be an interesting parameter to relate with thermal or mechanical properties"; aged (low-SSA) snow has high SGCA fraction.

**Claims it supports / contradicts / qualifies**
- Outline L799–802 ("SSA: how it is measured from φ ... possibly Flin et al. 2014 (check that the method matches)"): **does not match**. Flin et al. compute SSA on binary voxel images by ST/MC/VP; the manuscript computes SSA from a continuous phase field φ (presumably the length of the φ = 0.5 contour or ∫|∇φ| per ice mass). Do not cite Flin for the method; at most cite it for the definition "total ice–air surface area per unit ice mass" (p.3, Eq. 1 uses surface/mass). The phase-field estimator should be defined in the text without this citation, or cite Bouvet 2022 (who also compute SSA on phase-field output but after binarization, voxel projection).
- Note that SSA in 2D is a perimeter per unit area (m^-1) or per unit mass (m² kg^-1 with unit depth): units must be stated (outline L799–800 already asks).

**Candidate citations / claims**
- SGCA (p.6–7) is a direct conceptual link between grain contact area and thermal properties — could support one clause in Discussion 5.3 ("contact area, not total surface, is what conduction sees; Flin et al. (2011) proposed grain contact area as a parameter to relate to thermal properties"). Optional; it is a suggestion in their paper, not a result.
