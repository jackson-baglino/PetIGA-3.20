# flinComputationsSpecificSurface2014 — Flin, Lesaffre, Dufour, et al., in Physics and Chemistry of Ice 2010 (Hokkaido Univ. Press, 2011), pp. 321–328
No DOI (manifest confirms). Read in full (8 pp). NOTE: the PDF footer says "Physics and Chemistry of Ice 2010, Edited by Y. Furukawa et al., Hokkaido University Press, 2011" — the key/bib year 2014 should be checked (likely 2011).

## Summary
Compares three SSA estimators on binary µCT images of ten snow samples (fresh, rounded, faceted/depth hoar, melt forms; 4.9–8.6 µm voxels): stereology (ST), marching cubes (MC), voxel projection (VP); tests resolution sensitivity and the SSA REV; defines Specific Grain Contact Area (SGCA) via curvature-driven grain segmentation.

## Findings relevant to the manuscript
- At the finest voxel sizes "all the methods give similar results within a margin of about ±10 to 15%" (p.4/324); MC "systematically provides a higher SSA value" (p.4/324); VP sensitive to resolution (p.4–5).
- SSA REV: "the REV seems to be attained with a cube of 2.5 mm edge for all of the studied samples" (abstract, p.1/321; p.5/325).
- SSA defined as surface/mass (Eq. 1, p.3/323).
- SGCA: "This could be an interesting parameter to relate with thermal or mechanical properties" (p.7/327); aged (low-SSA) snow has large SGCA share (p.6/326).

## Claims supported / contradicted / qualified
- Outline L802 ("for computing SSA from a field, possibly Flin et al. 2014 (check that the method matches)"): the methods are for binary voxel images (ST, MC, VP); the manuscript computes SSA from a continuous phase field. Only a match if the manuscript extracts the φ = 0.5 contour and measures its length (2D analogue of MC/triangulation) — then Flin can be cited for the general approach and its ~10–15% inter-method spread. If SSA is computed as ∫|∇φ| dA (diffuse-interface estimator), Flin is not the right citation. Not a close match as it stands.
- Outline L852–853 (domain 2 mm; RVE 1–2 mm from Calonne 2011): Flin gives 2.5 mm for the SSA REV in 3D snow; the 2 mm cell is slightly below that, though Calonne's conductivity REV is the relevant one. A reviewer may raise it.
- Not cited in the draft.

## Candidate additions
- Possibly for SGCA as a bonding measure (p.7/327) in Section 5.3 discussion of what SSA does not capture (neck area vs. free surface) — loosely relevant; do not add unless the paper discusses contact area.
