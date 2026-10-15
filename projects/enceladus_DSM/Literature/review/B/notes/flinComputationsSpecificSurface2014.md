# flinComputationsSpecificSurface2014 — Flin, Lesaffre, Dufour, Gillibert, Hasan, Rolland du Roscoat, Cabanes & Pugliese, "On the computations of SSA and specific grain contact area from snow 3D images", in *Physics and Chemistry of Ice 2010* (Furukawa et al., eds.), Hokkaido Univ. Press, **2011**, pp. 321–328

**Summary.** Compares three SSA estimators for binary µCT snow images (stereology/intercepts, marching cubes, voxel projection) on ten samples (fresh, rounded, depth hoar, melt forms; 4.9–8.6 µm voxels, coarsened by factors of 2), REV tests, and a new specific grain contact area (SGCA) from curvature-driven grain segmentation.

**Findings relevant to the manuscript**
- p.4 (324): at the finest voxel size "all the methods give similar results within a margin of about ±10 to 15%"; MC "systematically" overestimates; VP sensitive to resolution.
- p.5 (325): "the REV seems to be attained with a sufficient precision before the side of the image cube reaches 2.5 mm."
- p.7 (327): "There is no ideal numerical method to estimate the SSA." Also "the SGCA is a parameter that provides scalar information on the bonding area between grains ... This could be an interesting parameter to relate with thermal or mechanical properties."
- SSA defined per unit mass (p.3 Eq. 1: SSA_tot = surface_tot/mass_tot).

**Claims supported / qualified**
- Outline L802 ("for computing SSA from a field, possibly Flin et al. 2014 (check that the method matches)"): **wrong-citation** — the method does not match. Flin's estimators operate on binarized voxel images. Our SSA comes from the diffuse field: monitoring.c integrates φ²(1−φ)² over the domain and divides by eps, and plot_ssa.py converts with the factor 6, which is exact for the equilibrium logistic profile (∫φ²(1−φ)² dx = eps/6 across one interface). The manuscript should state this formula itself. If an image-based cross-check is wanted, Flin 2011 is the right citation for "binary-image estimators agree within 10–15 %".
- Bib problem: the key says 2014 with no venue. The PDF is the *Physics and Chemistry of Ice 2010* proceedings, published 2011, pp. 321–328.
- Relevant to SSA units: Flin uses per-mass SSA (cm²/g), while our 2-D runs give interface length per area. The manuscript must state which.

**Candidate additions.** The p.7 SGCA sentence points out that SSA alone carries no bond-area information. That fits the manuscript's argument that k_eff is not a single-valued function of SSA, but it is only a suggestion ("could be"), so cite it at most as motivation.
