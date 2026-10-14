# Collection check: papers/ and text/

Checked 2026-10-14. This was an independent second pass after the collector's page-1 title check. All 111 files in `text/` were checked against their PDF (mutool) and against the Zotero `.bib` (`DrySnowMetamorphism.bib`). Nothing in `papers/` or `text/` was changed.

## Summary

- **One wrong attachment (blocking).** `libbrechtMeasurementsSurfaceAttachment2013.pdf` is a different paper: arXiv:1302.1231, Libbrecht, *Toward a Comprehensive Model of Snow Crystal Growth Dynamics: 2. Structure Dependent Attachment Kinetics near −5 C*. The cited paper is Libbrecht & Rickerby 2013, *J. Cryst. Growth* 377, 1–8. The correct PDF is already in the library at `Literature/1-s2.0-S0022024813002868-main.pdf` (8 pp, checked). The collector's check missed this because the two titles share many words.
- **Identity.** The other 110 files match the bib's title, first author and year on page 1 or 2. The journal was confirmed for every file with a journal field; for two scans (Auriault 1991, Swinkels 1981) the header OCR is too poor for the automatic match, but the page images agree. Where the PDF prints its own DOI, it matches the bib in every case except Byrne 2024 (known). The DOI is not printed in 25 older PDFs, so it cannot be checked there.
- **Completeness.** No paper is truncated. Every last page ends in references, acknowledgments, an appendix, an erratum notice or a publisher back page. Page-count differences come from cover sheets (IOP, AIP, JSTOR, HAL, HUSCAP), trailing Annual Reviews contents pages, Nature Methods and Extended Data pages, or bib page ranges that are wrong. Two uncited entries are arXiv preprints where the bib cites the journal version (Krishnamachari 1996, Echebarria 2004).
- **The text layer is unreliable for numbers in about 45 files.** The main risks are:
  1. **Minus signs are silently dropped.** In Wiley/AGU and Elsevier PDFs set in "Adv*"/MathematicalPi fonts (about 2002–2014), and in the AIP 1998 file, the text reads `−10 °C` as `10 °C`, `cm² g⁻¹` as `cm2 g1`, and `μm` as `mm`. These were confirmed against page images: Kaempfer 2007, "temperatures of −1.6, −8.3, −19.1, and −54 °C" → `1.6, 8.3, 19.1, and 54 °C`; Chen 2010, "−10 °C cold room" → `10 °C cold room`; Taillandier 2007, "849 cm² g⁻¹ (re = 39 μm)" → `849 cm2 g1 (re = 39 mm)`. In Colbeck 1998, `−` is encoded as `2` ("−3 °C" → `23 °C`) and `=` as `5`.
  2. **The × in powers of ten is dropped, and superscripts are flattened in every file.** Kaempfer & Plapp 2009 Table I gives `3.19×10⁻¹⁰ m` as `3.1910−10 m`, and Piqueux 2009a gives `4.6×10⁻¹⁰ m` as `4.6  1010 m`. Every power of ten in every text file must be checked against the PDF.
  3. **Equation glyphs are substituted.** `=` → `1⁄4`, `+` → `þ`, `( )` → `ð Þ`, `√` → `pffiffi`, decimal point → `:` (Elsevier: `5:6846`), and Greek letters become Latin (`ρ` → `r`, `Γ` → `G`). The ACS papers (Demmenie, Moure, Johnson) render brackets as `j k z {`.
  4. **OCR scans.** Presley 1997 and Colbeck 1983 are badly garbled (fused words, `μm` → `tam`, `β` → `/3`). Flatau 1992 has OCR digit errors in its coefficient tables. Kingery 1960, Wexler 1977 and Swinkels 1981 have OCR equations and numbers (`27 500` → `27 SOO`, `5.2` → `S.2`).
- **Tool caveat (affects every reviewer).** Opening a PDF page with the Read tool's `pages` parameter **fails in this environment**: `pdftoppm is not installed`. Workaround: render the page with mutool, then Read the PNG:
  `mutool draw -q -r 100 -o /path/to/scratch/<key>_p<N>.png Literature/papers/<key>.pdf <N>`
  For a text-only second opinion, use `mutool draw -F txt -o - <pdf> <N>`. It uses the same embedded font maps, so it does not fix dropped minus signs.

## Issues

Severity: **blocking** = a fact-check that uses this file will be wrong. **caution** = usable, but verify the specific item. **cosmetic** = fix in the bib or library for tidiness only.

| citekey | issue | severity | suggested fix |
|---|---|---|---|
| `libbrechtMeasurementsSurfaceAttachment2013` | Wrong paper attached: arXiv:1302.1231 (Libbrecht, "Toward a Comprehensive Model… 2. Structure Dependent Attachment Kinetics near −5 C", 13 pp), not Libbrecht & Rickerby 2013 *J. Cryst. Growth* 377:1–8. No DOI printed. Cited in outline prose. | **blocking** | Use `Literature/1-s2.0-S0022024813002868-main.pdf` (correct paper, 8 pp) for this key. Fix the Zotero attachment. |
| `colbeckSinteringDrySnow1998` | Text layer uses a legacy AIP encoding: `−`→`2`, `=`→`5`, `(`→`~`, `+`→`1`. "Sintered … at −3 °C" reads `23 °C`; "−24 °F" reads `224 °F`; Eq. 1 `Y = R sin(A/2)` reads `Y5R sin~A/2!`. Cited in draft. | **blocking** (numbers, eqs) | Take every number and equation from the PDF page image. |
| `kaempferObservationIsothermalMetamorphism2007` | Minus signs and μ dropped: all four test temperatures (−1.6, −8.3, −19.1, −54 °C) read as positive; "10 μm pixel" reads `10 mm`. Checked against page 1. Cited in outline. | **blocking** (numbers) | Take temperatures, sizes and exponents from the PDF. |
| `taillandierRateDecreaseSpecific2007` | Minus signs and μ dropped (`cm2 g1`, `re = 39 mm` for 39 μm, 42 unsigned °C), decimal points become `:`. Checked against the last page. Cited in outline. | **blocking** (numbers) | Take numbers from the PDF. |
| `chenStructuralEvolutionIcesphere2010` | Minus signs dropped: "−10 °C cold room", "−20 and −10 °C", "−180 °C" read as positive. Checked against page 2. Cited in draft. | **blocking** (numbers) | Take temperatures from the PDF. |
| `kaempferPhasefieldModelingDry2009` | Table I (parameter values) is scrambled: × is dropped and mantissa and exponent are fused (`3.1910−10 m`, `1.3810−23 J K−1`), the Greek symbol column is lost, and rows are out of order. Cited in draft and outline. | **blocking** (table values) | Read Table I and the equations from the PDF. |
| `presleyEffectBulkDensity1997` | Poor OCR: spaces lost (`theverylongtimeregime…`), `μm`→`tam`, `≈`→`•`, `−`→`ø`; 42 % dictionary-word rate. Not cited. | caution (unusable text) | Use page images only. |
| `colbeckTheoryMetamorphismDry1983` | OCR: fused words, `β`→`/3`, figure-axis debris mixed into the text; 77 % dictionary-word rate. Cited in draft. | caution | Read equations and numbers from the PDF. |
| `flatauPolynomialFitsSaturation1992` | OCR of the coefficient tables has digit errors (`0.I21487375E-04`, `0.17825 5421E-09`, `0.64708705 IE-11`); `a1`/`a7` read as `a,`; the text layer is duplicated. JSTOR cover page adds 1 page (8 vs 1507–1513). Cited in draft. | caution | Take the polynomial coefficients from the PDF. |
| `wexlerVaporPressureFormulation1977` | OCR: equations garbled, "81A" reads `B1A`, `the` reads `Ihe`. Bib `pages = {5}`; should be 5–20. Cited in draft. | caution | Take equations and coefficients from the PDF. Fix bib pages. |
| `kingery_regelation_1960` | OCR: `27 500` reads `27 SOO`, `5.2 kcal` reads `S.2 kcal`, `γ`→`y`, exponents lost in Eq. 7. AIP cover page adds 1 page, and page 2 opens with the end of the previous article (journal p. 833). Cited in outline prose. | caution | Take Eqs. 7–8 and numbers from the PDF. |
| `swinkelsSecondReportSintering1981` | Scan with OCR: Table 2 (rate equations) is unreadable in text; the year is not visible on the OCR'd page 1. Cited in draft and outline. | caution | Read Table 2 and the equations from the PDF. |
| `molaroMicrostructuralEvolutionWater2019` | Equation glyphs substituted (`=`→`1⁄4` ×53, `√`→`pffiffi`, `+`→`þ`); some exponent minus signs lost. Bib pages 243–277, but page 1 prints "124, 243–276" (34 pp, complete, ends with erratum note). PDF is the corrected version (10 Aug 2020). Cited in draft and outline. | caution | Read equations and Table 5 values from the PDF. Bib pages → 243–276. |
| `calonneNumericalExperimentalInvestigations2011` | Equation glyphs substituted (`=`→`1⁄4`, `+`→`þ`, `( )`→`ð Þ`). Bib `pages = {n/a--n/a}`; the article number is L23501. Cited in draft and outline. | caution | Read regression formulas from the PDF. Bib pages → L23501. |
| `legagneuxMeanFieldModel2005` | Minus signs dropped ("isothermal experiment at 2°C" = −2 °C; `cm2 g1`; `J mol1 K1`); eqs use `1⁄4`, `ffiffi`; 69 % of lines are 1–3-character equation fragments. Cited in outline. | caution | Take numbers and equations from the PDF. |
| `schneebeliTomographyTemperatureGradient2004` | Minus signs dropped (`about 8°C`, `K m1`, `m1`). Cited in outline. | caution | Take numbers from the PDF. |
| `wangEvolutionSpecificSurface2014` | Minus signs dropped ("mean temperature of 4°C" = −4 °C, `K m1` ×69). Not cited. | caution | Take numbers from the PDF. |
| `loweInterfacialStructuralRelaxations2011` | Mostly clean, but `( )` brackets around correlation-function averages are lost in Eq. 3 and some equations are split from their numbers. Cited in draft and outline. | caution | Read equations from the PDF. |
| `demmenieIceSinteringSublimation2025`, `mourePhaseFieldModelWet2024`, `johnsonMicrostructureIcyLunar2025` | ACS typesetting: brackets and operators become `j k z {` (Demmenie Eq. 3 reads `u i2 yn = Cn(T)j z u tk{`). Moure has `T0 = 273.15 °C`, probably −273.15. Demmenie and Moure are cited in draft and outline. | caution | Read equations from the PDF. |
| `grottThermalConductivityMartian2021` | Symbol-font private-use characters: `±`, `=`, Greek letters lost (`1211149 kg m−3`, `cp`). Not cited. | caution | Take values and uncertainties from the PDF. |
| `jansenGeneralizedaMethodIntegrating2000` | Symbol font mapped to Latin-1: `=`→`\x88`, `(`→`\x85`, Greek→Latin; equations unreadable. Not cited. | caution (unusable eqs) | Read from the PDF. |
| `krolUpscalingIceCrystal2018`, `kempfHowEnceladusDust2010`, `libbrechtGrowthRatesPrincipal2003`, `lehningPhysicalSNOWPACKModel2002`, `millerMicrostructuralApproachPredict2003`, `liMetamorphismObservationModel2022`, `piqueuxModelThermalConductivity2009`, `piqueuxModelThermalConductivity2009a`, `hamiltonObservationsPreliminaryScience2014`, `spencerEnceladusActiveCryovolcanic2009`, `ramirez-cabanasExploringGeneralChemistry2024`, `shimizuImpactComplexTopology2017`, `srivastavaObservationTemperatureGradient2010`, `matzlMeasuringSpecificSurface2006`, `flinThreedimensionalGeometricMeasurements2004`, `taillandier…` (above) | Same Adv-font family: `=`→`1⁄4`, `+`→`þ`, `( )`→`ð Þ`, `√`→`ffiffi`; exponent minus and × dropped (`m s1`, `kg s1`, `4.6  1010 m`, `1.38  1023`); decimal→colon (Kempf `5:6846`, Spencer 2009 ×25). Cited: Kempf, Spencer 2009, Matzl (outline); Lehning (outline prose); Flin 2004 (draft). | caution | Take any number with a sign, exponent or unit from the PDF. |
| `bouvetSnowEquiTemperatureMetamorphism2022`, `byrneRadarAttenuationSubsurface2024`, `ghanbarianThermalConductivityPorous2016` | 41–81 % of lines are 1–3-character fragments, so equations and tables are broken up. Bouvet 2022 is cited in draft and outline; Byrne in outline. | caution | Read equations and tables from the PDF. |
| `bretinPhasefieldModellingComputing2019`, `fourteauImpactWaterVapor2021`, `bouvetMultiscaleModelingHeat2024`, `vionnetDetailedSnowpackScheme2012`, `ettrichModellingTransientHeat2014`, `hemingwayInteriorEnceladus2018` | Large braces and matrix brackets in private-use characters (``…``). Equation content is otherwise readable. | cosmetic | Check systems of equations and matrices in the PDF. |
| `sommerfeldClassificationSnowMetamorphism1970` | Moderate OCR (`J fllIl"1Ja! (.f G/aciolvg)`; 89 % dictionary-word rate). 15 pp = journal pp. 3–17, complete; references end on p. 17. Bib says 3–18. | cosmetic | Text is fine for classification wording. Check the bib page range. |
| `byrneRadarAttenuationSubsurface2024` | PDF is the published JGR Planets article, doi:10.1029/2024JE008615; the bib has the ESS Open Archive preprint DOI and type `misc`. (Already known.) | caution (bib) | Update the bib DOI, journal and type. |
| `flinComputationsSpecificSurface2014` | Bib has year 2014 and no venue. The PDF is *Physics and Chemistry of Ice 2010*, eds. Furukawa et al., Hokkaido University Press, **2011**, pp. 321–328; Calonne 2011's reference list cites it the same way. Cited in outline. | caution (bib) | Fix the bib: `@inproceedings`, year 2011, booktitle, pages 321–328. |
| `krishnamachariGibbsThomsonFormulaSmall1996` | arXiv/REVTeX preprint (27 pp, regenerated 2018), with a different title wording ("Gibbs-Thomson formula at small island sizes - corrections for high vapour densities"); the bib cites *Phys. Rev. B* 54, 8899–8907. Not cited. | caution | If cited, check values against the published version or note that the preprint was used. |
| `echebarriaQuantitativePhasefieldModel2004` | arXiv preprint (51 pp); the bib cites *Phys. Rev. E* 70, 061604. Not cited. | caution | Same as above. |
| `libbrechtComprehensiveModelSnow2014` | arXiv preprint; the bib also cites arXiv (10.48550/arXiv.1407.0740). Consistent. | none | — |
| `auriaultHeterogeneousMediumEquivalent1991` | Acrobat Capture scan: journal header OCR'd as `ht. J. Engng Sci.`; body text OK (96 % dictionary-word rate). The DOI is not printed. | cosmetic | — |
| `bullardNumericalSimulationsTransientstage1997` | Scan with jumbled page-1 OCR (title split across lines, pages `12% 139`); body OK. | cosmetic | — |
| `noauthor_sintering_nodate` (Hobbs & Mason 1964) | HTML paywall stub only, so there is no full text. Cited in outline prose. (Already known.) | caution | Mark claims that rest on it `not-checkable`. |
| `libbrechtPhysicalDynamicsIce2017` | 27 pp vs 271–295 (25): 2 trailing Annual Reviews contents pages. The printed DOI wraps across a line (`10.1146/annurev-`). | cosmetic | — |
| `spencerEnceladusActiveIce2013` | 27 pp vs 693–717 (25): trailing contents pages. | cosmetic | — |
| `blackfordSinteringMicrostructureIce2007`, `flinFullThreedimensionalModelling2003`, `ettrichModellingTransientHeat2014` | IOP "To cite this article" cover page (+1 p). | cosmetic | PDF page N = article page N−1. |
| `calonneThermalConductivitySnow2019` | HAL cover page (+1 p); the following pages are the published GRL version. | cosmetic | — |
| `colbeckSinteringDrySnow1998`, `kingery_regelation_1960` | AIP cover page (+1 p). | cosmetic | — |
| `flatauPolynomialFitsSaturation1992`, `hemingwayInteriorEnceladus2018` | JSTOR cover page (+1 p); Hemingway also has a blank final page (p. 23). | cosmetic | — |
| `stehleMigrationBubblesIce1967` | HUSCAP cover page. Bib has no venue: *Physics of Snow and Ice: Proceedings* 1(1), 219–232 (1967). | cosmetic (bib) | Add the venue and pages to the bib. |
| `postberg_macromolecular_2018` | 26 pp vs 564–568 (5): includes Methods and Extended Data. Pages 11, 15 and 17 are figure-only with no text. | cosmetic | — |
| `prokschDensitySpecificSurface2015` | 18 pp vs 346–362 (17): an erratum note on the last page. | cosmetic | — |
| Wiley-downloaded files (Presley, Proksch, Molaro, Grott, McDowell 2023, Li, Wang 2017, etc.) | A "Downloaded from … Wiley Online Library" banner on every page is mixed into the text. | cosmetic | Ignore the banner lines. |
| Supplements | None attached (Choukroun 2020, Molaro 2019, Byrne 2024, McDowell 2023 and Demmenie 2025 all mention SI). | caution | Claims that rest on SI tables are `not-checkable`. |

## Papers whose equations or numbers must be read from the PDF

Use the mutool render workaround above. **Applies to all files:** superscripts and subscripts are flattened (`10⁴` → `104`, `m²` → `m2`), so check any power of ten, exponent or unit in the PDF.

**Text unusable or numbers silently wrong (always use the PDF for numbers, signs, units and equations):**
- `colbeckSinteringDrySnow1998`: draft
- `chenStructuralEvolutionIcesphere2010`: draft
- `kaempferPhasefieldModelingDry2009`: draft and outline (Table I, equations)
- `kaempferObservationIsothermalMetamorphism2007`: outline
- `taillandierRateDecreaseSpecific2007`: outline
- `legagneuxMeanFieldModel2005`: outline
- `schneebeliTomographyTemperatureGradient2004`: outline
- `colbeckTheoryMetamorphismDry1983`: draft
- `swinkelsSecondReportSintering1981`: draft and outline (Table 2)
- `flatauPolynomialFitsSaturation1992`: draft (coefficient tables)
- `wexlerVaporPressureFormulation1977`: draft
- `kingery_regelation_1960`: outline prose
- `presleyEffectBulkDensity1997`
- `jansenGeneralizedaMethodIntegrating2000`
- `grottThermalConductivityMartian2021`
- `wangEvolutionSpecificSurface2014`

**Equations garbled, prose and plain numbers mostly fine (use the PDF for equations, signed values and exponents):**
- `molaroMicrostructuralEvolutionWater2019`: draft and outline
- `calonneNumericalExperimentalInvestigations2011`: draft and outline
- `loweInterfacialStructuralRelaxations2011`: draft and outline
- `demmenieIceSinteringSublimation2025`: draft and outline
- `mourePhaseFieldModelWet2024`: draft and outline
- `bouvetSnowEquiTemperatureMetamorphism2022`: draft and outline
- `flinThreedimensionalGeometricMeasurements2004`: draft
- `matzlMeasuringSpecificSurface2006`, `kempfHowEnceladusDust2010`, `spencerEnceladusActiveCryovolcanic2009`, `byrneRadarAttenuationSubsurface2024`: outline
- `lehningPhysicalSNOWPACKModel2002`: outline prose
- `krolUpscalingIceCrystal2018`, `libbrechtGrowthRatesPrincipal2003`, `millerMicrostructuralApproachPredict2003`, `liMetamorphismObservationModel2022`, `piqueuxModelThermalConductivity2009`, `piqueuxModelThermalConductivity2009a`, `hamiltonObservationsPreliminaryScience2014`, `ramirez-cabanasExploringGeneralChemistry2024`, `shimizuImpactComplexTopology2017`, `srivastavaObservationTemperatureGradient2010`, `johnsonMicrostructureIcyLunar2025`, `ghanbarianThermalConductivityPorous2016`

**Braces and matrices only (cosmetic):** `bretinPhasefieldModellingComputing2019`, `fourteauImpactWaterVapor2021`, `bouvetMultiscaleModelingHeat2024`, `vionnetDetailedSnowpackScheme2012`, `ettrichModellingTransientHeat2014`, `hemingwayInteriorEnceladus2018`.

The other files have clean text layers: Copernicus/TC, recent AGU STIX, Elsevier CharisSIL/STIX, LaTeX/REVTeX, Nature, PNAS and Frontiers. Examples are Braun 2024, Calonne 2014, Fourteau 2021b, Granger 2021, Nicoli 2011 and Choukroun 2020 (a few `=`→`1⁄4`). Prose, plain numbers and equation structure in these files can be taken from the text, but the general rule on superscripts still applies.

## Method

- **Identity:** the normalized bib title (full, or the part before the colon, or its first 40 characters), the first author's surname and the year were searched for in PDF pages 1–3. The bib journal words were searched for on the same pages. DOIs were extracted from pages 1–3 and compared with the bib and manifest DOIs. The manifest's DOI, year and page count agree with the bib and with mutool for all entries.
- **Completeness:** the mutool page count was compared with the bib `pages` range. The last page of each text file was checked for references, acknowledgments or an appendix. Printed page-number sequences were checked for Molaro and Sommerfeld.
- **Text quality:** for each file, the following were measured: private-use and replacement characters, the dictionary-word rate (`/usr/share/dict/words`), fused words, the fraction of 1–3-character lines, glyph-substitution markers (`1⁄4`, `þ`, `ð Þ`, `ffiffi`, decimal-colon), dropped minus signs in units and °C values, and PDF producer and fonts (OCR software and Adv*/MathematicalPi fonts).
- **Spot check (15 papers):** page 1 and/or the last page were rendered with mutool and compared with the text: Abramov 2009, Calonne 2011, Chen 2010 (p. 2), Colbeck 1998, Demmenie 2025, Flatau 1992, Hemingway 2018, Jansen 2000, Kaempfer 2007, Kingery 1960, Libbrecht 2013, Molaro 2019, Presley 1997, Sommerfeld 1970, Taillandier 2007 and Wexler 1977. Identity and page ends agreed in all of them. The text-layer defects listed above (signs, μ, OCR) were confirmed in Kaempfer 2007, Chen 2010, Taillandier 2007, Presley 1997 and Kingery 1960. The wrong Libbrecht 2013 attachment was confirmed from its page 1 (arXiv stamp 1302.1231v1).
