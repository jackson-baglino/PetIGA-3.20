Run B1 as a bounded review of nine papers, with human approval of the evidence used in one decision memo.
Keep the existing PDF preparation; defer exhaustive extraction, discovery, embeddings, and library infrastructure.
Test conditions, attribution, and equation meaning against known answers—not just counts, agreement, and valid citations.
Use the existing calibration material as a diagnostic; fix the impossible per-kind gate before treating it as certification.
Set a short build budget and move to model work once the memo supports a concrete modeling decision.

1. **Build a small B1 path with a clear stopping point.**

   **Where:** §1–3; STATUS.md, planned order; PILLARS.md.

   **Problem:** The remaining build puts a reusable library system ahead of the first scientific result. Nine known papers do not require automated discovery, screening, hybrid retrieval, or exhaustive structured extraction before they can inform the model.

   **Failure:** The user reaches late November with well-tested infrastructure but no decision about the spurious air film or interface coupling.

   **Fix:**
   
   - **Build now:** a frozen B1 brief and paper list; a minimal resumable driver; one reader and a fresh checker; claim records tied to the exact PDF hash and page; citation validation; one short comparison memo. Reuse the existing Zotero pull, images, tiles, Marker output, and manifests.
   - **Check now:** finish a small T-res fixture, inspect representative layouts, and have the user approve every equation, condition, and result actually used to recommend a model option. Keep the complete PDFs and existing full-page extraction.
   - **Defer:** arXiv/JATS ingestion, whole-paper structured tables and figures, embeddings, the desk, automated acquisition and Zotero upload, general catalog recovery, and broad dimensional-analysis parsing. Read needed tables and equations directly from their images.
   - **Drop from B1:** discovery and screening of the fixed seed set; practice-library rechecking as a prerequisite; routine third and fourth transcriptions of facts the user will inspect; statistical auditing of unused claims.
   - Stop after all nine papers have been considered and the memo identifies supported options, conditions, and unresolved questions. A bounded review should make no field-wide completeness claim.

   **Effect on effort:** Saves large amounts. Estimate **25–40% of the remaining full-plan build**, with low precision because these files contain no implementation estimates. Timebox the small build to roughly 3–5 focused days; reduce automation if it exceeds that budget.

   **Severity:** critical.  
   **Confidence:** high on scope; low on effort estimate.

2. **Transcription agreement cannot verify the meaning of a prose claim.**

   **Where:** X1; §5, “Important claims”; §7 schemas.

   **Problem:** X1 applies to every memo-cited claim and every claim on a model-critical question, but its outputs contain only LaTeX or a quantity, value, and unit. A limitation, attribution, or applicability claim has no defined comparison target.

   **Failure:** “This method handles stationary sediment without an artificial film” receives perfect agreement on its associated equation, although the paper establishes the result only for a different boundary condition.

   **Fix:** Split verification by claim type. Use blind transcription for equations and values; use a semantic check for prose: what was established, by whom, under which assumptions, and with what evidence. For B1, the user should approve these fields for every statement supporting a model option.

   **Effect on effort:** Saves moderate implementation work by narrowing X1; costs a small amount of targeted human reading.

   **Severity:** major.  
   **Confidence:** high.

3. **The checker depends on the reader to identify missing conditions.**

   **Where:** E4; U3a; §5 recall checks.

   **Problem:** E4 receives condition pages the reader identified, plus selected overview pages. A condition omitted by the reader may live elsewhere and never enter the checker’s context. E1/E2 agreement and an audit of extracted claims cannot establish that an unextracted condition was found.

   **Failure:** Both readers report successful convergence but miss a restriction to smooth boundaries or a particular interface-width regime. Every extracted claim passes.

   **Fix:** Before seeing reader claims, have the checker inspect the relevant method section and produce a short applicability checklist: assumptions, boundary conditions, stationary/evolving fields, parameter limits, and validation regime. Test this with known examples whose limiting conditions appear on a different page from the result. Record shared omissions against gold answers, not just differences between readers.

   **Effect on effort:** Neutral to small cost; can replace a full E2 read.

   **Severity:** major.  
   **Confidence:** high.

4. **Several extraction checks compare outputs with the same potentially wrong representation.**

   **Where:** P1; tables; E3; T-table; P2/T-det.

   **Problem:** Two table extractions can agree on shifted headers, and a claim can equal the resulting wrong cell. Deterministic rebuilding proves consistency, not correct chunk boundaries, page attribution, or preservation of context. The layout-block tiler fix also depends on the accuracy of those blocks.

   **Failure:** A value is assigned to the wrong temperature column, a footnote excluding one regime disappears, or an indexed equation is attached to the next paragraph. Cell equality, source agreement, and rebuild hashes all pass.

   **Fix:** Use small image-checked fixtures containing complete headers, footnotes, page mappings, captions, and neighboring text. Include a multi-page table and a figure with labels outside the plot box. Check every memo-used table cell against its image, including headers and footnotes; defer automated extraction of unused tables. Whether Marker supplies sufficiently accurate blocks remains **unverified**.

   **Effect on effort:** Saves moderate work through selective extraction; costs a few fixture checks.

   **Severity:** major.  
   **Confidence:** high on the design weakness; tool behavior unverified.

5. **The proposed LaTeX canonicalizer can erase meaningful distinctions.**

   **Where:** §5, “Comparing LaTeX”; T-res.

   **Problem:** Dropping “similar font commands” is unsafe when fonts distinguish vectors, tensors, sets, or operators. A false-dispute rate measures harmless differences incorrectly rejected; it does not measure meaningful differences incorrectly accepted.

   **Failure:** A scalar and a bold vector canonicalize to the same expression, so a wrong transcription passes all comparisons.

   **Fix:** Preserve semantic font commands by default. Add paired tests with known outcomes: formatting variants that must match, and altered signs, indices, vector notation, derivatives, and grouping that must differ. Report false agreement as well as false disagreement. For B1, unresolved notation differences can go directly to the user.

   **Effect on effort:** Small cost; saves debugging and adjudication later.

   **Severity:** major.  
   **Confidence:** high.

6. **Dimensional consistency is being asked to establish too much.**

   **Where:** E3b; §5 units; T-units.

   **Problem:** Dimensionally consistent equations can still have wrong signs, coefficients, operators, or physical meanings. The description of replacing derivatives with dimensional effects is also too broad: spatial and temporal derivatives differ, and nondimensionalization requires the actual scales, not merely a reader’s label.

   **Failure:** A diffusion term has the wrong sign, or a nondimensional mobility is interpreted using the dimensional units table. Homogeneity passes while the model becomes unstable or physically different.

   **Fix:** For the few equations entering the model, keep a paper-specific symbol and scaling table and have the user verify it. Use simple dimensional checks where supported, plus known-answer checks such as conservation, equilibrium, or a limiting case when applicable. Test the parser on spatial versus temporal derivatives and dimensional versus scaled forms before expanding it.

   **Effect on effort:** Saves large build effort by deferring the general operator layer; costs targeted scientific review.

   **Severity:** major.  
   **Confidence:** high.

7. **Calibration is useful, but the current mixture does not justify a broad reliability claim.**

   **Where:** §5 C0; CALIBRATION_README.md; STATUS.md.

   **Problem:** The 40 gold statements and 68 approved plants are valuable existing assets; writing more is not automatically worthwhile. Mechanical mutations mainly test transcription and validation, while the most consequential failures are omitted conditions, incorrect attribution, and incorrect interpretation. Gold-statement recall also depends on whether the brief actually asks for those facts.

   **Failure:** The system catches many digit and unit changes and passes its pooled gate while missing the condition that makes a proposed method unsuitable. Conversely, it fails recall because it omitted an unrelated gold fact.

   **Fix:** Run the existing material now as a diagnostic. For B1, prioritize its 22 gold statements and the relevant judgment plants; explicitly label which brief question requires each gold fact. Include:
   
   - a condition stated away from the result;
   - a paper discussing another paper’s result;
   - a symbol-sensitive equation;
   - an entirely omitted important statement;
   - unchanged valid controls.

   Measure omission through a fresh reader run. Measure attribution and interpretation through the checker. Compare equations with verified image/source transcriptions. Defer additional B3 and scan calibration unless those inputs are needed for B1; do not commission more plants merely to reach a pooled count.

   **Effect on effort:** Saves moderate user effort and avoids another calibration-writing phase.

   **Severity:** major.  
   **Confidence:** high.

8. **The revised per-kind calibration gate is still mathematically impossible at some sizes.**

   **Where:** §5 C0, held-out pass rules.

   **Problem:** A kind with eight held-out items must meet an 85% one-sided 95% lower bound, but even 8/8 gives only about **68.8%**. No sample size from 8 through 18 can meet that gate, even with zero misses. The 11 user-written plants alone likewise cannot establish 85% sensitivity: 11/11 gives about **76.2%**.

   **Failure:** Adding an eighth example of an error kind converts a reportable result into an impossible gate, stalling the branch.

   **Fix:** Remove the eight-item threshold. Report per-kind counts and intervals without certification at small sizes; require every existing sign/exponent plant to be caught as a fixture rule. If statistical certification remains a later goal, predefine its population and adequate sample size before collecting more items. Do not interpret the aggregate E4 score as demonstrated sensitivity for each judgment kind.

   **Effect on effort:** Saves potentially days of unnecessary plant creation and failed reruns.

   **Severity:** major.  
   **Confidence:** high.

9. **The calibration protocol can exhaust itself before the pipeline stabilizes.**

   **Where:** §5 held-out reuse and reruns; CALIBRATION_README.md split instructions; U3b remedy.

   **Problem:** Fresh held-out sets after failures require a replenishment process and enough independent examples, neither of which is specified. Plants derived from the same gold statement are grouped in the split, which prevents leakage, but they remain correlated observations for a binomial bound. Repeated attempts and combined audit samples also do not automatically retain a nominal 95% guarantee.

   **Failure:** A prompt repair consumes the remaining judgment holdout, and the user must write another answer key before B1 can resume. Multiple mutations of a few facts produce an overly reassuring interval.

   **Fix:** Use tuning fixtures during development and one frozen final evaluation for B1. Report results by source statement and paper as well as by plant; describe intervals as conditional on this constructed test set. After a real failure, repair and inspect the affected claim family rather than automatically rerunning the entire wave and drawing another 60. Any later statistical assurance needs a separately specified sampling and stopping protocol.

   **Effect on effort:** Saves moderate to large recurring effort.

   **Severity:** major.  
   **Confidence:** high.

10. **Merge roles that duplicate judgment; preserve separation where the inputs differ.**

   **Where:** E1–E5, X1, Y1–Y3; §4 agents.

   **Problem:** B1 could accumulate two readers, a two-step checker, a blind reader, math OCR, revisions, two memo-checking roles, and repeatability calls. These calls add orchestration and disputes without evidence here that each improves decisions enough to justify its cost.

   **Failure:** Most user time goes to resolving transcription differences in claims that never affect the model.

   **Fix:**
   
   | Work | B1 implementation |
   |---|---|
   | File identity, hashes, schemas, citations, status, resume | Scripts |
   | Extract question-relevant evidence | One reader call per paper |
   | Inspect assumptions, then verify claims | Fresh checker calls with staged inputs |
   | E2 and X1 | Omit routine E2; use targeted extra reads for unresolved equations |
   | E5 revision | Same reader role in a fresh call; no separate role machinery |
   | Y1 | One synthesis call over approved evidence |
   | Y3 memo checker and skeptic | One fresh skeptical checker covering support and omissions |
   | Final applicability and model choice | Human review |
   | Planter | Stop generating plants for B1 |

   Keep authoring and checking in separate contexts. A script should not decide scientific equivalence or whether an inference is justified.

   **Effect on effort:** Saves large orchestration and runtime effort.

   **Severity:** major.  
   **Confidence:** medium; marginal quality gains from the extra calls have not been measured.

11. **Keep capability isolation, but simplify the main-session boundary.**

   **Where:** §4 security model; reference desk; `wf_launch`.

   **Problem:** Restricting paper readers’ tools, network, and filesystem is proportionate. Requiring the working session to receive only counts and IDs makes evidence discussion and debugging depend on a separate reading interface. The benefit of that additional boundary is unmeasured here.

   **Failure:** The user becomes the manual bridge between the memo, paper, checker, and modeling conversation, or routinely pastes passages back and defeats the intended boundary.

   **Fix:** For B1, use a dedicated review session that may read the selected evidence but has no external-action tools, broad private-file access, or arbitrary shell. Keep acquisition and publishing as separate scripts. Retain read-only sources, controlled output validation, and one writer; defer the counts-only relay and custom desk interface.

   **Effect on effort:** Saves large setup and usability costs.

   **Severity:** major.  
   **Confidence:** medium; whether Claude Code and the macOS sandbox enforce this simplified profile is **unverified** and requires a live refusal test.

12. **The security argument covers exfiltration better than scientific contamination.**

   **Where:** §4 “lethal trifecta”; L1; T-trifecta.

   **Problem:** Removing network access does not prevent a paper from steering a reader toward false, schema-valid claims. Also, “no network” really excludes the model-provider connection: the paper contents are deliberately sent to Claude. Neither schema validity nor the absence of all three trifecta legs establishes scientific integrity.

   **Failure:** A planted instruction produces a plausible false claim with a real quote and a supported verdict. No forbidden tool is called, so the security test passes.

   **Fix:** Extend the existing planted-instruction test to check the resulting scientific claims and publication decision. Use one fixture with a known fact and an instruction to report its opposite; the result must preserve the fact or remain blocked. Treat publisher/checker outputs as untrusted data and keep user approval for model-critical evidence. Check launcher run IDs and paths for containment as well as script-name allow-listing.

   **Effect on effort:** Small cost for a concrete end-to-end fixture.

   **Severity:** major.  
   **Confidence:** high on the gap; actual tool containment is **unverified**.

13. **B1 has no single authoritative roster, tier, or run contract.**

   **Where:** §2; §9 item 38; STATUS.md; U0.

   **Problem:** The request and STATUS.md describe nine B1 papers, §9 lists seven seeds, and STATUS.md reports six acquired seeds. The B1 brief remains open. Tier S permits only five papers and disables E2, while branch calibration requires E1 and E2 together to recover every gold statement.

   **Failure:** The driver cannot determine which papers to run or which acceptance rules apply, or waits for irrelevant discovery and screening stages.

   **Fix:** Before more implementation, write one B1 run contract: all nine paper IDs and copy hashes, acquisition status, core/context designation, questions, required outputs, and explicit stage exceptions. Allow a nine-paper fixed-seed pilot without forcing tier M discovery. Resolve the calibration E2 requirement independently of the production tier.

   **Effect on effort:** Small cost, likely under an hour once the paper list and questions are known.

   **Severity:** major.  
   **Confidence:** high.

14. **Publishing and verification have unresolved ordering and ownership conflicts.**

   **Where:** P1, E1, E3b, L1, Y1–Y2; §6 content layout.

   **Problem:** P1 publishes immutable content, yet that content includes the reader’s symbol list and E3b’s inferred units, which do not exist until later. X1 selection also depends on which claims the memo cites, while Y1 is described as receiving verified claims; the provisional-versus-final status is not explicit. U3a’s newly missed claims are routed through E3/E4 without explicitly inheriting the other applicable checks.

   **Failure:** An immutable folder needs mutation, a memo enters a verification loop, or a human-discovered equation bypasses required checks.

   **Fix:** Keep extraction content immutable and store reader symbols and inferred units as separate versioned analysis records. Define a sequence: checked candidate claims → draft memo → verification of its cited claims → final citation check → approved memo. Put all new or corrected claims through one shared validation function that computes required checks from claim type and intended use.

   **Effect on effort:** Small design cost; saves moderate rework.

   **Severity:** major.  
   **Confidence:** high.

15. **Archive checks do not yet prove recovery or a reproducible public deliverable.**

   **Where:** §6 archive; T-backup; §8; PILLARS.md milestones.

   **Problem:** T-backup tests storage and hashes, but not reconstruction of a usable review from an empty directory. The public repository omits the papers, claims, and memos, so a pinned environment alone cannot reproduce the demonstrated result.

   **Failure:** The local review works, but the course evaluator cannot reproduce its key behavior, or a recovered archive lacks a dependency needed to rebuild the memo’s evidence chain.

   **Fix:** For B1, use a simple snapshot and test restoring one paper, its claims, and its memo citations. Prepare one small redistributable or synthetic fixture with expected outputs for the public repository. Defer the general content-addressed recovery system; preserve the real run privately with a manifest explaining unavailable inputs.

   **Effect on effort:** Saves moderate infrastructure work; costs a small restore test and public fixture.

   **Severity:** major.  
   **Confidence:** high on the missing recovery test; iCloud eviction and recovery behavior is **unverified**.

Questions these files do not answer:

- What are the authoritative nine B1 papers, and which exact copies are already available?
- Which one or two modeling decisions must the first memo support, and what evidence would be enough to start Pillar 2?
- How many hours can the user allocate to pipeline building and evidence review before the conference?
- How many held-out plants actually reach E4, by error kind, source statement, and paper? How many unchanged controls exist?
- Does B1 require any scanned papers, crucial supplements, or corrections?
- Which permission and launcher refusal tests have passed live, beyond the preparation sandbox already reported?
- What public fixture and reproducibility evidence will satisfy the course requirements?