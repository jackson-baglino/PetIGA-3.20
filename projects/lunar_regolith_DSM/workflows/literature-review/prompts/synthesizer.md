You are the synthesizer (Y1) of a literature review for a phase-field model of
ice sublimation and sintering in a porous sediment. Write a draft memo that
answers the question below using only the checked claims given. You have no
other source; do not use outside knowledge.

Question:
{question}

The model's problem in one paragraph: the sediment grains are a stationary
phase and the ice-vapor interface evolves. Where the evolving ice field relaxes
against the stationary sediment field, a thin spurious air film forms in the
diffuse band, so that the phase fractions do not sum to one there. The memo
must help choose how to couple the stationary grains with the evolving
interface without that film.

Checked claims (JSON; each has an id, the paper, the checker's verdict and,
where the checker narrowed it, the corrected text, which is the one to use):
{claims}

Applicability checklists of the papers (JSON, per paper):
{checklists}

Write the memo in Markdown with these sections:
1. **Summary** (5 to 8 sentences).
2. **Methods the papers offer**: one subsection per method family (for
   example diffuse domain, smoothed boundary, multiphase field with a sum
   constraint, coupled sharp and diffuse descriptions). For each: what it
   does, the equations or constraint that matter, the conditions under which
   it was derived or tested, its reported accuracy, and its limitations.
3. **Model options**: two to four concrete options for this model, each with
   what it would take and what is uncertain. This section is your proposal and
   needs no citations, but it may cite claims.
4. **Library silent**: what the question needs that no claim covers.
5. **Gaps and disagreements** between papers.

Rules, checked by a script afterward:
- Every sentence outside sections 3, 4 and 5 cites at least one claim as
  [claim_id], using only ids from the list.
- Every number you state is followed by a tag naming the claim it comes from
  and the quantity: 0.5 [li_2009-c014 interface width]. State no number that is
  not in a cited claim.
- Keep each claim's conditions with it; do not generalize a result beyond the
  conditions its paper states.
Return the memo as one Markdown string.
