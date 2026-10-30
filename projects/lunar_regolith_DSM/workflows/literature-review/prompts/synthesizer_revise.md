You are the synthesizer (Y1) of a literature review for a phase-field model of
ice sublimation and sintering in a porous sediment. You wrote the draft memo
below; a script (Y2) checked its citations and found the failures listed
after it. Revise the memo so that it passes, using only the checked claims
given. You have no other source; do not use outside knowledge.

Question:
{question}

Your draft:
{memo}

Failures found by the script, each with the sentence as the script read it:
{failures}

Checked claims (JSON; each has an id, the paper, the checker's verdict and,
where the checker narrowed it, the corrected text, which is the one to use):
{claims}

Applicability checklists of the papers (JSON, per paper):
{checklists}

How to fix each failure:
- uncited: if a claim in the list states what the sentence says, cite it as
  [claim_id]. If the sentence is your own inference or proposal, move it to
  Model options or Gaps and disagreements. If no claim supports it and it is
  not needed, delete it. Never cite a claim that does not state the point.
- tag_mismatch or untagged_number: state the number exactly as the cited
  claim gives it (with its precision and qualifiers, such as "of order"),
  followed by [claim_id quantity]; if no claim gives the number, remove it.

Change only what the failures need; keep every other sentence, the section
structure and the citations as they are. The rules of the draft still hold:
every sentence outside Model options, Library silent and Gaps and
disagreements cites at least one claim; every number is followed by a tag
naming the claim it comes from and the quantity; each claim keeps its
conditions. Return the whole revised memo as one Markdown string.
