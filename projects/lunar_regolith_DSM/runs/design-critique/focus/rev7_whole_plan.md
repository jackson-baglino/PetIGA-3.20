This is revision 7 plus later edits. A previous critique round found 25
findings; they are applied. Do not repeat them unless the fix itself is
flawed. In addition to your usual review, answer these:

1. Weakest points of the whole plan. Name the three places where the design
   is most likely to fail in practice (produce wrong records that pass every
   check, stall, or cost far more than it returns), and why.
2. Review mode and the lethal trifecta (section 4). Is the claim that the
   main session holds no way out and no untrusted text during a wave
   actually enforced by the settings described? Name any path that the
   design misses: a tool, a script output, a file the session can read, a
   host on the allow-list, or the user's own actions.
3. The two-layer library (section 6, edit of 2026-10-28): each paper's whole
   extracted content is stored, never edited, and only checked claims are
   cited. Where could unchecked content leak into a memo or a decision as if
   it were checked? Is anything a later question would need still missing
   from what is stored?
4. Calibration (section 5). Are the pass rules, the held-out split, the
   Clopper-Pearson bounds and the sample sizes sound, and do they measure
   what the acceptance test needs? Is anything measured on a population
   other than the one it is applied to?
5. Simplicity. Which stages, scripts or rules could be removed without
   losing a check the acceptance test needs?
