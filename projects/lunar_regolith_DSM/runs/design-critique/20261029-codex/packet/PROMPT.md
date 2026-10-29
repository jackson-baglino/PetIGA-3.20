You are reviewing the design of a literature-review pipeline before most of it
is built. The files in this folder are everything you get:

- WORKFLOW.md: the design (revision 9). Section 11 lists four earlier
  critique rounds; their findings are applied. Do not repeat them unless the
  fix itself is flawed.
- STATUS.md: what is built, what the first real runs taught us, and what is
  left. Read it first; it corrects WORKFLOW.md where they differ.
- CALIBRATION_README.md: the calibration files the user writes by hand.
- PILLARS.md: the project's goals and schedule.
- RETROSPECTIVE.md: what has gone wrong so far and why.

Read only these files. You have no web access; where a finding depends on how
a tool behaves (Claude Code, Marker, Zotero, macOS sandbox), say so and mark
it "unverified". Do not edit anything.

The goal is a better plan, not a longer one. The user's main work is a
numerical model with a conference in about five weeks, so every hour spent on
this pipeline is an hour not spent on the model. A finding that removes work
is worth as much as one that adds a check. Be specific and skeptical; do not
assume the design is good because it has been reviewed four times.

Answer these, in order:

1. The fastest trustworthy path. What is the smallest version of this
   pipeline that would let the user run the first branch (B1, 9 papers) and
   trust what comes out? List what to build, what to defer and what to drop,
   with the reason for each. Estimate the build effort relative to the full
   plan.
2. Checks that would pass on wrong output. The retrospective shows that
   coverage and "no failed pages" passed while the tiles were unreadable,
   and an equation count was wrong while every test passed. Where else in
   the plan does a check measure something easy instead of what matters?
   For each, name the concrete wrong output that would pass and a check that
   would catch it, preferably one with a known answer.
3. Calibration. Is the calibration (section 5 and CALIBRATION_README.md)
   worth its cost at this scale: 40 gold statements, 11 hand-planted and 57
   script-planted errors on 4 papers? Would it detect the failures that
   matter most (a reader that misses a condition, attributes a result to the
   wrong paper, or misreads an equation), and are the pass rules achievable
   with these sample sizes? What would you change?
4. Agents versus scripts. Which remaining steps should be a script, a single
   model call, or human review instead of what the plan says? Which agent
   roles could be merged?
5. Security model. The design keeps paper text away from the main session
   (prompt injection through papers) and separates reading, network access
   and writing. Is the protection proportionate for a single user reading
   published papers? Name any gap that matters, and any rule whose cost
   clearly exceeds its benefit.
6. Anything else that would make the first branch fail or stall: an ordering
   problem, a missing input from the user, a step with no owner, a format
   two stages disagree on.

Output format. Start with a five-line summary of your main recommendations.
Then the findings, numbered, most important first. For each:
- **Where:** section or stage (e.g. "§5, core calibration" or "P2").
- **Problem:** what is wrong, in two or three sentences.
- **Failure:** a concrete thing that would go wrong.
- **Fix:** what to do instead, specific enough to act on.
- **Effect on effort:** saves, costs or neutral, with a rough size.
- **Severity:** critical, major or minor.
- **Confidence:** high, medium or low, and "unverified" where it applies.

End with the questions you could not answer from these files. Use plain,
concise American English.
