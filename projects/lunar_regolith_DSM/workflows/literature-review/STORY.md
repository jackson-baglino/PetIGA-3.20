# The literature review, as a story

Version 2 of the design, agreed with the user on 2026-11-03 (planning notes 19 to 22 in `docs/planning/`). This is
the plain-language version, written for you. The same design for agents, with every stage, schema and rule, is
`SPEC.md` in this folder. `WORKFLOW.md` still describes the pilot line as it is built today, until the new design
replaces it (section 9).

Stage codes are written as **code: name**, for example **R3: check facts**. The letter is the phase:
- **F** find;
- **G** get;
- **R** read and check;
- **L** library;
- **M** memo;
- **K** keep;
- **C** calibrate.

The number is the order within the phase.

---

## 1. What this is for

The project needs a library of facts from the literature that it can trust, for three uses:
- choosing the model (Pillar 2);
- checking the implementation (Pillar 3);
- writing the paper.

"Trust" has a narrow meaning here. Every fact points to:
- the paper, page and paragraph it came from;
- the words the paper uses;
- the conditions under which the paper says it holds;
- a check by something other than the thing that extracted it.

The review is one part of the project, not the whole of it. It is built so it can run in the background while
model work goes on, and it never waits on a person except at a few named points.

## 2. Who does what

- **You** set the questions, choose the seed papers, approve each discipline's brief and shortlist, settle the facts
  the checks disagree on, and approve each memo. Nothing else needs you.
- **Scripts** do everything a script can do: search the databases, follow citations, download, check that a file
  is the paper it claims to be, parse, number paragraphs, check quotes, publish, back up. A script's result is the
  same every time.
- **Model calls** make the judgments: write search queries, decide relevance, extract facts, check facts, reread
  equations blind, revise, draft and criticize the memo. Each one is a single request with the inputs chosen by
  the script and an answer that must fit a fixed format (a schema). The model can't open other files, run commands
  or wander off.
- **One agent, the scout**, gets a real tool loop: it searches the web for work the databases miss, decides what to
  look at next, and proposes papers. Scripts check everything it proposes.
- **The stage runner** (`run_b1.py`, renamed when the new design lands) decides what runs next, in which order, and
  when to pause. A rule in code, not an agent, does this.
- **The foreman** is the session you talk to while a review runs. It starts the runner, reads its reports, explains
  them, and gives you the exact commands for your decisions. It can't change anything. A hook stops a real run from
  being started in the main session, so you can't forget to switch.
- **The main session** (the one writing this) designs and builds. It does not run real reviews.

Why so few agents: every model step must name its reason for being separate. There are four reasons:
- **information:** it must not see something, as the checker must not see the reader's reasoning;
- **permissions:** it needs different rights, as the scout needs the web;
- **judgment:** it needs a different kind of judgment;
- **breadth:** the work divides into independent parts.

A step with no reason is a script. A step gets a tool loop only if its next input depends on what it just found.
Every extra agent adds randomness and uses more of your allowance.

## 3. The pilot first: the nine-paper methods review

Before any of the new design is built, the nine-paper methods review finishes on the line we have, as the
**pilot**. It answers the question model development needs now: how to couple the stationary grains with the
evolving ice–vapor interface without a spurious air film. Its results become the **baseline** that the new design
must match or beat on the same papers.

What you do, in order:
1. **Redo the claim you mislabeled** (c006). Why: the library must not contain a claim you now know is wrong. Run
   in your own terminal:
   `venv_workflows/bin/python scripts/literature_review/u3_review.py --run runs/literature-review/b1/dev-20261101-kp --redo kaempfer_plapp_2009-c006`
2. **Rule on the nine flagged control statements.** Why: the checker is judged by how often it flags statements
   that are correct. Until each flag is settled as "the gold statement was wrong" or "the checker was wrong", that
   rate is unknown. All nine involve a dropped condition, so your rulings also test the new condition rule
   (section 5).
3. **Approve the freeze.** Why: from then on, nothing that affects the results can change during the run, so the
   result belongs to one known version of the line.
4. **The nine papers run**, with your gates along the way: copy decisions, the calibration result, the claim
   review, the memo.

## 4. A review in the new design, from start to finish

### Find (F): which papers?

**F1: approve the discipline brief (you).** A review covers one discipline at a time, in this order:
1. numerical methods;
2. planetary observations and setting;
3. experiments (sublimation, ice–sediment, ice sintering and metamorphism);
4. sublimation physics;
5. material properties, looked up when the model needs a property rather than swept.

The brief holds the discipline's questions, its seed papers (3 to 7, at least one review, at least one you know
well), its scope and its inclusion rules. You also list about 20 papers the search **must** find. That is the
recall test: if the search misses them, it isn't good enough yet. Why you: the questions and the scope decide
everything after this point.

**F2: write queries.** One model call turns the brief into search queries.

**F3: search and chain citations.** A script runs the queries against OpenAlex and Semantic Scholar, then follows
the seed papers' references and the papers that cite them, one hop at a time. It stops when a hop adds less than
about 5–10 % new relevant papers. Every query and its results are logged.

**F4: scout the web.** The scout looks for what the databases miss: NASA technical reports, theses, conference
abstracts. It proposes papers with its reasons and adds nothing to the library itself.

**F5: screen relevance.** One model call per candidate paper reads the title and abstract against the brief and
decides: include, exclude with a reason, or belongs to another discipline. **Every candidate gets exactly one
recorded decision**, and a script lists any candidate without one. That record is the coverage ledger.

**F6: approve the shortlist (you).** You see:
- the shortlist;
- a random sample of the rejects;
- the queries and what each one found;
- the recall test's result.

Why you: this is the last cheap point to say "you're missing a whole area" or "this is out of scope".

### Get (G): the papers themselves

**G1: fetch.** A script tries your Zotero library first, then open-access copies, then arXiv only if the paper is
nowhere else. Paywalled papers it can't get go on a list for you to fetch with the Zotero Connector.

**G2: check copy.** A script checks that the file is the paper: the DOI's metadata against the title and authors on
page 1, the retraction status, and the file's hash. A mismatch comes to you.

**G3: prepare paper.** A script renders every page as an image and parses the text, equations (in LaTeX) and
layout with MinerU, on this Mac, offline. It **numbers every paragraph** (`p7.3` is page 7, paragraph 3). Papers are
prepared ahead of reading, for example overnight, so reading never waits for them.

**G4: index paragraphs.** Every paragraph of every collected paper goes into a search index, by keyword and by
meaning. This is **layer 1** of the library: nothing is thrown away, so a question we think of later can still find
text we never extracted. No model is involved.

### Read and check (R): the facts

These steps run for the **core papers** (the seeds and the highest-ranked shortlist papers), and later for any
paragraph a new question needs. Every other paper stays in layer 1 until then.

**R1: extract facts.** One model call reads a group of pages, as images plus numbered paragraphs, and extracts
**every** fact on them, not only the ones our questions ask about. That way there are no blind spots. Each fact has:
- its statement;
- its kind (equation, value, relation, method, observation);
- whose result it is (this paper's or a cited work's);
- a verbatim quote and the paragraph it is in;
- its conditions, each with its own quote and paragraph.

**R2: quote check.** A script checks that each quote is in the paragraph it names, word for word. A quote that
isn't there fails, however plausible it looks.

**R3: check facts.** A separate model call, which never sees the reader's reasoning, first writes a list of every
assumption the whole paper makes. Then, page group by page group, it judges each fact and each condition. Each
condition gets one of three labels:
- **necessary:** the paper says or shows the fact fails without it;
- **scope:** the paper claims the fact only under it and says nothing beyond;
- **unnecessary:** the paper itself shows the fact holds without it.

Necessary and scope conditions both stay with the fact, because the paper's evidence stops there. A fact is
**complete** when every necessary or scope condition the paper attaches to it is listed, including ones stated
elsewhere in the paper (an earlier section's assumptions, a figure caption, a table header). Missing conditions
come back with a quote and a paragraph, so a script can confirm them.

**R4: revise facts.** A fact that fails is revised by one model call and checked again, at most twice. Then it goes
to you.

**R5: blind reread.** For every equation and value, a separate call is told only "equation 13, page 4". It
transcribes what the page prints, and a script compares that with the extracted fact. Agreement there means
something, because the second reader never saw the first one's answer.

### Library (L)

**L1: publish.** A script publishes every fact that passed into the library's **layer 2**: verified facts, each
stored once and never edited. A correction is a new version that points to the old one.

**L2: your fact review (you).** You see each disputed fact with:
- the paper, page and paragraph number;
- the quote and the conditions;
- what the checker said.

You decide right, wrong or unsure. Why you: these are exactly the cases where two independent checks disagree, so
no script can settle them.

**L3: publish your decisions.** A script applies your decisions. A fact you reject is marked rejected and never
deleted.

### Memo (M)

**M1: draft memo.** For each discipline, one call writes a memo answering the brief's questions. Every sentence
cites fact IDs. The memo ends with what the papers support, under which conditions, and what they don't answer.

**M2: citation check.** A script checks every citation, number and unit in the memo against the cited fact.

**M3: skeptic review.** A separate call looks for any statement the cited facts don't support.

**M4: approve memo (you).**

### Keep (K)

**K1: backup.** A script copies the library, with a hash manifest, to the course archive, and a restore test proves
the copy can be read back.

### Calibrate (C)

**C1: calibrate the checker.** Before any real run, the checker is scored on statements with known answers:
- your gold statements, which are correct;
- planted errors, which are wrong on purpose.

It must catch the planted errors and flag at most 10 % of the correct ones. This runs again after any change to a
prompt, a model or the engine.

## 5. What runs where, and the limits

- **Your login, nothing else.** Every model call runs through Claude Code on your enterprise seat, drawing on the
  five-hour usage window. Nothing is billed anywhere else. There is no API key.
- **Usage cap.** Every call reports how much of the five-hour window is used. At about 90 %, the runner starts no
  new work, lets running items finish, and writes a pause report: what ran, where every paper stands, what is
  left, and when the window resets. It then stops, or in "keep running" mode it waits for the reset and continues.
  The last 10 % stays free for your own work.
- **Speed.** Many papers run at once, and papers are prepared ahead of time. Single calls replace the page-by-page
  reading loop. Time per paper is measured in the pilot and again on the new engine.
- **Isolation.** Workers start with no CLAUDE.md, no memory and no settings of ours. This was tested with planted
  marker words on 2026-11-03 and fixed in the same session.
- **Untrusted text.** Paper text is never treated as instructions. A paper that contains planted instructions can
  at most produce a wrong answer, and the checks catch wrong answers.

## 6. What you get

1. **The library.** Layer 1: every paragraph of every collected paper, searchable. Layer 2: verified facts with
   conditions, quotes and locations, each checked by something that didn't write it.
2. **One memo per discipline**, every sentence cited, checked by script and by a skeptic, approved by you.
3. **A record of how each fact was made:** which version of the line, which model, which check, which of your
   decisions.

## 7. When we change the line

During a run, the definitions are frozen. An idea goes to the backlog and is decided before the next run, not in
the middle of one. A real error found by your fact review, the blind reread or the skeptic stops the run. It is
fixed, the checker is calibrated again, and you approve a new freeze.

## 8. How we know it works

- **The pilot baseline:** the new design reruns the nine pilot papers and must match or beat the pilot.
- **A same-effort single-agent comparison** on the first discipline. The team setup has to earn its extra parts.
- **A hand check:** you check 30 facts against the PDFs after the first discipline's pilot of about 40 papers.
- **Calibration** before every frozen run (C1).
- **The recall test** for every search (F1, F6).
- **Repeat runs** of one paper, compared fact by fact, to measure how much the answers vary.

## 9. The order of work

1. **Now to about 8 November: the pilot.** The paused items, the freeze, the nine-paper run, the methods memo.
   Model development (Pillar 2) starts from that memo, for Lab Report 3 on 15 November.
2. **Alongside model work: build the new pieces.**
   - the single-call module and the usage cap;
   - paragraph numbering;
   - the foreman guard;
   - the Find and Get scripts and the index;
   - the new stage codes.

   Then rerun the pilot papers on the new engine and compare.
3. **Then one discipline at a time**, numerical methods first. Each starts with about 40 papers and a hand check
   before it grows.

You decide the dates. This order only makes sure that the review never holds up the model work.
