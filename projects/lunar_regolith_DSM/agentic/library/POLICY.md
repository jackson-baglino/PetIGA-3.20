# Library policy

**Purpose.** One place for everything we learned from the literature (papers, their extracted and checked records, search sets, memos, audits, errata), kept safe from accidental or careless change.

**Who may do what**
- **Curator** (the program `agentic.library.store.Curator`, with the `librarian-curator` agent proposing catalog entries): the only writer. It publishes files that are hashed, versioned (old versions go to `history/`), made read-only and logged as events. Entries are published only after the deterministic checks of the intake pipeline (quote matching, independent claim check, page-image check) have run.
- **Reference desk** (`agentic.library.desk`, with the `librarian-reference` agent composing answers): the only way any other agent gets information out. Every request needs a stated purpose, is checked against a policy table (who may ask for what), is answered from hashed files (a file that no longer matches its hash is refused) and is logged.
- **Auditor** (`agentic.library.verify`, with the `librarian-auditor` agent explaining the result): re-hashes everything and reports problems to the human. It cannot change anything.
- **Everyone else** (extractors, checkers, synthesizers, reviewers, drafters, deriving agents): read through the desk, never write. Corrections arrive as new versions, errata or status changes made by the curator, with a reason, and appear in the event log.
- **The human** may ask for anything and decides on status changes such as "retracted".

**How it is enforced, and where it is not**
1. The headless workers have no file tools at all; they only see text the driver hands them, so they cannot touch the library.
2. Project permission rules deny the Edit/Write and Read tools on the library folder for Claude Code sessions (including subagents), so a session must go through the desk command.
3. Every file and directory is read-only at the operating-system level, and the catalog is opened read-only by readers.
4. Hashes and an event log make any change detectable (`verify`).
These are strong against mistakes and against an agent following its instructions. They are **not** a wall against a deliberate shell command run as the same user (for example a script that changes file modes first): the operating system cannot tell the curator from another process of the same user. The optional Claude Code sandbox (`denyWrite`/`denyRead`) closes most of that gap for shell commands. Detection (hashes, events) is the backstop.

**Layout**: `catalog.sqlite` (papers, claims with full-text search, flags, candidates, searches, memos, events, access log); `papers/<id>/` (record.json, record.md, verification.json, text.txt, variants/); `memos/`; `searches/` (query sets, wave manifests); `history/` (superseded versions); PDFs stay in the Literature folders and are recorded by path and hash.
