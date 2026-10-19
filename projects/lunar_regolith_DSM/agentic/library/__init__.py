"""The library: one place for everything we have learned from the literature, with a strict split between who may write and who may read.
- Curator (program, with the librarian-curator agent proposing entries): the ONLY writer. Publishes immutable, hashed, versioned files.
- Reference desk (program, with the librarian-reference agent composing answers): the ONLY way other agents get information out. Every request is logged.
- Auditor (program, with the librarian-auditor agent explaining): verifies hashes and consistency; reports to the human.
See research/Library/POLICY.md (copied from agentic/library/POLICY.md) for the rules and their limits."""
