# What lives in .claude/

- `agents/`: one instruction file per agent (name, model tier, tools, description, instructions). Catalogue: `agentic/WORKFLOWS.md`.
- `settings.json`: shared project permissions. It denies the Edit and Read tools on the library folder, so a session reads the library only through the reference desk (`python -m agentic.library.desk`).
- `settings.local.json`: personal permission rules (not shared).
- `ACTIVITY_LOG.md`: running log of sessions.
- Tools are not kept here: Claude Code has no recognised tools folder. See `agentic/tools/REGISTRY.md`.
