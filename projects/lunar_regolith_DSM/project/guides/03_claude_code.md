# Using Claude Code well

## Settings files

| File | Applies to | In git | Typical content |
|---|---|---|---|
| managed policy (set by Caltech) | everyone on the account | — | org rules; overrides everything |
| `~/.claude/settings.json` | you, every project, this machine | no | status line, theme, effort level, plugins |
| `.claude/settings.json` | everyone in this project | yes | pinned model, shared deny rules, shared hooks |
| `.claude/settings.local.json` | you, this project, this machine | no | personal allow rules, extra directories, personal hooks |

Precedence: managed, then command-line flags, then local, then project, then
user. Permission rules from all files combine, and a deny anywhere wins.
"Local" means one person on one machine in one project: the file sits in the
project folder but is git-ignored, so it never reaches collaborators or other
clones. Clicking "always allow" at a permission prompt writes a rule into
`settings.local.json`; that is how it fills up with one-off commands.

What researchers usually keep in the local file: allow rules for their own
safe, frequent commands (build, test, read-only git, the project Python),
`additionalDirectories` for data folders outside the repo, and personal hooks
such as desktop notifications. Shared rules (anything that protects the
project) belong in `settings.json`, so they hold for every user and every
clone.

## Working habits

- **Plan first** for anything bigger than a one-sentence change: Shift+Tab to
  plan mode, or `claude --permission-mode plan`. Ctrl+G edits the plan.
- **Give Claude a check it can run** (a test, a script) and ask for its output
  as evidence.
- **Context:** `/context` shows usage; `/compact <focus>` summarizes;
  `/clear` starts fresh (run `wrap-up` first); `/btw` asks a side question
  without adding it to the history.
- **Undo:** Esc Esc or `/rewind`. Commit often; git is the real undo.
- **Shell from inside Claude:** prefix a line with `!`.
- **Headless runs:** `claude -p "<task>" --permission-mode dontAsk
  --output-format stream-json`. Anything not allowed is refused, so this is
  the mode for workflow stages.
- **Inspect what is loaded:** `/memory`, `/hooks`, `/permissions`, `/status`.

## Sessions

- Every session has an ID. The SessionStart hook writes it to
  `.claude/sessions.tsv` with time, source and branch.
- `claude --resume <id>` reopens a session; `claude --continue` reopens the
  latest; `/resume` lists past sessions to pick from. `/rename` gives a
  session a readable name. The VS Code extension also lists past sessions in
  its panel and can reopen them.

## Worktrees for parallel work

A worktree is a second checkout of the same repository, on another branch, in
its own folder. All worktrees share one `.git`, so commits made in one are
visible from the others. This lets two sessions work on two branches at once
without touching each other's files.

```bash
claude -w d0-check        # new worktree .claude/worktrees/d0-check on branch worktree-d0-check
```

- A worktree is a fresh checkout. Untracked and ignored files (venv, `obj/`,
  the built binary, `library/`) are not there; rebuild or point at them.
- When the branch is merged or archived, remove it:
  `git worktree remove .claude/worktrees/<name>`.
- Never make sibling folders by hand (the old `PetIGA-3.20-*`).

## Housekeeping

- Stop background shells and agents when their task ends.
- Close idle background sessions; they keep their context and can still write
  files.
- One task per session where possible; `wrap-up` before closing.

## Anthropic references

- Best practices: https://code.claude.com/docs/en/best-practices
- Memory and CLAUDE.md: https://code.claude.com/docs/en/memory
- Settings and permissions: https://code.claude.com/docs/en/settings, https://code.claude.com/docs/en/permissions
- Hooks: https://code.claude.com/docs/en/hooks-guide
- Checkpointing: https://code.claude.com/docs/en/checkpointing
- Sessions: https://code.claude.com/docs/en/sessions
- Worktrees: https://code.claude.com/docs/en/worktrees
