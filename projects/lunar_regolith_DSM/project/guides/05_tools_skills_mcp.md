# Tools, skills, MCP servers and plugins

## Definitions

- **Tool**: a function an agent can call. Built in: Read, Edit, Write, Bash,
  Grep, Glob, WebFetch, WebSearch. Others: your scripts (through Bash) and MCP
  tools. Restrict each agent to the tools it needs.
- **Skill**: a folder `.claude/skills/<name>/` with a `SKILL.md` and any
  scripts or reference files. The name and description are always in context;
  the body loads only when the skill is used; scripts run without their source
  entering context. A skill is a written procedure plus the code it needs.
  Invoke it by task, or as `/<name>`.
- **MCP server**: a program that connects Claude Code to other software
  through the Model Context Protocol (a database, GitHub, a reference
  manager, an instrument). It exposes tools; Claude Code is the client.
- **Plugin**: a package of skills, agents, hooks and MCP servers, installed
  from a marketplace with `/plugin`.

## Skill file format

```markdown
---
name: comp-eps
description: <what it does and when to use it; this line decides when Claude picks it>
allowed-tools: Bash(venv_lunar/bin/python preprocess/comp_eps.py *)
---
<when to use it, the exact commands, what to report, what it may not do>
```

Useful frontmatter: `allowed-tools` (pre-approves exactly these commands),
`disable-model-invocation: true` (only the user can call it), `context: fork`
(runs in a subagent), `paths` (activates only for matching files). Keep
`SKILL.md` under 500 lines; move detail into sibling files.

Skills in this project: `comp-eps` (eps and mesh), `wrap-up` (save state
before closing or clearing).

## MCP

```bash
claude mcp add --transport http <name> <url>              # remote server
claude mcp add --transport stdio <name> -- <command> ...  # local program
claude mcp list
```

Scopes: `local` (you, this project; the default), `project` (committed
`.mcp.json`, shared), `user` (you, all projects). `/mcp` shows status in a
session. Servers run with your permissions: read a server's code before
installing it. Caltech's policy blocks some claude.ai connectors.

Candidates for this project, to evaluate before installing:

| Need | Candidate |
|---|---|
| Paper metadata and DOI checks | a Crossref, OpenAlex or Semantic Scholar server |
| Citation context (supporting/contrasting) | Scite connector (offered on the account) |
| Reference manager | a Zotero server |
| Run database | a SQLite server, or plain scripts |
| GitHub issues and PRs | GitHub's official MCP server, or the `gh` CLI |

A plain script called through a skill is often simpler and more reproducible
than an MCP server. Use MCP when the software is interactive or remote.

## Where to find them

- Skills: https://github.com/anthropics/skills (examples, template, spec).
  Engineering post: https://www.anthropic.com/engineering/equipping-agents-for-the-real-world-with-agent-skills
- Plugins: `/plugin` → Discover. Official marketplace
  `anthropics/claude-plugins-official`; community marketplace
  `anthropics/claude-plugins-community` (not reviewed by Anthropic).
- MCP servers: registry https://modelcontextprotocol.io/registry/about;
  reference servers https://github.com/modelcontextprotocol/servers (examples,
  not production); Anthropic directory https://claude.ai/directory.
- Writing tools well: https://www.anthropic.com/engineering/writing-tools-for-agents
  (few, consolidated tools; unambiguous parameter names; error messages that
  say how to fix the call).
