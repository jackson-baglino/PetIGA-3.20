"""Run a project agent headless for the units-verification workflow, and audit what it did.

call_agent() takes the agent's definition from .claude/agents/<name>.md (model, tools, body as system prompt) and
runs it with scripts/agentlib.py, which copies the given images into a fresh empty folder outside the project (so
no CLAUDE.md, memory or units table is in reach) and runs

    claude -p <task> --model <model> --tools <tools> --permission-mode dontAsk --setting-sources ""
           --strict-mcp-config --no-session-persistence --disable-slash-commands --system-prompt <body>
           --output-format stream-json --verbose --json-schema <schema>

from that folder. Reading outside the folder is refused by the harness in dontAsk mode. The transcript is then
audited: the session model must be the pinned one, and every tool call must be a Read of one of the given images.
The structured output is validated against the schema. A failed call is retried once, then reported.
Each attempt is a row in <run>/runs.db; each transcript is <run>/transcripts/<stage>/<label>_try<k>.jsonl.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import agentlib  # noqa: E402
from agentlib import TIMEOUT_S, AgentError, validate  # noqa: E402,F401

import uvlib  # noqa: E402

AGENTS_DIR = uvlib.ROOT / ".claude" / "agents"


def load_agent(name):
    return agentlib.load_agent(AGENTS_DIR / f"{name}.md", name)


def _audit(events, agent, clean, image_names):
    return agentlib.audit(events, agent, clean, image_names)


def call_agent(rdir, stage, agent_name, label, task, images, schema, check=None, retries=1, allow=()):
    """Run one agent task. images: {name in the clean folder: source path}. check(output) -> extra problems.

    allow: tools to pre-approve with --allowedTools. dontAsk mode refuses any tool needing permission, such as
    WebSearch (tested 2026-10-23); Read inside the clean folder needs none.

    Returns (output, info). Raises AgentError when every attempt failed; info lists the problems.
    """
    return agentlib.call_agent(rdir, stage, load_agent(agent_name), label, task, schema, files=images, check=check,
                               retries=retries, allow=allow)
