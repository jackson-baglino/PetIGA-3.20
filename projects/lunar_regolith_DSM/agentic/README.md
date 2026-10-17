# agentic/: tools, schemas and tests for the agent layer

Everything an agent can call or must obey lives here, versioned with the solver. Agent definitions are in `.claude/agents/`, reusable procedures in `.claude/skills/`; neither contains logic, only instructions that point at the tools below.

```
agentic/
  tools/            one directory (or module) per tool: pure, typed, deterministic, JSON in / JSON out
    mesh_params/    eps, mesh and derived phase-field parameters (was preprocess/comp_eps.py)
    spec/           tool specs (JSON schema + description) exported for Claude / MCP clients
  lit/              literature pipeline scripts (registry, discovery, extraction driver, verification)
  schemas/          record and hand-off schemas (what workers must return)
  tests/            unit tests plus golden fixtures that pin tools to known numbers
  mcp/              (planned) MCP servers that wrap tools, including the guarded solver/HPC tool
```

## Rules for a tool
1. **One job**, named by what it returns. The docstring and `TOOL_SPEC` description are the specification the agent reads (course Lecture 5).
2. **Units in argument names** (`domain_x_m`, `temperature_c`), validated up front; all input errors reported at once.
3. **Pure**: no printing, no file writes inside the library function. File changes happen in a CLI flag or in the caller, so an agent cannot modify inputs by calling the tool.
4. **Structured results**: numbers plus a `checks` list (`ok | note | warn | fail`) instead of prose warnings.
5. **Pinned**: a test freezes the numbers of any tool replacing older code; a drift test compares constants with the C solver.
6. **Stdlib only** for tools that may run on the HPC login node.

## Running
```bash
python -m agentic.tools.mesh_params --Lx 6e-4 --Ly 6e-4 --T0 -10 --alpha 1e-3        # readable summary
python -m agentic.tools.mesh_params --json '{"domain_x_m": 6e-4, "alpha_c": 1e-3}'   # JSON in/out (agents)
python preprocess/comp_eps.py ...                                                     # old command, forwards to the tool
python -m unittest discover -s agentic/tests -v
```
