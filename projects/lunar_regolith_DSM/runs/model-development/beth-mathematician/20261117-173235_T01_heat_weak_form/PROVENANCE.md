# PROVENANCE: beth-mathematician on T01_heat_weak_form (run 20261117-173235_T01_heat_weak_form)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-17T17:32:35-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 5f08282a3bbd)
- task: T01_heat_weak_form (task.json sha256 99e68c5399f3, grader sha256 a1dad86f3681)
- git commit: ec1ffc1f24; uncommitted definition files: ['M projects/lunar_regolith_DSM/scripts/agentkit/launch.py', ' M projects/lunar_regolith_DSM/workflows/model-development/agents/beth-mathematician/env.allow']
- Claude Code: 2.1.326 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261117-173235_T01_heat_weak_form/workspace; settings: runs/model-development/beth-mathematician/20261117-173235_T01_heat_weak_form/settings.json
- file hashes: see manifest.json

`chain: 5161939b8d0c0653f7e0ea8d6983819e2d446e45cd98d3769f1d1179a2cfb146`

## 2. agent call

- started 2026-11-17T17:32:35-0800, 32.6 s, exit 0, turns 6
- tokens: input 10, output 3749
- tool calls: {'Read': 1, 'Write': 2, 'Bash': 1}
- transcript: transcript.jsonl

`chain: 6d329ba41132904ddfec8f4582f17cae5c7893caa33e047144b0ce2d33597d22`

## 3. barrier events

- none

`chain: e20ff1538488e9f5de9c39ba469e2ad4dba91ed1ac176b421f524eecc0c7c8bd`

## 4. answer

- status: done
- result_expr: Integral(rho*c*w(x)*Derivative(T(x, t), t), (x, 0, L)) + Integral(k(x)*Derivative(w(x), x)*Derivative(T(x, t), x), (x, 0, L)) - Integral(w(x)*Q(x), (x, 0, L)) - Subs(w(x)*k(x)*Derivative(T(x, t), x), x, L) + Subs(w(x)*k(x)*Derivative(T(x, t), x), x, 0)
- steps: 7; questions: []
- answer.json sha256 04c984e8fb23

`chain: e4a38091cab1bb7a521a89a64f33f1c708de5a0b61422233e5ac2208f60ebc8e`

## 5. grade

- pass: True
- details: equals 1 times the weak form on 6 random samples

`chain: 8609dd3a7d6c3ce0b84010184925a9f4163f2eb41ed8336c9f290aeaff7c6e62`

## 6. verdict

- **PASS**

`chain: 7757b1d47260a51c5680f911b8e3c562e6446f93dd44797970d9d1a7b4dde203`
