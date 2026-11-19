# PROVENANCE: beth-mathematician on T05_vapor_balance (run 20261118-123403_T05_vapor_balance)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T12:34:03-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 62aec57faad3)
- task: T05_vapor_balance (task.json sha256 f4aaa549e323, grader sha256 341e5610f7f0)
- git commit: 645fce867a; uncommitted definition files: none
- Claude Code: 2.1.326 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-123403_T05_vapor_balance/workspace; settings: runs/model-development/beth-mathematician/20261118-123403_T05_vapor_balance/settings.json
- file hashes: see manifest.json

`chain: bb7d7eb522ce62413c562a355f5d387db4b8600b12239b0b541fa3f33922971e`

## 2. agent call

- started 2026-11-18T12:34:03-0800, 57.8 s, exit 0, turns 6
- tokens: input 10, output 5668
- tool calls: {'Read': 1, 'Write': 2, 'Bash': 1}
- transcript: transcript.jsonl

`chain: efc73b01b4b86156470e4eaf15cbae477901d2d96a88339d85d28f3c5905d3d2`

## 3. barrier events

- none

`chain: bc4bae5a768fae8c943fed661528c06a78812d73f9af9d60e1fb0ed056039985`

## 4. answer

- status: done
- result_expr: Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + rho_i*Derivative(phi(x, y, z, t), t) - Derivative((1 - phi(x, y, z, t))*D(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative((1 - phi(x, y, z, t))*D(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative((1 - phi(x, y, 
- steps: 6; questions: []
- answer.json sha256 cb13bfc032bf

`chain: 896a8bc45729102af226e11ce13930a4d51acabe48bb6bd324596785ae9d5db7`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units kg / m ** 3 / s'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform rho_v and T and a static phi (3 draws, 2 points)'}, {'check': 'exchange', 'pass': True, 'details': 'deposition consumes vapor at one ice-like density (draw 2: 1840 kg m^-3) for every phi and rho_v: water is conserved'}, {'check': 'diffusion', 'pass': True, 'details': 'a vapor peak decays (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'transport integrates to zero over a period (3 draws): the flux is a divergence'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['Derivative(-phi(x, y, z, t)*rho_v(x, y, z, t) + rho_v(x, y, z, t), t)'], 'diffusion': ['-Derivative(-D(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), x) + D(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x)', '-Derivative(-D(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), y) + D(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y)', '-Derivative(-D(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), z) + D(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)'], 'phase change': ['rho_i*Derivative(phi(x, y, z, t), t)'], 'other': []}, 'current model (src/assembly.c:188-191)': {'storage': ['(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), t)'], 'diffusion': ['-(Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z))'], 'phase change': ['(rho_i - rho_v(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i', 'meaning': 'density of ice, taken as constant', 'units': 'kg m^-3', 'typical': 920}, {'name': 'D(T)', 'meaning': 'diffusivity of water vapor in the pore gas', 'units': 'm^2 s^-1', 'typical': 1e-05}], 'agent left out': [{'effect': 'Thermal expansion of ice (rho_i depending on T)', 'why': 'The volumetric expansion coefficient is about 1.5e-4 K^-1, so rho_i changes by less than 1% over tens of K'}, {'effect': 'Advection by bulk gas flow, including Stefan flow', 'why': 'Near 200 K the vapor mole fraction is about 1e-6 to 1e-5, so the flow that sublimation drives is negligible next to diffusion'}, {'effect': 'Soret (thermal) diffusion', 'why': 'The thermal diffusion factor of water in air is about 0.1 or less, so the term is a small correction to the flux driven by concentration gradients'}, {'effect': 'Adsorption on regolith grains and grain motion', 'why': 'Outside this model of the pore space; the grains are treated as inert boundaries'}, {'effect': 'Knudsen and rarefied-gas transport', 'why': 'Folded into the effective D(T); if it matters, D would also depend on pore size'}]}}

`chain: c23fb74fde1eb2182a127aa54f0277039c896b18155051dcafc4b409ef4390f8`

## 6. verdict

- **PASS**

`chain: 015b6cdced04d0e76551d535b9d71a54b2972d613c926cdcc8b139577a599e38`
