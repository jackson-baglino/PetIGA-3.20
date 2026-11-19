# PROVENANCE: beth-mathematician on T05_vapor_balance (run 20261118-162055_T05_vapor_balance)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T16:20:55-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T05_vapor_balance (task.json sha256 a8bd65485dca, grader sha256 341e5610f7f0)
- git commit: 4386cd5971; uncommitted definition files: none
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-162055_T05_vapor_balance/workspace; settings: runs/model-development/beth-mathematician/20261118-162055_T05_vapor_balance/settings.json
- file hashes: see manifest.json

`chain: 55ee2a7d4071db69804a8bd00c5971d39a12665b9a9c17769d0d49f71395fb59`

## 2. agent call

- started 2026-11-18T16:20:55-0800, 65.5 s, exit 0, turns 5
- tokens: input 8, output 6586
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 1}
- transcript: transcript.jsonl

`chain: 2b399052bf4c663382faa4fb11cfef9a06962041574744176dc0da3587072367`

## 3. barrier events

- none

`chain: 383ce4bdecb855a9a0e3ecc2d5e6812d5ed2bff2b90625e5ffa49652444c6ca4`

## 4. answer

- status: done
- result_expr: Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + Derivative(phi(x, y, z, t)*rho_i(T(x, y, z, t)), t) - Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative(D(phi(x, y
- steps: 6; questions: []
- answer.json sha256 3b40cdf82255

`chain: 75c129b35b0cd2bc6bd1bde928b5eaaf1105ea27c3f6e82bf31f68a93d202ff8`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units kg / m ** 3 / s'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform rho_v and T and a static phi (3 draws, 2 points)'}, {'check': 'exchange', 'pass': True, 'details': 'deposition consumes vapor at one ice-like density (draw 2: 1897.2 kg m^-3) for every phi and rho_v: water is conserved'}, {'check': 'diffusion', 'pass': True, 'details': 'a vapor peak decays (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'transport integrates to zero over a period (3 draws): the flux is a divergence'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['Derivative(-phi(x, y, z, t)*rho_v(x, y, z, t) + rho_v(x, y, z, t), t)'], 'diffusion': ['-Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x)', '-Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y)', '-Derivative(D(phi(x, y, z, t), T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)'], 'phase change': ['Derivative(phi(x, y, z, t)*rho_i(T(x, y, z, t)), t)'], 'other': []}, 'current model (src/assembly.c:188-191)': {'storage': ['(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), t)'], 'diffusion': ['-(Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z))'], 'phase change': ['(rho_i - rho_v(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i(T)', 'latex': '\\rho_i(T)', 'meaning': 'density of ice (mass of ice per unit volume of ice)', 'units': 'kg m^-3', 'typical': 930}, {'name': 'D(phi, T)', 'latex': 'D(\\phi, T)', 'meaning': 'effective diffusivity of water vapor through the pore gas, per unit total cross-section (includes the gas fraction, tortuosity and any Knudsen correction; tends to 0 in the ice)', 'units': 'm^2 s^-1', 'typical': 0.0001}, {'name': 'm_w', 'latex': 'm_w', 'meaning': 'total water mass (ice plus vapor) per unit total volume', 'units': 'kg m^-3', 'typical': 1}, {'name': 'j_v', 'latex': '\\mathbf{j}_v', 'meaning': 'mass flux of water vapor per unit total area', 'units': 'kg m^-2 s^-1', 'typical': 1e-06}, {'name': 'Omega', 'latex': '\\Omega', 'meaning': 'arbitrary fixed control volume with boundary $\\partial\\Omega$ and outward unit normal $\\mathbf{n}$', 'units': 'm^3', 'typical': 1e-15}, {'name': 'sigma', 'latex': '\\sigma', 'meaning': 'mass rate of sublimation (ice to vapor) per unit total volume', 'units': 'kg m^-3 s^-1', 'typical': 1e-06}], 'agent left out': [{'effect': 'Advection of vapor by bulk gas flow, including Stefan flow caused by sublimation', 'why': 'Vapor is a trace species or the gas is rarefied; the Stefan velocity is about $\\sigma/\\rho_{gas}$ and gives a Peclet number $vL/D\\ll 1$ at pore scale ($L\\sim10^{-5}$ m).'}, {'effect': 'Thermal diffusion (Soret effect)', 'why': 'The Soret coefficient for water in air is about $10^{-3}$ K$^{-1}$ or smaller, so it is small next to the Fickian flux from the much steeper dependence of saturation density on $T$.'}, {'effect': 'Surface diffusion of adsorbed water on grains and ice', 'why': 'Adsorbed mass on the grains is not part of the model; its flux is small next to gas-phase transport for pore sizes much larger than the adsorbed layer.'}, {'effect': 'Motion and creep of the ice', 'why': 'Ice creep velocities at 200 K are many orders of magnitude slower than sublimation fronts.'}, {'effect': 'Temperature dependence of ice density', 'why': 'Kept in the result through $\\rho_i(T)$; ice thermal expansion is about $1.5\\times10^{-4}$ K$^{-1}$, so setting $\\rho_i$ constant changes the result by under 0.1% per 10 K if a later task wants it.'}, {'effect': 'Water bound in or adsorbed on regolith grains', 'why': 'Outside the ice-gas pore model by the statement.'}]}}

`chain: b98fc377d7e1f5acf9ac7e338497f1a4c1490a8e518c350d56800a206e9ef1d4`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S2: possibly undefined in the equation: n
- format: S3: possibly undefined in the equation: D
- format: S4: possibly undefined in the equation: D
- format: S5: possibly undefined in the equation: D
- format: S6: possibly undefined in the equation: D, R

`chain: 5f8ecf6fbd06e823503259f3871f81112d8959120b8bb92d81cd2bb0143cf24c`

## 7. verdict

- **PASS**

`chain: adde46d3a976a1b5ec02f088b81fd146a7a71116d79c0d2bfa68a69e2fcce339`
