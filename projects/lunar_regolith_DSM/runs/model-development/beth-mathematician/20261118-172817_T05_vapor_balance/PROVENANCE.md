# PROVENANCE: beth-mathematician on T05_vapor_balance (run 20261118-172817_T05_vapor_balance)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T17:28:17-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T05_vapor_balance (task.json sha256 a8bd65485dca, grader sha256 341e5610f7f0)
- git commit: affd1b2f2b; uncommitted definition files: ['?? projects/lunar_regolith_DSM/workflows/model-development/agents/beth-mathematician/tutorial/README.md']
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-172817_T05_vapor_balance/workspace; settings: runs/model-development/beth-mathematician/20261118-172817_T05_vapor_balance/settings.json
- file hashes: see manifest.json

`chain: aee1b885dc47d22d4219942e3361a1fb590f89e6548ecb59159803352ebf554f`

## 2. agent call

- started 2026-11-18T17:28:17-0800, 58.6 s, exit 0, turns 5
- tokens: input 8, output 6143
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 1}
- transcript: transcript.jsonl

`chain: 58fbae3cae006834e904b6402d0a02d256612b517da402679601e05fe9676c9a`

## 3. barrier events

- none

`chain: 99c1afc87fca68e308f20df17157dabe3d3d01ec0ff75a6eee529832032f39e4`

## 4. answer

- status: done
- result_expr: Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + Derivative(rho_i(T(x, y, z, t))*phi(x, y, z, t), t) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivat
- steps: 6; questions: []
- answer.json sha256 aebd35b68752

`chain: 09437e37f2c68fe7dd4694a692081c6c3f332ebc5d7b3d0ed9f3d78b1f7f9d89`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units kg / m ** 3 / s'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform rho_v and T and a static phi (3 draws, 2 points)'}, {'check': 'exchange', 'pass': True, 'details': 'deposition consumes vapor at one ice-like density (draw 2: 1897.2 kg m^-3) for every phi and rho_v: water is conserved'}, {'check': 'diffusion', 'pass': True, 'details': 'a vapor peak decays (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'transport integrates to zero over a period (3 draws): the flux is a divergence'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['Derivative(-phi(x, y, z, t)*rho_v(x, y, z, t) + rho_v(x, y, z, t), t)'], 'diffusion': ['-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), x) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x)', '-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), y) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y)', '-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), z) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)'], 'phase change': ['Derivative(phi(x, y, z, t)*rho_i(T(x, y, z, t)), t)'], 'other': []}, 'current model (src/assembly.c:188-191)': {'storage': ['(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), t)'], 'diffusion': ['-(Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z))'], 'phase change': ['(rho_i - rho_v(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i(T)', 'latex': '\\rho_i(T)', 'meaning': 'density of ice (mass of water per unit volume of ice)', 'units': 'kg m^-3', 'typical': 930}, {'name': 'D_v(T)', 'latex': 'D_v(T)', 'meaning': 'diffusivity of water vapor in the pore gas', 'units': 'm^2 s^-1', 'typical': 2e-05}, {'name': 'm', 'latex': 'm', 'meaning': 'total water mass (ice plus vapor) per unit total volume', 'units': 'kg m^-3', 'typical': 900}, {'name': 'J', 'latex': '\\mathbf{J}', 'meaning': 'vapor mass flux per unit total area', 'units': 'kg m^-2 s^-1', 'typical': 1e-08}, {'name': 'V', 'latex': 'V', 'meaning': 'fixed control volume in the pore space, with boundary \\partial V and outward unit normal \\mathbf{n}', 'units': 'm^3', 'typical': 1e-15}, {'name': 'R', 'latex': 'R', 'meaning': 'residual of the water mass balance', 'units': 'kg m^-3 s^-1', 'typical': 0}], 'agent left out': [{'effect': 'Advection of vapor by bulk gas flow, including Stefan flow caused by sublimation', 'why': 'The Stefan velocity is about $\\rho_v/\\rho_{gas}$ times the diffusive velocity. Near 200 K, $\\rho_v \\sim 10^{-5}$ kg m$^{-3}$, which is small, and no imposed pressure gradient drives flow in the pores.'}, {'effect': 'Thermal diffusion (Soret effect) and pressure diffusion', 'why': 'The Soret coefficient for water vapor in air is of order 0.01-0.1 of the Fick flux per relative temperature change; with pore-scale $\\Delta T/T \\ll 1$ it is a small correction.'}, {'effect': 'Knudsen or transition-regime transport in near-vacuum lunar pores', 'why': 'This is not strictly neglected: it changes the value and the dependence of $D_v$, which can depend on pore size as well as $T$. The form of the balance stays the same, but $D_v$ must be chosen for the gas-pressure regime.'}, {'effect': 'Vapor adsorption on regolith grain surfaces and surface diffusion', 'why': 'The work order defines only ice and gas in the pore space, so grains are outside the domain; adsorbed water would need its own balance.'}, {'effect': 'Ice motion (creep, sintering flow) and density change with temperature', 'why': 'Ice velocities are negligible on the time scales of sublimation. The thermal expansion of ice is about $10^{-4}$ K$^{-1}$; $\\rho_i(T)$ is kept in the result anyway, and it can be set constant.'}, {'effect': 'Anisotropy or tortuosity of diffusion in the diffuse interface beyond the $(1-\\phi)$ weighting', 'why': 'These affect only the interface layer, which is thin compared with the pore size, so the error is of the order of the interface width over the pore size.'}]}}

`chain: dee9c8a0fabb8f3c8db876b136ce21fa9e0361e67e79292a7f1f2af6569dc6d8`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S6: possibly undefined in the equation: \xi

`chain: e948406b701b5ff55cb3d0f416ec35362ef2281009c89b1d492bc499ec8e3570`

## 7. verdict

- **PASS**

`chain: 3fbfc515d84d6ca9f972e8a9a101eda3d003d06715bf0d9da6f303a4212a7a48`
