# PROVENANCE: beth-mathematician on T05_vapor_balance (run 20261118-173247_T05_vapor_balance_baseline)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T17:32:47-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T05_vapor_balance (task.json sha256 a8bd65485dca, grader sha256 341e5610f7f0)
- git commit: affd1b2f2b; uncommitted definition files: ['?? projects/lunar_regolith_DSM/workflows/model-development/agents/beth-mathematician/tutorial/README.md']
- Claude Code: 2.1.327 (Claude Code)
- MODE: single-prompt baseline (--no-tools): no tools, work order in the prompt; tools: none; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-173247_T05_vapor_balance_baseline/workspace; settings: runs/model-development/beth-mathematician/20261118-173247_T05_vapor_balance_baseline/settings.json
- file hashes: see manifest.json

`chain: 6645b2b90c982f5478c4e1b6c58e55d5778342f09e15205c69c9f9808fd962fa`

## 2. agent call

- started 2026-11-18T17:32:47-0800, 45.9 s, exit 0, turns 2
- tokens: input 2, output 5148
- tool calls: none
- transcript: transcript.jsonl

`chain: 594e253da02b0094ee4f6382284125a9e67b40c5d9ad05dbc5fa8b174db7418b`

## 3. barrier events

- none

`chain: c88f56b3b8b95c88d2d9b48e7b738f5ab199a5cf7b0571e803bdaee51216502f`

## 4. answer

- status: done
- result_expr: Derivative((1 - phi(x, y, z, t))*rho_v(x, y, z, t), t) + rho_i*Derivative(phi(x, y, z, t), t) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative((1 - phi(x, y, z, t))*D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative((1 - phi(x,
- steps: 5; questions: ['Should the vapor flux carry the gas-fraction weighting $(1-\\phi)$, as used here so that water mass is conserved exactly? The alternative is the common simplified form $\\partial_t\\rho_v = \\nabla\\cdot(D_v\\nabla\\rho_v) - \\rho_i\\partial_t\\phi$.']
- answer.json sha256 c6d9e80a04fe

`chain: f6c1078ceca585c28de2ccd968f7aefcf958a386b926b6e9f75e11cf595d261a`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units kg / m ** 3 / s'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform rho_v and T and a static phi (3 draws, 2 points)'}, {'check': 'exchange', 'pass': True, 'details': 'deposition consumes vapor at one ice-like density (draw 2: 1834 kg m^-3) for every phi and rho_v: water is conserved'}, {'check': 'diffusion', 'pass': True, 'details': 'a vapor peak decays (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'transport integrates to zero over a period (3 draws): the flux is a divergence'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['Derivative(-phi(x, y, z, t)*rho_v(x, y, z, t) + rho_v(x, y, z, t), t)'], 'diffusion': ['-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), x) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x)', '-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), y) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y)', '-Derivative(-D_v(T(x, y, z, t))*phi(x, y, z, t)*Derivative(rho_v(x, y, z, t), z) + D_v(T(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z)'], 'phase change': ['rho_i*Derivative(phi(x, y, z, t), t)'], 'other': []}, 'current model (src/assembly.c:188-191)': {'storage': ['(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), t)'], 'diffusion': ['-(Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) + Derivative(D(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z))'], 'phase change': ['(rho_i - rho_v(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i', 'latex': '\\rho_i', 'meaning': 'density of ice (mass of ice per unit volume of ice), taken constant', 'units': 'kg m^-3', 'typical': 917}, {'name': 'D_v(T)', 'latex': 'D_v(T)', 'meaning': 'diffusivity of water vapor in the pore gas', 'units': 'm^2 s^-1', 'typical': 2e-05}, {'name': 'm_w', 'latex': 'm_w', 'meaning': 'water mass (ice plus vapor) per unit total volume', 'units': 'kg m^-3', 'typical': 917}, {'name': 'j_v', 'latex': '\\mathbf{j}_v', 'meaning': 'vapor mass flux per unit total area', 'units': 'kg m^-2 s^-1', 'typical': 1e-07}, {'name': 'Omega', 'latex': '\\Omega', 'meaning': 'fixed control volume in the pore space, with boundary \\partial\\Omega and outward normal \\mathbf{n}', 'units': 'm^3', 'typical': 1e-12}], 'agent left out': [{'effect': 'Advection of vapor by bulk gas flow (including the Stefan flow at the interface)', 'why': 'The vapor is very dilute near 200 K (vapor pressure about 0.16 Pa, so $\\rho_v$ is about 1e-6 kg m^-3), so the Stefan velocity is about $D_v\\nabla\\rho_v/\\rho_{gas}$ and is negligible compared with diffusion. There is no imposed pressure gradient across the pores.'}, {'effect': 'Thermal diffusion (the Soret effect) and cross-diffusion', 'why': 'The thermal diffusion factor for water in air is about 0.1 or less. Over pore-scale temperature differences of about 1 K at 200 K, it contributes a flux of order 1e-3 of the Fickian flux.'}, {'effect': 'Knudsen (rarefied-gas) transport', 'why': 'This is valid only while the mean free path is smaller than the pore size. Under lunar near-vacuum it is not, and then $D_v$ must be replaced by a Knudsen diffusivity. The form of the equation stays the same.'}, {'effect': 'Variation of ice density with temperature, and ice motion', 'why': 'The thermal expansion of ice is about 5e-5 per K, so changes in $\\rho_i$ are about 1e-4 of $\\rho_i$. Grains are assumed fixed.'}, {'effect': 'Liquid water and adsorbed water on regolith grains', 'why': 'Near 200 K there is no liquid phase. Adsorption on the grain surfaces is outside the ice/gas model.'}]}}

`chain: 600de2230d928b27601077c584933da05d10f1e22119ae9b38bb61ec0b7a8967`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S5: possibly undefined in the equation: R, x_{k}
- format: symbol Omega: symbol names in prose outside $...$: Omega

`chain: 7a3907c1c31802e6bb71e7c463c7e9ffc3a3be7a32a8b196d685b3f542da35ca`

## 7. verdict

- **PASS**

`chain: 0fa7b9b3c90e745bc2d29c59ba07702e853323376f11df29a3b9c2ff6b1304a6`
