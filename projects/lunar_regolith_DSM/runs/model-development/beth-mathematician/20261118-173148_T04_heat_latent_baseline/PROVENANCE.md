# PROVENANCE: beth-mathematician on T04_heat_latent (run 20261118-173148_T04_heat_latent_baseline)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T17:31:48-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T04_heat_latent (task.json sha256 24b13f19bf2a, grader sha256 67bdb32e5d38)
- git commit: affd1b2f2b; uncommitted definition files: ['?? projects/lunar_regolith_DSM/workflows/model-development/agents/beth-mathematician/tutorial/README.md']
- Claude Code: 2.1.327 (Claude Code)
- MODE: single-prompt baseline (--no-tools): no tools, work order in the prompt; tools: none; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-173148_T04_heat_latent_baseline/workspace; settings: runs/model-development/beth-mathematician/20261118-173148_T04_heat_latent_baseline/settings.json
- file hashes: see manifest.json

`chain: ef32d5a377ccd93ca9684a42a5091608e5f95321c2adde7cc44d3200b5157849`

## 2. agent call

- started 2026-11-18T17:31:48-0800, 58.6 s, exit 0, turns 2
- tokens: input 2, output 6768
- tool calls: none
- transcript: transcript.jsonl

`chain: 4a44adc905af1ee0ad67f7e68efb06169636e56f620425909c853aee2be5604f`

## 3. barrier events

- none

`chain: 0773c9ca63ebb7559f8f4876b38a0b312701c4fdb53c0e3399e1b3eef74fd4eb`

## 4. answer

- status: done
- result_expr: C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), 
- steps: 8; questions: []
- answer.json sha256 1467296a08d3

`chain: 889397a19d7beb353fcbbbce3826af105af33d990bd1a77bcb6ab63b51aecb29`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units W / m ** 3'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform T and a static phi (3 draws, 2 points)'}, {'check': 'latent sign', 'pass': True, 'details': 'deposition warms and, by linearity, sublimation cools (3 draws)'}, {'check': 'conduction', 'pass': True, 'details': 'a hot spot cools at its peak (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'conduction integrates to zero over a period (3 draws): it is the divergence of a flux'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z)'], 'phase change': ['-rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, 'current model (src/assembly.c:183-185)': {'storage': ['rho(phi(x, y, z, t))*cp(phi(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-(Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), x), x) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), y), y) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), z), z))'], 'phase change': ['-rho_i*L*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i', 'latex': '\\rho_i', 'meaning': 'density of ice (taken constant, see A4)', 'units': 'kg m^-3', 'typical': 920}, {'name': 'rho_v', 'latex': '\\rho_v', 'meaning': 'mass concentration of water vapor in the gas', 'units': 'kg m^-3', 'typical': 2e-06}, {'name': 'h_i(T)', 'latex': 'h_i(T)', 'meaning': 'specific enthalpy of ice', 'units': 'J kg^-1', 'typical': -300000}, {'name': 'h_v(T)', 'latex': 'h_v(T)', 'meaning': 'specific enthalpy of water vapor', 'units': 'J kg^-1', 'typical': 2500000}, {'name': 'c_i(T)', 'latex': 'c_i(T)', 'meaning': 'specific heat of ice, dh_i/dT', 'units': 'J kg^-1 K^-1', 'typical': 1400}, {'name': 'c_v(T)', 'latex': 'c_v(T)', 'meaning': 'specific heat of water vapor at constant pressure, dh_v/dT', 'units': 'J kg^-1 K^-1', 'typical': 1850}, {'name': 'H', 'latex': 'H', 'meaning': 'enthalpy per unit volume of the ice-gas mixture', 'units': 'J m^-3', 'typical': 1000000}, {'name': 'q', 'latex': '\\mathbf q', 'meaning': 'conductive heat flux', 'units': 'W m^-2', 'typical': 10}, {'name': 'J_v', 'latex': '\\mathbf J_v', 'meaning': 'diffusive mass flux of water vapor', 'units': 'kg m^-2 s^-1', 'typical': 1e-08}, {'name': 'm_dot', 'latex': '\\dot m', 'meaning': 'mass of vapor deposited as ice per unit volume and time (negative for sublimation)', 'units': 'kg m^-3 s^-1', 'typical': 1e-06}, {'name': 'k(phi, T)', 'latex': 'k(\\phi,T)', 'meaning': 'effective thermal conductivity of the ice-gas mixture', 'units': 'W m^-1 K^-1', 'typical': 2.5}, {'name': 'C(phi, T)', 'latex': 'C(\\phi,T)', 'meaning': 'volumetric heat capacity of the mixture, phi rho_i c_i + (1-phi) rho_v c_v (+ air if present)', 'units': 'J m^-3 K^-1', 'typical': 1300000}, {'name': 'L(T)', 'latex': 'L(T)', 'meaning': 'specific latent heat of sublimation, h_v - h_i', 'units': 'J kg^-1', 'typical': 2840000}, {'name': 'R', 'latex': 'R', 'meaning': 'residual of the energy equation', 'units': 'W m^-3', 'typical': 1000}], 'agent left out': [{'effect': 'Enthalpy carried by vapor diffusion, $c_v\\mathbf J_v\\cdot\\nabla T$', 'why': '$\\mathbf J_v\\sim D\\nabla\\rho_v$ with $\\rho_v\\sim10^{-6}$ kg m^-3 near 200 K. The ratio to conduction is about $c_v D\\rho_v/(k\\ell)\\cdot\\ell \\sim 1850\\cdot10^{-4}\\cdot10^{-6}/2.5\\sim10^{-7}$, and its gradient factor makes it smaller still.'}, {'effect': 'Gas convection (bulk flow of the gas in the pores)', 'why': 'The pores are small (tens of micrometers) and the pressures are low, so the Peclet number is far below 1.'}, {'effect': 'Pressure work and viscous dissipation', 'why': 'Pressure changes are tiny compared with $\\rho c T$, and the flow is negligible.'}, {'effect': 'Radiative heat transfer across pores', 'why': 'With pore size $d\\sim10^{-4}$ m, $4\\sigma T^3 d\\sim2\\times10^{-4}$ W m^-1 K^-1, much less than $k$ of ice ($\\sim2.5$).'}, {'effect': 'Interfacial (surface) energy and Gibbs-Thomson change of $L$', 'why': '$\\gamma/(\\rho_i L r)\\sim0.1/(920\\cdot2.8\\times10^6\\cdot10^{-6})\\sim4\\times10^{-5}$ relative to the latent heat.'}, {'effect': 'Thermal expansion of ice (variable $\\rho_i$)', 'why': 'The relative change is about $10^{-4}$ per 1 K, a negligible correction to $\\rho_i\\partial_t\\phi$.'}, {'effect': 'Non-equilibrium between ice and gas temperatures', 'why': 'The thermal diffusion time across a pore, about $d^2/\\alpha\\sim10^{-8}/10^{-6}=10^{-2}$ s, is far shorter than sublimation time scales.'}]}}

`chain: e57dd08055b1b4e17cdca8c0e56ba99e97d46010f17072ad0ec4090ceac41af4`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: symbol C(phi, T): symbol names in prose outside $...$: phi, rho_i, rho_v

`chain: 34dc776f53ba5501fd1cbdd067b4edf7468be189d49cfb48664182db7bad725a`

## 7. verdict

- **PASS**

`chain: e33283e370b4c91d3a39306761c109db7c4e4c574669baecfb7614c4dc11bb4c`
