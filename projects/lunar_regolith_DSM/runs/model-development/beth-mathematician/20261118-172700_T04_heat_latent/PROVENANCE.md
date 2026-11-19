# PROVENANCE: beth-mathematician on T04_heat_latent (run 20261118-172700_T04_heat_latent)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T17:27:00-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T04_heat_latent (task.json sha256 24b13f19bf2a, grader sha256 67bdb32e5d38)
- git commit: affd1b2f2b; uncommitted definition files: none
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-172700_T04_heat_latent/workspace; settings: runs/model-development/beth-mathematician/20261118-172700_T04_heat_latent/settings.json
- file hashes: see manifest.json

`chain: e38cea0f22fa18e9d59f5f9e1f787d21264d570bdbe183a94249dce9a672280c`

## 2. agent call

- started 2026-11-18T17:27:00-0800, 76.1 s, exit 0, turns 5
- tokens: input 8, output 8819
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 1}
- transcript: transcript.jsonl

`chain: 6cc4f6504c5f8a2dd21b84bad286afac13231c53a8c4bb804c9f31899dcbdce7`

## 3. barrier events

- none

`chain: 8d676a30cb069362f5a399a7a72ae96f56140c8d3aea23019bc900ef7715ac04`

## 4. answer

- status: done
- result_expr: C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), 
- steps: 6; questions: []
- answer.json sha256 0947efda403d

`chain: 26dca0d0c41baf902130de9390d6e741569126d10099ccb5bb2db00467b0b298`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units W / m ** 3'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform T and a static phi (3 draws, 2 points)'}, {'check': 'latent sign', 'pass': True, 'details': 'deposition warms and, by linearity, sublimation cools (3 draws)'}, {'check': 'conduction', 'pass': True, 'details': 'a hot spot cools at its peak (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'conduction integrates to zero over a period (3 draws): it is the divergence of a flux'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z)'], 'phase change': ['-rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, 'current model (src/assembly.c:183-185)': {'storage': ['rho(phi(x, y, z, t))*cp(phi(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-(Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), x), x) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), y), y) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), z), z))'], 'phase change': ['-rho_i*L*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i', 'latex': '\\rho_i', 'meaning': 'density of ice (constant, A3)', 'units': 'kg m^-3', 'typical': 920}, {'name': 'rho_a', 'latex': '\\rho_a', 'meaning': 'partial density of dry air in the gas', 'units': 'kg m^-3', 'typical': 0.001}, {'name': 'rho_v(x, y, z, t)', 'latex': '\\rho_v', 'meaning': 'partial density of water vapor in the gas', 'units': 'kg m^-3', 'typical': 1e-06}, {'name': 'h_i(T)', 'latex': 'h_i(T)', 'meaning': 'specific enthalpy of ice', 'units': 'J kg^-1', 'typical': -280000}, {'name': 'h_a(T)', 'latex': 'h_a(T)', 'meaning': 'specific enthalpy of air', 'units': 'J kg^-1', 'typical': 200000}, {'name': 'h_v(T)', 'latex': 'h_v(T)', 'meaning': 'specific enthalpy of water vapor', 'units': 'J kg^-1', 'typical': 2550000}, {'name': 'c_i(T)', 'latex': 'c_i(T)', 'meaning': 'specific heat of ice, dh_i/dT', 'units': 'J kg^-1 K^-1', 'typical': 1400}, {'name': 'c_a(T)', 'latex': 'c_a(T)', 'meaning': 'specific heat of air at constant pressure, dh_a/dT', 'units': 'J kg^-1 K^-1', 'typical': 1000}, {'name': 'c_v(T)', 'latex': 'c_v(T)', 'meaning': 'specific heat of water vapor at constant pressure, dh_v/dT', 'units': 'J kg^-1 K^-1', 'typical': 1850}, {'name': 'L(T)', 'latex': 'L(T)', 'meaning': 'latent heat of sublimation, h_v - h_i', 'units': 'J kg^-1', 'typical': 2830000}, {'name': 'C(phi, T)', 'latex': 'C(\\phi,T)', 'meaning': 'effective volumetric heat capacity, \\phi\\rho_i c_i + (1-\\phi)(\\rho_a c_a + \\rho_v c_v)', 'units': 'J m^-3 K^-1', 'typical': 1300000}, {'name': 'k(phi, T)', 'latex': 'k(\\phi,T)', 'meaning': 'effective thermal conductivity of the ice-gas mixture', 'units': 'W m^-1 K^-1', 'typical': 2}, {'name': 'H', 'latex': 'H', 'meaning': 'enthalpy per unit volume of the mixture', 'units': 'J m^-3', 'typical': -260000000}, {'name': 'q', 'latex': '\\mathbf{q}', 'meaning': 'conductive heat flux', 'units': 'W m^-2', 'typical': 1}, {'name': 'j_v', 'latex': '\\mathbf{j}_v', 'meaning': 'mass flux of water vapor', 'units': 'kg m^-2 s^-1', 'typical': 1e-07}, {'name': 'j_a', 'latex': '\\mathbf{j}_a', 'meaning': 'mass flux of air', 'units': 'kg m^-2 s^-1', 'typical': 1e-08}], 'agent left out': [{'effect': 'Advection of sensible heat by vapor and air fluxes', 'why': 'Péclet number $|\\mathbf{j}|c\\ell/k \\approx 10^{-7}\\cdot1850\\cdot10^{-4}/0.02 \\approx 10^{-6}$.'}, {'effect': 'Pressure work and changes in total pressure', 'why': 'Total pressure in the pore is nearly constant on the time scales of sublimation. $\\partial_t p$ is many orders of magnitude smaller than $C\\,\\partial_t T$.'}, {'effect': 'Kinetic energy and viscous dissipation of the gas', 'why': 'Gas speeds are tiny (mm/s or less), so $\\rho u^2/2$ is negligible next to the enthalpies (about $10^6$ J/kg).'}, {'effect': 'Thermal expansion of ice (temperature-dependent $\\rho_i$)', 'why': 'The expansion coefficient is about $10^{-4}$ K$^{-1}$, so the extra terms are of relative order $10^{-4}\\,\\Delta T$.'}, {'effect': 'Interfacial (surface) energy and curvature (Gibbs-Thomson) correction to the latent heat', 'why': '$\\gamma/(\\rho_i L\\,r) \\approx 0.1/(920\\cdot2.8\\times10^6\\cdot10^{-6}) \\approx 4\\times10^{-5}$ for curvature radius $r \\approx 1\\,\\mu$m. The phase-field equation for $\\phi$ is a later task.'}, {'effect': 'Thermal radiation across pores', 'why': 'Effective radiative conductivity $4\\sigma T^3\\ell \\approx 4\\cdot5.7\\times10^{-8}\\cdot8\\times10^6\\cdot10^{-4} \\approx 2\\times10^{-4}$ W/(m K), small next to the conductivity of ice (about 2-4 W/(m K)). It matters only in nearly empty, large pores.'}, {'effect': 'Dufour effect and Soret effect', 'why': 'These cross effects are second order in dilute vapor diffusion.'}, {'effect': 'Heat exchange with the regolith grains', 'why': 'The task restricts the domain to the pore space. The grains enter only through boundary conditions.'}]}}

`chain: 9529e98f34e6434c8cee36f689a5cee36119dfd3655758c85f475139da15aeb9`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S6: possibly undefined in the equation: R
- format: symbol C(phi, T): symbol names in prose outside $...$: phi, rho_a, rho_i, rho_v

`chain: 0006c4290915cf02836b9226efc395f47d0a1cd0df8bc38bc8e7e08ea093c571`

## 7. verdict

- **PASS**

`chain: e03d983f8c70f16b7fcceb483d36f91caf3c726ecba68d90d684f4451da6b304`
