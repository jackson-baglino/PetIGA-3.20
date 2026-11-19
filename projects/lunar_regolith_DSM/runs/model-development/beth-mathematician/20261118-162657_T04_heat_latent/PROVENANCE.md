# PROVENANCE: beth-mathematician on T04_heat_latent (run 20261118-162657_T04_heat_latent)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T16:26:57-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T04_heat_latent (task.json sha256 24b13f19bf2a, grader sha256 67bdb32e5d38)
- git commit: 2219d818d6; uncommitted definition files: none
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-162657_T04_heat_latent/workspace; settings: runs/model-development/beth-mathematician/20261118-162657_T04_heat_latent/settings.json
- file hashes: see manifest.json

`chain: f50953d510bc59895787349a731ee50a7356ebed8d83d975377526f5e38037ce`

## 2. agent call

- started 2026-11-18T16:26:57-0800, 78.5 s, exit 0, turns 5
- tokens: input 8, output 8443
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 1}
- transcript: transcript.jsonl

`chain: eeae4361ef20ed923e72312079089db38b28c526c0199f63eade4677a0082be4`

## 3. barrier events

- none

`chain: 5d6a323da8f1f92730eebc4c34d2dfff34f804adb775abec968dd8baea9f12a8`

## 4. answer

- status: done
- result_expr: C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), 
- steps: 6; questions: []
- answer.json sha256 b2b10ec254b4

`chain: 9b8fbba0a4db7ec30b30c38b7630087a146d3175d16673e51bbe660e9f261047`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units W / m ** 3'}, {'check': 'equilibrium', 'pass': True, 'details': 'R = 0 for uniform T and a static phi (3 draws, 2 points)'}, {'check': 'latent sign', 'pass': True, 'details': 'deposition warms and, by linearity, sublimation cools (3 draws)'}, {'check': 'conduction', 'pass': True, 'details': 'a hot spot cools at its peak (3 draws)'}, {'check': 'conservation', 'pass': True, 'details': 'conduction integrates to zero over a period (3 draws): it is the divergence of a flux'}], 'comparison (for the user)': {"agent's terms by role": {'storage': ['C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y)', '-Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z)'], 'phase change': ['-rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)'], 'other': []}, 'current model (src/assembly.c:183-185)': {'storage': ['rho(phi(x, y, z, t))*cp(phi(x, y, z, t))*Derivative(T(x, y, z, t), t)'], 'conduction': ['-(Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), x), x) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), y), y) + Derivative(k(phi(x, y, z, t))*Derivative(T(x, y, z, t), z), z))'], 'phase change': ['-rho_i*L*Derivative(phi(x, y, z, t), t)'], 'other': []}, "agent's symbols": [{'name': 'rho_i', 'latex': '\\rho_i', 'meaning': 'density of ice, taken constant (A3)', 'units': 'kg m^-3', 'typical': 920}, {'name': 'h_i(T)', 'latex': 'h_i(T)', 'meaning': 'specific enthalpy of ice', 'units': 'J kg^-1', 'typical': 300000}, {'name': 'h_v(T)', 'latex': 'h_v(T)', 'meaning': 'specific enthalpy of water vapor', 'units': 'J kg^-1', 'typical': 3100000}, {'name': 'h_a(T)', 'latex': 'h_a(T)', 'meaning': 'specific enthalpy of the non-condensable gas (air)', 'units': 'J kg^-1', 'typical': 200000}, {'name': 'rho_v', 'latex': '\\rho_v', 'meaning': 'mass of water vapor per unit gas volume', 'units': 'kg m^-3', 'typical': 1e-07}, {'name': 'm_a', 'latex': 'm_a', 'meaning': 'mass of air per unit total volume, $(1-\\phi)\\rho_a$', 'units': 'kg m^-3', 'typical': 0.01}, {'name': 'j_v', 'latex': '\\mathbf{j}_v', 'meaning': 'diffusive mass flux of water vapor', 'units': 'kg m^-2 s^-1', 'typical': 1e-09}, {'name': 'E', 'latex': 'E', 'meaning': 'enthalpy per unit total volume', 'units': 'J m^-3', 'typical': 300000000}, {'name': 'q', 'latex': '\\mathbf{q}', 'meaning': 'total energy flux (conduction plus enthalpy carried by vapor)', 'units': 'W m^-2', 'typical': 1000}, {'name': 'k(phi, T)', 'latex': 'k(\\phi,T)', 'meaning': 'effective thermal conductivity of the ice/gas mixture', 'units': 'W m^-1 K^-1', 'typical': 2}, {'name': 'c_i(T)', 'latex': 'c_i(T)', 'meaning': 'specific heat of ice, $dh_i/dT$', 'units': 'J kg^-1 K^-1', 'typical': 1500}, {'name': 'c_v(T)', 'latex': 'c_v(T)', 'meaning': 'specific heat (constant pressure) of water vapor, $dh_v/dT$', 'units': 'J kg^-1 K^-1', 'typical': 1850}, {'name': 'c_a(T)', 'latex': 'c_a(T)', 'meaning': 'specific heat of air, $dh_a/dT$', 'units': 'J kg^-1 K^-1', 'typical': 1000}, {'name': 'C(phi, T)', 'latex': 'C(\\phi,T)', 'meaning': 'effective volumetric heat capacity, $\\rho_i\\phi c_i + m_a c_a + (1-\\phi)\\rho_v c_v$', 'units': 'J m^-3 K^-1', 'typical': 1400000}, {'name': 'L(T)', 'latex': 'L(T)', 'meaning': 'specific latent heat of sublimation, $h_v(T)-h_i(T)$', 'units': 'J kg^-1', 'typical': 2830000}, {'name': 'R', 'latex': 'R', 'meaning': 'residual of the energy equation', 'units': 'W m^-3', 'typical': 0}], 'agent left out': [{'effect': 'Enthalpy carried by vapor diffusion, $c_v\\mathbf{j}_v\\cdot\\nabla T$', 'why': 'Its ratio to conduction is about $10^{-11}$ at pore scale for $\\rho_v\\sim10^{-7}$ kg m$^{-3}$ near 200 K (check_size).'}, {'effect': 'Pressure work and kinetic energy of the gas', 'why': 'The gas is at rest and at very low pressure. $Dp/Dt$ work is of order $\\rho_v$ times the vapor heat content, which is far below $\\rho_i c_i \\partial_t T$.'}, {'effect': 'Advection by bulk gas flow, and air mass change as pores open (if $m_a$ is not constant)', 'why': '$m_a\\ll\\rho_i$ (lunar pore gas is near vacuum), so any $h_a\\partial_t m_a$ term is smaller than the latent term by $m_a h_a/(\\rho_i L)\\lesssim10^{-6}$.'}, {'effect': 'Thermal expansion of ice (a variable $\\rho_i$)', 'why': 'The relative density change is about $10^{-4}$ per K times $\\Delta T$, which is small over the temperature range.'}, {'effect': 'Interface (surface) energy and curvature (Gibbs-Thomson) contribution to the energy', 'why': 'The surface energy per unit volume, $\\gamma/\\ell\\sim0.1/10^{-4}=10^3$ J m$^{-3}$, is far below the latent heat $\\rho_i L\\sim2.6\\times10^9$ J m$^{-3}$.'}, {'effect': 'Radiative heat transfer across pores', 'why': 'At about 200 K, radiative conductance $4\\epsilon\\sigma T^3\\ell\\sim2\\times10^{-4}$ W m$^{-1}$ K$^{-1}$ for $\\ell=10^{-4}$ m. This is far below $k$ of ice ($\\sim2$). It can be included in $k$ if needed.'}]}}

`chain: a631ed080ae5a9858649d3c1f02ac325b68538687fc9fc002c08fad180ee7a00`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: no problems

`chain: bd117c321ebee014eb24fdc04a5593bcae5be3b12ba285fdb623873a7ff67522`

## 7. verdict

- **PASS**

`chain: f874bf9ba393d67b4486ba9872c4d0398815f90d0e4ccbf3f55beea5f348bf60`
