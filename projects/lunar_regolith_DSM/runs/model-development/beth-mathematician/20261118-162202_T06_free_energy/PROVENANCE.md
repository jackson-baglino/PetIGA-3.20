# PROVENANCE: beth-mathematician on T06_free_energy (run 20261118-162202_T06_free_energy)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T16:22:02-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T06_free_energy (task.json sha256 47e48537f869, grader sha256 ca28d2ef724f)
- git commit: 4386cd5971; uncommitted definition files: none
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-162202_T06_free_energy/workspace; settings: runs/model-development/beth-mathematician/20261118-162202_T06_free_energy/settings.json
- file hashes: see manifest.json

`chain: 54b056a467d25fd8cf11702f7d6b3e6422e910edc9cb9c36d53e336bf119d46b`

## 2. agent call

- started 2026-11-18T16:22:02-0800, 70.4 s, exit 0, turns 5
- tokens: input 6, output 7788
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 1}
- transcript: transcript.jsonl

`chain: c79efdca0097110a35d825bf4413487c44cadf26f154ce5ce2e09370507bfb8d`

## 3. barrier events

- none

`chain: d7eae36063bedd3dd16e72a72d47d40d9a3802d1ad87e58871272c7f18945e95`

## 4. answer

- status: done
- result_expr: 6*gamma/eps*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2 + 3*gamma*eps/2*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2)
- steps: 7; questions: []
- answer.json sha256 3a1f251274e5

`chain: e3351aed62e6ddb33fb97d90ca58d6ac8435773c8dd78e24ee2a79dcdbccde02`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'f has units J / m ** 3 (needs an energy per volume)'}, {'check': 'wells', 'pass': True, 'details': 'minima at 0 and 1 with equal depths and a barrier between (3 draws)'}, {'check': 'gradient', 'pass': True, 'details': 'isotropic, quadratic and penalizing: kappa(p) |grad phi|^2 / 2 with kappa > 0 (3 draws)'}, {'check': 'interface', 'pass': True, 'details': 'the flat interface has energy gamma (3 draws); its 10-90 % width is 2.197 eps'}], 'comparison (for the user)': {"agent's f": '6*gamma/eps*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2 + 3*gamma*eps/2*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2)', 'implied by the current model (src/assembly.c:177-181)': '3*gamma*eps*(Derivative(phi(x, y, z, t), x)**2 + Derivative(phi(x, y, z, t), y)**2 + Derivative(phi(x, y, z, t), z)**2) + 3*gamma*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2/eps', "agent's symbols": [{'name': 'W', 'latex': 'W', 'meaning': 'height coefficient of the double-well (bulk) term; fixed to 6 gamma/eps in S5', 'units': 'J m^-3', 'typical': 600000}, {'name': 'kappa', 'latex': '\\kappa', 'meaning': 'gradient-energy coefficient; fixed to 3 gamma eps in S5', 'units': 'J m^-1', 'typical': 1.5e-09}, {'name': 'g', 'latex': 'g(\\phi)', 'meaning': 'dimensionless double-well function phi^2 (1-phi)^2', 'units': '1', 'typical': 0.0625}, {'name': 'X', 'latex': 'X', 'meaning': 'coordinate normal to a flat interface, with X=0 at phi=1/2', 'units': 'm', 'typical': 1e-08}, {'name': 'ell', 'latex': '\\ell', 'meaning': 'intrinsic decay length of the equilibrium profile, sqrt(kappa/(2W))', 'units': 'm', 'typical': 5e-09}, {'name': 'phi_0', 'latex': '\\phi_0(X)', 'meaning': 'equilibrium profile across a flat interface', 'units': '1', 'typical': 0.5}, {'name': 'sigma', 'latex': '\\sigma', 'meaning': 'excess free energy per unit area of a flat equilibrium interface', 'units': 'J m^-2', 'typical': 0.1}, {'name': 'delta', 'latex': '\\delta', 'meaning': 'interface width, defined as 1 / max |d phi_0/dX| (the distance over which the tangent at the midpoint rises from 0 to 1)', 'units': 'm', 'typical': 2e-08}], 'agent left out': [{'effect': 'Driving force for sublimation and deposition (tilting term, for example proportional to $(T-T_m)$ or vapor supersaturation, times $h(\\phi)$)', 'why': 'The work order is at coexistence, so this term is zero by assumption A1.'}, {'effect': 'Anisotropy of the ice surface energy (hexagonal ice facets)', 'why': 'Taken as isotropic (A2). The anisotropy of ice-vapor surface energy is a few percent, which is small compared with the geometric effects of pore confinement.'}, {'effect': 'Regolith grain surfaces: wetting and the ice-grain and gas-grain energies', 'why': 'The work order treats the pore space as ice plus gas only. Contact angles would enter through boundary terms, not through $f$.'}, {'effect': 'Coupling to temperature and vapor density fields (latent heat, vapor diffusion)', 'why': 'At constant uniform temperature and coexistence they do not change the equilibrium interface. They belong to the dynamics, not to this functional.'}, {'effect': 'Elastic strain energy in ice and density change on phase change', 'why': 'Ice in open pores against low-pressure gas is nearly stress-free. Elastic energy in the bulk is much smaller than $\\gamma/\\varepsilon$ for nanometer-to-micron $\\varepsilon$.'}, {'effect': 'Curvature (Gibbs-Thomson) effects in the calibration', 'why': 'Calibration uses a flat interface. Curved interfaces get corrections of relative order $\\varepsilon/R$, which are small when $\\varepsilon$ is much smaller than the pore radius $R$.'}, {'effect': 'Higher-order gradient terms (for example $(\\nabla^2\\phi)^2$)', 'why': 'They are excluded by the deliverable. They are of relative order $(\\varepsilon/L)^2$ on scales $L\\gg\\varepsilon$.'}]}}

`chain: d5bcdfdccfc49806b20b11181c6d213d01a4300740e0e5304144bd598d8db02a`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S1: possibly undefined in the equation: f
- format: S2: possibly undefined in the equation: g
- format: S3: possibly undefined in the equation: e, g
- format: S7: possibly undefined in the equation: f
- format: symbol W: symbol names in prose outside $...$: eps, gamma
- format: symbol kappa: symbol names in prose outside $...$: eps, gamma
- format: symbol g: symbol names in prose outside $...$: phi
- format: symbol X: symbol names in prose outside $...$: phi
- format: symbol ell: symbol names in prose outside $...$: kappa

`chain: 00f3c28db1d8981b44d333c716a940567826d92c29c01b97bff70247ff6f23e8`

## 7. verdict

- **PASS**

`chain: 3894d96265b1981eafc598efc807cad627521c7e86323ed98d59f5569c1fb1e4`
