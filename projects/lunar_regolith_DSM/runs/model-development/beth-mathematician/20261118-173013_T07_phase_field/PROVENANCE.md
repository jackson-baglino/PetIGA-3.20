# PROVENANCE: beth-mathematician on T07_phase_field (run 20261118-173013_T07_phase_field)

Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json.

## 1. run started

- time: 2026-11-18T17:30:13-0800
- agent: beth-mathematician, model claude-opus-5-5, workflow model-development
- definition: workflows/model-development/agents/beth-mathematician.md (sha256 6e3a126ffb15)
- task: T07_phase_field (task.json sha256 c328b5b87ff4, grader sha256 a6333c7fee7e)
- git commit: affd1b2f2b; uncommitted definition files: ['?? projects/lunar_regolith_DSM/workflows/model-development/agents/beth-mathematician/tutorial/README.md']
- Claude Code: 2.1.327 (Claude Code)
- tools: Read, Write, Edit, Bash; disallowed: WebFetch, WebSearch, Agent, NotebookEdit
- working folder: runs/model-development/beth-mathematician/20261118-173013_T07_phase_field/workspace; settings: runs/model-development/beth-mathematician/20261118-173013_T07_phase_field/settings.json
- file hashes: see manifest.json

`chain: e2f69f5785db591b3e54c6aca75a23d4f4353411f3011c814b4f76edbb332132`

## 2. agent call

- started 2026-11-18T17:30:13-0800, 92.7 s, exit 0, turns 7
- tokens: input 10, output 9464
- tool calls: {'Read': 1, 'Write': 1, 'Bash': 2, 'Edit': 1}
- transcript: transcript.jsonl

`chain: 2c984f2a70e59e5118067bd5744073f0e40bca6ab163afff44a54b5f21edebaa`

## 3. barrier events

- none

`chain: 68a670ebd97c9331c109a2879794bbda33ceb2345ceeb6a79a1d0b6b6bc3e99f`

## 4. answer

- status: done
- result_expr: Derivative(phi(x, y, z, t), t) + M*(6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t)) - 6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2)))) - alpha*30*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v
- steps: 6; questions: ['The interface function $p(\\phi)=30\\phi^2(1-\\phi)^2$ is my choice. Any localized $p$ with $\\int_0^1p\\,d\\phi=1$ gives the same sharp-interface speed. If your model uses a different one, such as $\\phi^2(1-\\phi)^2$ without the 30, $\\alpha$ has to be rescaled to match.', 'Should $M$ be tied to $\\alpha$ through the Gibbs-Thomson relation, so that curvature-driven sublimation agrees with the Kelvin equation? Or should it stay a free parameter?']
- answer.json sha256 162081c6fbd7

`chain: eb759730ac8f95a01cce63df94b8102ac745a0c6e68640ca4adb6ead283f1b2c`

## 5. grade

- pass: True
- details: {'checks': [{'check': 'units', 'pass': True, 'details': 'every term has units 1 / s'}, {'check': 'bulk', 'pass': True, 'details': 'pure gas and pure ice stay put above and below saturation (3 draws)'}, {'check': 'growth', 'pass': True, 'details': 'supersaturation deposits ice and undersaturation sublimates it (3 draws)'}, {'check': 'relaxation', 'pass': True, 'details': 'at saturation phi_t = -M dF/dphi with M > 0, so F decreases (3 draws, phi = 0.3 and 0.7)'}], 'comparison (for the user)': {"agent's terms by role": {'rate': ['Derivative(phi(x, y, z, t), t)'], 'relaxation': ['M*(-6*eps*gamma*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2))) + 6*gamma*(1 - 2*phi(x, y, z, t))*(1 - phi(x, y, z, t))*phi(x, y, z, t)/eps)'], 'driving force': ['-30*alpha*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))*phi(x, y, z, t)**2'], 'other': []}, 'current model (src/assembly.c:177-181)': {'rate': ['Derivative(phi(x, y, z, t), t)'], 'relaxation': ['- 3*M*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2))) + 3*M/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t))'], 'driving force': ['-alpha/rho_i*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))'], 'other': []}, "agent's symbols": [{'name': 'M', 'latex': 'M', 'meaning': 'interface mobility (Allen-Cahn kinetic coefficient): how fast $\\phi$ relaxes down the gradient of $F$; with it a curved interface moves at $v_n=-6\\gamma\\varepsilon M\\kappa$', 'units': 'm^3 J^-1 s^-1', 'typical': 1e-06}, {'name': 'alpha', 'latex': '\\alpha', 'meaning': 'deposition/sublimation rate coefficient; a flat interface advances into the gas at $v_n = 6\\varepsilon\\alpha(\\rho_v-\\rho_{vs})$, so $6\\varepsilon\\alpha = k/\\rho_i$ with $k$ a kinetic mass-transfer velocity (Hertz-Knudsen: $k=a_s\\sqrt{k_BT/2\\pi m}$) and $\\rho_i$ the ice density', 'units': 'm^3 kg^-1 s^-1', 'typical': 1000}, {'name': 'p(phi)', 'latex': 'p(\\phi)', 'meaning': 'interface-localizing function $30\\phi^2(1-\\phi)^2$, zero in bulk ice and gas, $\\int_0^1 p\\,d\\phi = 1$', 'units': '1', 'typical': 1}, {'name': 'mu(phi)', 'latex': '\\mu', 'meaning': 'variational derivative $\\delta F/\\delta\\phi$', 'units': 'J m^-3', 'typical': 1000}, {'name': 'phi_0(n)', 'latex': '\\phi_0(n)', 'meaning': 'equilibrium flat-interface profile along the normal coordinate $n$', 'units': '1', 'typical': 0.5}, {'name': 'n', 'latex': 'n', 'meaning': 'coordinate normal to a flat interface, pointing into the ice', 'units': 'm', 'typical': 1e-06}, {'name': 'v_n', 'latex': 'v_n', 'meaning': 'normal velocity of the interface, positive when ice grows', 'units': 'm s^-1', 'typical': 1e-09}, {'name': 'kappa', 'latex': '\\kappa', 'meaning': 'sum of principal curvatures of the interface (positive for a convex ice grain)', 'units': 'm^-1', 'typical': 10000}, {'name': 'Delta', 'latex': '\\Delta', 'meaning': 'supersaturation $\\rho_v-\\rho_{vs}(T)$', 'units': 'kg m^-3', 'typical': 1e-08}], 'agent left out': [{'effect': 'Anisotropy of ice surface energy and kinetics (faceting)', 'why': 'Taken as isotropic: $\\gamma$, $M$ and $\\alpha$ are scalars. This is acceptable for rounded pore-scale grains; anisotropy would add orientation-dependent factors.'}, {'effect': 'Temperature dependence of $M$ and $\\alpha$', 'why': 'Kept constant for simplicity. Their Arrhenius or $\\sqrt{T}$ variation can be restored by writing $M(T)$ and $\\alpha(T)$ without changing the form of the equation.'}, {'effect': 'Nonlinear kinetics in large supersaturation', 'why': 'The pore vapor is close to saturation ($|\\rho_v-\\rho_{vs}|/\\rho_{vs}\\ll1$), so the linear term is the leading order and the dropped terms are $O((\\Delta/\\rho_{vs})^2)$.'}, {'effect': 'Curvature dependence of the saturation density (explicit Gibbs-Thomson term in $\\rho_{vs}$)', 'why': 'It is represented by the $M$-term ($v_n\\propto-\\kappa$). Matching $6\\gamma\\varepsilon M=6\\varepsilon\\alpha\\rho_{vs}\\gamma/(\\rho_i R_vT)\\cdot$(molar factors) would make it quantitative; that choice is left to the calibration of $M$. The size is $\\rho_{vs}\\gamma/(\\rho_iR_vT\\,r)\\sim10^{-4}\\rho_{vs}$ for $r\\sim1\\,\\mu$m.'}, {'effect': 'Mechanical stress, elasticity and grain motion of ice and regolith', 'why': 'The regolith grains and ice are treated as rigid; stresses are outside this model.'}, {'effect': 'Thermal noise / fluctuations', 'why': 'Deterministic model; noise is needed only for nucleation, which is excluded.'}, {'effect': 'Wetting of regolith grains by ice (contact angle at grain walls)', 'why': 'Boundary conditions at the grain surfaces are taken as natural (no-flux), as needed to drop the boundary term in S1; a wetting energy would add a surface term.'}]}}

`chain: e5eba4270bcc1d5a61b6097e599d88d62af27183e2caac5e46d10948647aa107`

## 6. derivation.md (written by render_derivation.py from answer.json; format problems do not change the verdict)

- format: S2: possibly undefined in the equation: a, e, l, r
- format: S3: possibly undefined in the equation: a, e, h, s
- format: S4: possibly undefined in the equation: R
- format: S5: possibly undefined in the equation: e

`chain: 0cdf78b328eab6d3dff3af3e96ae2648e81652a7800d27b7fdc2c4073e917e05`

## 7. verdict

- **PASS**

`chain: 4e5ffc177f695f4ed9ac3ec389b29c853a75dd0967b22cc77534960c5a62a772`
