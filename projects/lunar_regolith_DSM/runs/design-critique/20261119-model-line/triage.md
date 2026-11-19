# Triage of the beth-critic review of the model line (run 20261119-model-line)

Date: 2026-11-19. Critic: beth-critic, `claude-opus-5-5`, 272 s, one attempt. Input: the bundle
`runs/design-critique/bundles/model_line_20261119.md` (38 files; the critic read all 6170 lines) and the
focus questions `runs/design-critique/focus/model_line_20261119.md`. Output: `critiques/beth-critic.json`
(21 findings, 7 unclear).

The main session checked each finding against the files. Verdicts:
- **confirmed**: the check below reproduces the claim;
- **partly**: the core holds, but part of the claim is wrong;
- **inferred**: a judgment I agree with that no script can check.

The two URLs come from real web searches. The transcript has four WebSearch calls with results; the
manifest's `web_search_requests: 0` is wrong. Neither page was opened by me.

| # | Sev. | Verdict | My check | Who acts |
|---|---|---|---|---|
| F1 | crit | partly | Runs 124517, 162816 and 173013 of T07 all PASS with different couplings: $\alpha \propto \rho_iR_vT/\rho_{vs}$ (local $T$), constant $\alpha$, and $30\alpha$. The growth, bulk and relaxation checks hold $T$ uniform (`T07/grade.py`, `phi_rate` default `T_body`). Wrong part: D1's source run *is* recorded (`studies/derived_model/README.md:19`, run 124517). | grader fix; T09 |
| F2 | crit | confirmed | The work orders, graders, planted errors, D1 and doc 30 were all written by the main session (Opus 5.5), and each `key_check` requires the code's row to pass. Doc 29 §3 item 10 already names the shared blind spot; nothing removes it. | user (sign-off); script (systematic mutations) |
| F3 | crit | confirmed | `check.printed` is the agent's own text. Nothing compares it with what `check.py` printed (`render_derivation.py:101,191` only checks that the function exists). Lab Report 3's "20 of 25 steps carry a run check" rests on it. | harness fix; recheck LR3 claim from transcripts |
| F4 | crit | confirmed | Plan item 3 ("the same agent writes code") contradicts doc 27 line 55 ("The math agent cannot write C") and `WORKFLOW.md` §M5 ("written … by an agent that never sees the code"). An FD Jacobian and an MMS built from the same equations cannot catch a wrong equation. | user decision |
| F5 | crit | partly | Wrong part: the approved design *does* lock physics tests (`WORKFLOW.md` §M5: closed-box water, free-energy decay, 1D front speed, 90° contact, separate refinement). True part: plan item 3 drops them, and no list has Kelvin equilibrium ($\sigma = d_0/R$) or an analytic planar front speed. | user (test list) |
| F6 | crit | confirmed | `t08_sizes.py` and the G3 prototype are in the session scratchpad. The reference run is in `~/SimulationResults/.../scratch/`. Doc 31 tags their results [checked]. | me: commit them |
| F7 | major | confirmed | Doc 31 §6 is wrong. T08 *gives* the vapor equation, so it cannot catch an error the agent makes in T05. | me: correct doc 31; new T05 option |
| F8 | major | confirmed | `WORKFLOW.md` M1 and G:Scope ("air or lunar laws") were never held. The code's $\rho_{vs}$ is the air-at-1-atm mixing form (`material_properties.c:180`), 4 % below pure vapor at −20 °C (8.49e-4 against 8.84e-4 kg m$^{-3}$, Murphy & Koop 2005, script in this session). | user decision |
| F9 | major | inferred | The plan is serial; the DDM comes last. | user decision |
| F10 | major | confirmed | The definition says "a grader you never see checks `result_expr` against a known answer" (`beth-mathematician.md:23`). The T07 work order gives F, the driving variable and no nucleation. The no-tool baseline passes 4/4. | wording (LR4); DDM work orders |
| F11 | major | inferred | The map $M = M_{\text{code}}/(2\gamma)$ absorbs the thin-interface correction (on by default; `outp.txt` "tau_sub terms … -thin_iface_corr 1"). "Same coupling [checked]" holds only under the map at $T_0$. | me: downgrade the claim in doc 30 |
| F12 | major | confirmed | `blind.py`: `T_SCALE = 200`, `RHOV_SCALE = 1e-6`. | grader fix |
| F13 | major | confirmed | `stop_check.py:11` "Blocks at most once"; it requires only `scratch/check.py`. | harness fix |
| F14 | minor | confirmed | `launch.py:441-453` records uncommitted files and does not refuse to run. | harness fix |
| F15 | major | confirmed | Grader stderr is used only on failure (`launch.py:300-305`). No SymPy or pint version and no lock hash in the manifest. | harness fix |
| F16 | major | confirmed | The audit reads the model from the init event only (`launch.py:253,281`). `modelUsage` is not checked, and that is where the 2026-11-16 Fable fallback showed. | harness fix |
| F17 | major | confirmed | A flaw in my T08 proposal: a group can hide inside a "dimensionless function", and G5 is "at most", not "equal to", the Buckingham count. | me: fold into doc 31 |
| F18 | minor | confirmed | `docs/material_parameters.md` §2.2: $\rho_{vs}(200\,\mathrm K) \approx 4\times10^{-8}$. Murphy & Koop gives $p = 0.163$ Pa, so $\rho_{vs} = 1.76\times10^{-6}$ kg m$^{-3}$, 43× higher. The code's $p$ agrees with Murphy & Koop to 4 digits. | registry |
| F19 | minor | confirmed | The Claude Code version is in the barrier fingerprint (`launch.py:99`) and is not pinned. | harness |
| F20 | minor | confirmed | `launch.py` has no retry, no usage-limit exit and no resume. | harness |
| F21 | minor | confirmed | The T07 PROVENANCE lists "S2: possibly undefined … a, e, l, r". These look like lint false positives (letters of a LaTeX word). Not checked. | lint fix |

Unclear items:
1. **Which document is the approved design.** `WORKFLOW.md` is out of date: it names beth-modeler, beth-verifier and the foreman, while what was built is beth-mathematician under docs 27–29. Real; for the user.
2. **Doc 27 against plan item 3.** Same as F4.
3. **README label.** The README calls T07 "known answer" (`beth-mathematician/README.md:90`), while elsewhere it is key-free. The label comes from `agent_readme.py`. Real; minor.
4. **The 8-block cap on the Stop hook.** Not checked.
5. **Read isolation.** The critic says Claude Code's OS sandbox covers Bash and its children, not the built-in Read tool, so Read isolation rests on permission denies and `guard.py`. Plausible; not checked by me against the documentation.
6. **$\sigma_0$ in doc 31.** $d_0/R_g = 6.8\times10^{-5} > \sigma_0 = 10^{-5}$. $\sigma_0$ was inferred from experiment names and should be removed or defined. Agreed.
7. **Doc 26's status.** Minor.
