# Work order: M5:Test lock (beth-verifier)

Run folder: {run}
Feature worktree: {worktree} (no new code exists yet; do not read its src/ or include/ after this stage)
Read `{run}/derivation/spec.md`, `code_map.json`, `mms/`, and `{run}/decisions.json`.
{feedback}
Write the tests in `{tests_dir}`: `tests.json` ({{"tests": [...]}}, the format in
`{root}/scripts/model_development/run_test.py`) and one check script per test. Each test also has
"tolerance" and "tolerance_reason" (for the user to approve), "layer" (from the spec's `## Layers`) and
"baseline": "main" when it must already pass on main.
{test_hint}
Run each baseline test on the worktree, which is still identical to main at this stage, with
`{root}/scripts/model_development/run_test.sh --tests {tests_dir}/tests.json --id <id> --project {main_project} --out {run}/tests --draft`
(one at a time; each obeys the run caps or run_test refuses it). Save the reference outputs a regression test
compares against (numbers, not whole run folders) in `{tests_dir}/reference/`, with the command that made them.

Done criterion: every test in the list has a check script; every baseline test passed on main; every test
fits the caps.
