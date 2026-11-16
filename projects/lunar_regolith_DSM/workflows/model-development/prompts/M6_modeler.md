# Work order: M6:Implementation, layer "{layer}" (beth-modeler)

Run folder: {run}
Feature worktree (your only write area): {worktree}
Read `{run}/derivation/spec.md` and `code_map.json`; the tests of this layer are listed in
`{tests_dir}/tests.json` (read only; you may not edit tests).
{feedback}
Implement layer "{layer}" behind the runtime option the code map names, default off. Then run, one at a
time, the tests of this layer and of every earlier layer:
`{root}/scripts/model_development/run_test.sh --tests {tests_dir}/tests.json --id <id> --project {worktree_project} --out {run}/tests`
Tests for this layer: {layer_tests}. Earlier layers: {earlier_tests}.

Done criterion: the code builds and every listed test passes. If a test fails after your best fix, return
`blocked` with the evidence; if you believe a test is wrong, return `dispute` with the evidence.
