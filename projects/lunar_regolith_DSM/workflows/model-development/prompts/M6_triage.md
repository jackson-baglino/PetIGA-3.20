# Work order: M6 triage (beth-model-foreman)

Run folder: {run}
Test {test_id} still fails after {rounds} fix rounds. Read: the test entry in `{tests_dir}/tests.json` and
its check script, the latest results in `{run}/tests/{test_id}/`, the spec `{run}/derivation/spec.md`, the
code diff `{run}/triage/diff_{n}.patch`, the modeler's last result `{run}/triage/modeler_{n}.json`, and
earlier triage records in `{run}/triage/`.
{dispute}
Write `{run}/triage/{n}.json` with the same content as your JSON answer (route, reason, evidence,
next_work_order, alternatives).
