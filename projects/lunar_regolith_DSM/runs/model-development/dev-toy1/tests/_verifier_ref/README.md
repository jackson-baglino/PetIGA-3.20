# Verifier's calibration projects (M5, DEV TOY DRY RUN)

`ok/toy/heat1d.py` is the verifier's own implementation of spec.md, written only to check that the M5 tests
pass on a correct implementation (no test is vacuously failing) and to reproduce the spec's scratch numbers
independently. Each `mut_*` project loads `ok` and plants one error, to check that the tests can fail.
It is kept outside the feature worktree so that the M6 implementation stays independent. Do not copy it.
