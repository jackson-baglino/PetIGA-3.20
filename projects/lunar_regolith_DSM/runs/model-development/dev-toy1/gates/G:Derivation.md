# G:Derivation, M3 try 3: toy `heat1d-fd-be` (DEV TOY DRY RUN, not citable)

## 1. What M3 produced
- `derivation/ledger.json`: 57 steps, made up of 8 assumptions (S1–S8), 1 choice (S19 = C1) and 48 derived steps. No step ID was added, removed or renumbered since try 2. The try-3 edits are:
  - S1 now states that T_ref is constant in x and t, and S13 cites S1 for this.
  - S12 now lists S5 as a premise.
  - S21 and S22 now list S6.
  - S28 now lists S21 and S22. Nobody asked for this change.
- `derivation/derivation.md`: the derivation in prose. §0 maps each try-3 change to the step it affects.
- `derivation/scratch/m3_series.py`, `m3_try3.py` and `m3_ledger_check.py` were run at try 3. The modeler says every number is quoted from their output.
- Main results are unchanged:
  - S15: ρc_p T_t = k T_xx + q.
  - S25: θ_τ = θ_ξξ + 1, with zero data (S26).
  - S44: T_s = T_0 + qx(L−x)/(2k).
  - S45: F_0 = F_L = qL/2.
  - Uniqueness (S36), θ ≥ 0 (S50) and the series solution (S54).
  - S32 makes claims only for 0 < τ ≤ 0.1 and τ ≥ 1.

## 2. Checks (`checks/results_M3.json`, 2026-11-16 18:36)
- `all_pass: true`. All 12 scripts exit 0 with empty stderr.
- The file has 178 claim records and none fails. These are 119 checks plus 59 planted errors, and the self-test caught all 59. The verifier's summary says "118 claims". Excluding the environment probe from the 119 checks would give 118, but this is not stated.
- All 48 derived steps are covered by at least one check.
- The two failures from try 2 now pass:
  - `required_premises`
  - `tref_constant_for_S13`
- Independent numbers:
  - min of g on (0, 0.1] is 0.4744875;
  - the τ = 1 bounds are 6.5856e-5 and 6.6726e-6;
  - the series and the finite-difference (FD) solution differ by at most 4.7e-8;
  - the FD steady state differs from S44 by 9.9e-14 (relative).
- **Failures: none.**
- Small oddity: inside the passing T_ref check, the detail field "S13 operation attributes T_ref constancy to S1" reads `false`. S13's operation text does say "(all three by S1)". This looks like a text-matching problem in the script, and it does not affect the pass.

## 3. Claims from papers to confirm
None. The ledger cites no paper, as the brief requires.

## 4. Decisions for you
- **D-e. Pass G:Derivation, or send M3 back?**
  - Options: (a) pass and go on to M4; (b) send back.
  - The modeler made no recommendation. No check fails, and no open item blocks the gate.
- **D-f. Keep the extra premises S21 and S22 on S28?** The note did not ask for them.
  - Options: (a) keep; (b) revert to the try-2 list.
  - The modeler made no recommendation. It added them so S28 matches how S34 cites S9 ("as in" references). The verifier rated this pattern informative only.
- **D-g. Should the line record the earlier decisions in `decisions.json`?** The ledger cites D2(a), D3(a), D4(a) and the note "C1, C2 and the notation stand" from try 2, but `decisions.json` has no entry for any of them.
  - Options: (a) add entries so those citations can be checked; (b) leave as is.
  - The modeler asks this question but makes no recommendation.
- **C1, C2, D-b(a) and D-c(a)** were decided by the main session "for plumbing only (not a user decision)". Do you ratify them for the toy?
- **D-d, a reminder for M4:** S57 shows that the run to τ_end = 1 can still differ from steady state by up to 6.7e-6 in θ (5.3e-5 relative to the peak value 1/8). Any τ_end-vs-S44 tolerance tighter than about 1e-5 in θ would fail for this reason alone. Please confirm that M4 tolerances will respect this limit.

## 5. Not verifiable from the files
- D-d is not recorded in `decisions.json`. The only G:Derivation entry there is the try-3 send-back.
- The last section of derivation.md still says the verifier "has not yet been rerun on this revision". `results_M3.json` is timestamped after the try-3 send-back and says it was run against revision 3, so that sentence is out of date.
- The `planted` array at the end of `results_M3.json` lists only 3 entries, while the self-test section lists 59.
- S54 relies on three standard theorems that are named but not proved here: the Weierstrass M-test, term-by-term differentiation, and uniform convergence of Fourier series.
- S56's sign pattern for g′ is proved analytically. The scripts only sample it at 4 points.
