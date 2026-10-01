# DoF/core scaling test (2026-09-29)

**Setup.** 30 mid-run steps, restarted from the 3a φ 0.325 −20 °C run at
~10 d, with k_eff every 3 steps and `-log_view`. 6 targets × 2 repeats.
Data: `~/SimulationResults/HPC_results/enceladus_DSM/batch_2026-09-29__scaling/`.
The four batch folders there holding only `inputs/` are empty earlier
submissions: no jobs ran.

![scaling](scaling.png)

| ranks | DoF/core | nodes | phase-field s/step | k_eff s/sample | $/run −20 °C | $/run −5 °C |
|---|---|---|---|---|---|---|
| 121 | 198k | 4 | 23.4 | 5.2 | 3.85 | 12.1 |
| 161 | 149k | 6 | 19.6 | 4.4 | 4.28 | 13.5 |
| 241 | 100k | 8 | 19.8 | 75.3 | 16.7 | 38.0 |
| 401 | 60k | 13 | 15.8 | 94.5 | 30.4 | 64.5 |
| 601 | 40k | 19 | 18.4 | 102.9 | 50.4 | 108 |
| 1201 | 20k | 38 | 16.4 | 148.6 | 131 | 263 |

Costs per run use the SSA-trigger cadence: 179 / 308 samples at −20 / −5 °C.

**Findings.**

1. **The phase-field step does not strong-scale** here. From 121 to 1201
   ranks it stays at 16–23 s, so its cost in core-seconds grows almost
   linearly with rank count. The fewest ranks is cheapest.
2. **The k_eff solve falls off a cliff between 6 and 8 nodes**: 4–5 s per
   sample below, 75–149 s above, at the same iteration count. In
   `-log_view` the extra time is all communication:
   - CG global reductions (`VecTDot`): 12 s → 174 s, for the same number of
     calls;
   - neighbour exchange (`VecScatterEnd`): 94 s → 632 s;
   - multigrid restriction (`MatMultTranspose`): 7 s → 357 s.
3. **The likely cause is network placement.** `run_enceladus.sh` requests
   `--constraint='icelake|skylake|cascadelake'`, which lets SLURM mix node
   types in one job. The bracketed form, `[icelake|skylake|cascadelake]`,
   requires a single type.
   - This is a hypothesis; the test below settles it.
   - It would also explain batch 2's 10× swings within a run, and rev64
     being fast at 20 nodes.

**Decision.** The cheapest measured point is 198k DoF/core (121 ranks, 4
nodes): 3.4× cheaper than 60k at −20 °C and 5× cheaper at −5 °C. The −5 °C
run takes ~8.3 h.

**Follow-up result (2026-10-01): the cost valley is flat.**

| ranks | DoF/core | s/step | k_eff s/sample | $/run −20 °C | $/run −5 °C |
|---|---|---|---|---|---|
| 61 | 394k | 42.1 | 8.7 | 3.46 | 10.95 |
| 81 | 296k | 38.3 (1 repeat) | 6.8 | 4.13 | 13.13 |
| 121 | 198k | 23.4 | 5.2 | 3.85 | 12.12 |
| 161 | 149k | 19.6 | 4.4 | 4.28 | 13.48 |

- The 61/81-rank jobs needed 4 GB/core; at 1G the 61-rank ones were
  OOM-killed (see `scripts/lib/alloc.sh` `mem_per_cpu`).
- Refit over all eight points: t ≈ 12.7 s + 1765/P per step. The parallel
  share is larger than the first two-point fit (15.6 s + 940/P) suggested.
- Cost is flat to ±10% from 61 to 161 ranks. 61 ranks saves ~10% but makes
  every run ~1.8× longer (−5 °C: ~15 h vs ~8.3 h).
- **Kept 200k DoF/core (121 ranks, 4 nodes)**: bottom of the valley, with
  about twice the throughput. 400k (2 nodes) is the fallback if queueing
  becomes the bottleneck.
- The bracketed-constraint test was cancelled after a day pending on
  Priority. It is moot at ≤ 6 nodes; revisit only if a run must exceed 6.

**Follow-ups as submitted** (~$5 total):
- **Node-type constraint:**
  `./scripts/HPC/submit_scaling_test.sh --run <run> --targets 100000 --tag-suffix bracket --sbatch "--constraint=[icelake|skylake|cascadelake]"`
  If k_eff drops to ~5 s at 241 ranks, adopt the bracketed constraint in
  `run_enceladus.sh`; larger jobs are then safe too.
- **Fewer ranks still:**
  `./scripts/HPC/submit_scaling_test.sh --run <run> --targets "300000 400000"`
  81 and 61 ranks. Memory was 21.6% of the reservation at 121 ranks.
