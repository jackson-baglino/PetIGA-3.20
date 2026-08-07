# Determining the REV size for k_eff

## Why

`k_eff` is computed by homogenization, which assumes the domain is a
Representative Elementary Volume: large enough that the effective property has
stopped depending on domain size or on which random realization you drew. Run
below the REV and you get a number that looks fine, varies wildly seed to seed,
and is not a material property.

The current working size, 2 mm at `R_ave = 100 um` (`L/R_ave = 20`), was
inherited from Calonne et al. (2011) rather than measured on *these* packings,
which are gravity-deposited rather than jammed and have different throat
statistics. This determines it directly.

**REV size grows with porosity** — more void means larger pores and more
room for the microstructure to vary between realizations. So the most porous
case bounds the others, and that is where the effort goes.

## Scope

Measure on **phi = 0.400** (the most porous case in the study). If it converges
at some `L/R_ave`, every lower porosity is covered. Confirm with a single
spot-check at phi = 0.250 rather than repeating the sweep.

Work in **`L/R_ave`, never in metres.** The whole problem is scale-invariant:
`L = 50 um` at `eps = 45 nm` needs exactly the same 1572 elements as `L = 2 mm`
at `eps = 1.8 um`. A result in `L/R_ave` transfers to any grain size; a result
in millimetres does not.

---

## Stage 1 — Geometric convergence (no solver, ~minutes)

Porosity variance versus window size is the classical REV criterion and costs
nothing. It gives a lower bound on the REV and a sanity check before spending
HPC time.

For `L/R_ave` in {5, 10, 15, 20, 30, 40, 60, 80}, generate 12 seeds each at
phi = 0.400, `--periodic xy`:

```bash
python3 preprocess/generate_packing_gravity.py \
    --Lx <L> --Ly <L> --porosity 0.400 --mean-r 2.5e-6 --sigma-ln 0.5 \
    --periodic xy --seed <s> --no-preview \
    --out studies/rev_sizing/packings/LR<n>_seed<s>
```

Record from `metadata.json`, across seeds at each size: mean and standard
deviation of `porosity_achieved`, `local_density_cv`,
`max_void_radius_per_mean_r`, `coordination_number`, `half_domain_asymmetry`,
and the throat percentiles.

**Read:** the size at which the seed-to-seed standard deviation of porosity
falls below ~1% absolute, and `coordination_number` stops drifting. Expect this
to land *below* the true k_eff REV — porosity self-averages faster than a
transport property, which is why Stage 1 is a lower bound and not the answer.

**Turn off the acceptance gate for this stage** (`--max-void-ratio 99
--max-density-cv 99 --max-asymmetry 99 --no-percolation-gate`). The gate exists
to reject unrepresentative packings, so leaving it on censors exactly the
tail whose shrinkage with `L` is what is being measured.

## Stage 2 — Resolve the throat question first

Do not start Stage 3 before settling this, because it shifts every `k_eff`
number and would invalidate the sweep.

These packings have a much narrower throat tail than the jammed ones:

| | p25/R_ave | p50/R_ave |
|---|---|---|
| gravity, L = 2 mm, R = 100 um | 0.040 | 0.706 |
| gravity, L = 100 um, R = 2.5 um | 0.027 | 0.403 |
| jammed reference, L = 2 mm | 0.258 | 0.602 |

Drop-and-roll leaves many grains in near-tangency that the jamming relaxation
would have opened. The 1–99% diffuse band, `9.2*eps = 0.41 um`, is ~5.7x the
p25 throat, where the jammed reference had it at 0.64x. Where the band is wider
than a throat, phi never reaches 0 inside it, `ThermalCond` returns well above
`k_air` in a channel that ought to be air, and **k_eff reads high**.

Decide between:

1. **Accept it.** Argue the near-tangent throats are physically sintered necks,
   so a solid bridge there is right. Then quantify the bias: run `k_eff` at
   two `eps` values a factor of 2 apart on the *same* packing. If `k_eff` moves
   materially, it is resolution, not material, and this option is dead.
2. **Open the tail.** Add a `--min-throat` option that pushes apart pairs whose
   surface gap is below a threshold (a short `packing_lib.relax` with inflated
   radii). Set the threshold from the band, e.g. `gap >= 9.2*eps`, which is the
   condition for phi to actually reach 0 in the channel.

Option 2 is the recommendation: it makes the geometry honest at the eps in use
rather than asking the reader to trust an argument. It slightly perturbs
contact tangency, so report `n_contacts` and `coordination_number` before and
after.

## Stage 3 — k_eff convergence (HPC)

For each `L/R_ave` in {10, 15, 20, 30, 40, 60}, 8–12 seeds at phi = 0.400,
`--periodic xy` (homogenization wants a periodic cell), `T = -20 C`,
`eps` from `comp_eps.py` with `R_feat = R_ave/50` held **fixed in units of
R_ave** across all sizes — otherwise the sweep confounds domain size with
resolution.

Run through `-keff 1` at t = 0 only. No time integration is needed: this is a
property of the initial microstructure, and it makes each case cheap.

Sizing, from the scale-invariant mesh rule `Nx = ceil(sqrt(2)*L/eps)`:

| L/R_ave | Nx = Ny | elements | vs L/R = 20 |
|---|---|---|---|
| 10 | 786 | 0.6 M | 0.25x |
| 20 | 1572 | 2.5 M | 1x |
| 40 | 3143 | 9.9 M | 4x |
| 60 | 4715 | 22 M | 9x |

Cost grows as `(L/R_ave)^2`, so ~10 seeds at `L/R = 60` dominates the budget.
Start at the small end and stop as soon as the criterion below is met.

**REV criterion.** Report `k_eff` as the mean over seeds with its standard
deviation, both components of the diagonal (`k_xx`, `k_yy`) separately. The REV
is the smallest `L/R_ave` where **all three** hold:

- relative scatter `std(k_eff)/mean(k_eff)` across seeds `< 5%`;
- the mean has stopped drifting — it changes by less than its own scatter
  between consecutive sizes;
- `k_xx` and `k_yy` agree within their scatter. Gravity deposition is the one
  process here with a preferred direction, so a persistent gap is real
  anisotropy and means an isotropic `k_eff` is the wrong summary, not that the
  domain is too small.

The third is the reason to keep both components rather than averaging early.

## Stage 4 — Confirm and write down

- Spot-check phi = 0.250 at the chosen size and one size below. If the scatter
  criterion holds at both, the phi = 0.400 result is conservative and covers
  the study.
- Record the outcome as `L/R_ave`, with the metre value for the current grain
  size as a derived convenience.
- Commit `rev_sweep.csv` (size, seed, phi, k_xx, k_yy, and the geometric
  metrics), the plotting script, and a `README.md` with the convergence figure
  and the number, per the results rule in `CLAUDE.md`.

---

## Verification

- Every packing passes its own gate before entering Stage 3 (Stage 1 runs
  ungated on purpose; Stage 3 does not).
- `eps/R_ave` and `R_feat/R_ave` are identical across every size in the sweep —
  assert this from the `metadata.json` and the geometry `.opts` rather than
  trusting the driver.
- One case run at two `eps` a factor of 2 apart, to show the plateau is a
  property of the domain size and not of the resolution.
- `k_eff` at phi -> 1 (no grains) returns `k_air`, and a single-crystal domain
  returns `k_ice`. Cheap, and catches a mis-scaled corrector before the sweep.

## What is already done

The 12 sample packings at `L/R_ave = 40`, phi = 0.325, R_ave = 2.5 um, one per
periodicity mode x 3 seeds, are in `inputs/packings/periodic_*/`. They are the
2x-REV-guess build, and the `--periodic xy` three are directly reusable as the
`L/R_ave = 40` point of Stage 3 once Stage 2 is settled.

Runs go to HPC via `scripts/HPC/submit_batch.sh`, not local `mpiexec`.
