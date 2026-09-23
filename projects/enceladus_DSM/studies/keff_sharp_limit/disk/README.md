# k_eff tensor-law verification: disk array

Tests `-keff_interp tensor` on a curved interface, against a sharp closed form.

A single ice disk (`R = 250 µm`) in a 1 mm periodic cell is a square array of
disks at area fraction `f = 0.19635`. Its sharp conductivity is Rayleigh's
(Perrins, McKenzie & McPhedran 1979): `k = 0.0295684 W/m/K`, and Maxwell-Garnett
alone is within 1.8e−4 of that, so the reference is not the limiting error.

## Why this case

- **It exercises the whole tensor.** The laminate's normal is exactly `ŷ`, so
  there the tensor law is diagonal and an error in the `n⊗n` off-diagonals would
  pass unnoticed. Here the normal takes every direction.
- **It is where the arithmetic law is worst.** A conducting disk in an
  insulating matrix sends its flux *normal* to the interface. The first-order
  prediction (dipole field, no fitted constants, `keff_disk_analytic.py`) is a
  +15% bias at `eps = L/100` and +1.9% at `L/800`.

## Running it

```bash
cd /Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM
./studies/keff_sharp_limit/disk/verify_keff_disk.sh
```

Four rungs (`eps = L/100 … L/800`, `Ny = 512 … 4096`, `eps/dy = 5.12` held
fixed as in the laminate ladder), each run under **both** laws: 8 `-keff_only`
solves, no time integration. `--rungs`, `--dry-run` as in the laminate driver.

The ladder starts at L/100 on purpose. At L/50, `eps/R = 0.08` and the visible
band (~9·eps) spans most of the radius, so the thin-interface expansion being
tested does not hold: the first run measured +46% for arith against +30% first
order, and +6.7% for tensor.

Writes `keff_disk.csv`, `keff_disk.log`, `keff_disk.png` here.

## The gates

Compared quantity: `(k_00 + k_11)/2` against the sharp value at the *nominal*
`f`. (On a curved interface the diffuse disk really holds `π³eps²/3` more ice;
that is part of the `O(eps²)` residual, not an IC error to divide out.)

| gate | law | pass condition |
|---|---|---|
| 1 | both | `k_00 = k_11`, `k_01 = k_10 = 0` to 1e−3 |
| 2 | arith | quadratic-fit intercept within 0.5% of Rayleigh; linear coefficient within 15% of the dipole prediction |
| 3 | tensor | error ≤ 0.2 × arith error at every rung; ≤ 0.5% at the finest |
| 4 | tensor | quadratic-fit intercept within 0.5% of Rayleigh; linear coefficient below 10% of arith's |

The dipole slope ignores the neighbouring disks' higher multipoles, so gate 2's
slope tolerance is loose on purpose. Gates 3 and 4 are the tensor-law test.

## Result (2026-09-23): all gates pass

| eps | arith err (pred. first order) | arith err (meas.) | tensor err (meas.) | tensor err, IC offset removed |
|---|---|---|---|---|
| L/100 | +15.07% | +18.46% | +1.377% | +1.219% |
| L/200 | +7.53% | +8.37% | +0.283% | +0.205% |
| L/400 | +3.77% | +3.99% | +0.073% | +0.034% |
| L/800 | +1.88% | +1.95% | +0.026% | +0.006% |

- **Arith:** fitted slope 440.9 vs predicted 445.5 W m⁻² K⁻¹ (−1.0%); fitted
  intercept within 0.04% of Rayleigh. The first-order theory is confirmed on a
  curved interface, with no fitted constants.
- **Tensor:** fitted first-order coefficient −2.5% of arith's; intercept within
  0.05% of Rayleigh. Error at most 0.075× arith's at every rung.
- **The IC offset.** `IC_COORD_UNIFORM` places the periodic mesh's N nodes at
  `L·i/(N+1)` instead of `L·i/N`, stretching the disk by `(N+1)/N`. Its measured
  area excess over the exact diffuse disk is `2/N` to four digits at every
  rung, which adds `(dk/df)·f·2/N` to both laws. With that removed, the tensor
  residual falls ~6× per halving of eps (observed order ≈ 2.5): the first-order
  term is gone. The packing ICs use `IC_COORD_GREVILLE` and are not affected.

