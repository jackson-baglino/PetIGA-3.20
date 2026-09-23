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
  +30% bias at `eps = L/50` and +3.8% at `L/400`.

## Running it

```bash
cd /Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM
./studies/keff_sharp_limit/disk/verify_keff_disk.sh
```

Four rungs (`eps = L/50 … L/400`, `Ny = 256 … 2048`, `eps/dy = 5.12` held
fixed as in the laminate ladder), each run under **both** laws: 8 `-keff_only`
solves, no time integration. `--rungs`, `--dry-run` as in the laminate driver.

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
