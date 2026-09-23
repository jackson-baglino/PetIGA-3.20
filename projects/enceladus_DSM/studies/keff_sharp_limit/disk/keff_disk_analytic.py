"""Sharp reference and first-order bias for a square array of ice disks.

Single source of truth for gate_keff_disk.py. The theory is in

    projects/effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex

section 5 (the surface excesses) and its tensor-law corollary.

THE SHARP REFERENCE
-------------------
A square array of disks (conductivity K_i) in a matrix (K_a), area fraction f.
Rayleigh's multipole solution, in the form of Perrins, McKenzie & McPhedran
(1979), with T = -1/beta and beta = (K_i - K_a)/(K_i + K_a):

    k/K_a = 1 - 2f / ( T + f - 0.305827 f^4 T / (T^2 - 1.402958 f^8)
                         - 0.013362 f^8 )

Dropping everything past T + f leaves Maxwell-Garnett, k/K_a =
(1 + beta f)/(1 - beta f), which is exact to O(f^4). At the benchmark's
f = 0.196 the f^4 correction is 1.8e-4 relative, so the reference is good to
far better than any gate that uses it.

THE FIRST-ORDER BIAS OF THE ARITHMETIC LAW
------------------------------------------
To first order the band is a sharp interface carrying the excesses
Sigma_t = 0 and Sigma_n = -eps*C, C = (1/K_a - 1/K_i) ln(K_i/K_a). For a unit
mean gradient the effective conductivity then shifts by

    dk = -Sigma_n * (1/|Y|) INT_Gamma q_n^2 dGamma      (tangential term is 0)

where q_n is the normal flux of the SHARP solution on the interface. In the
dipole (Maxwell-Garnett) field the disk sees a local far field
E0 = 1/(1 - beta f), and the flux through its surface is
q_n = K_a (1 + beta) E0 cos(theta). Integrating cos^2 over the circle gives

    dk/deps = C * (pi R / L^2) * [K_a (1 + beta) E0]^2

No fitted constants. The dipole field neglects the neighbours' multipoles, the
same O(f^4) as above, so this slope is good to a few percent, not to 1e-4.

Under the tensor law both excesses vanish and the first-order slope is ZERO.
What remains is O(eps^2): among other terms, the diffuse disk holds
pi^3 eps^2 / 3 more ice than the sharp one (2% of f at eps = L/50).
"""

import math

K_ICE_DEFAULT = 2.29
K_AIR_DEFAULT = 0.02


def beta(k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    return (k_ice - k_air) / (k_ice + k_air)


def area_fraction(R, L):
    return math.pi * R * R / (L * L)


def k_maxwell_garnett(f, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    b = beta(k_ice, k_air)
    return k_air * (1.0 + b * f) / (1.0 - b * f)


def k_sharp(f, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """Rayleigh / Perrins et al. (1979), square array of disks."""
    T = -1.0 / beta(k_ice, k_air)
    denom = (T + f - 0.305827 * f**4 * T / (T * T - 1.402958 * f**8)
             - 0.013362 * f**8)
    return k_air * (1.0 - 2.0 * f / denom)


def excess_constant(k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    return (1.0 / k_air - 1.0 / k_ice) * math.log(k_ice / k_air)


def arith_slope(R, L, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """dk/deps of the arithmetic law, first order, dipole field. W m^-2 K^-1."""
    f = area_fraction(R, L)
    b = beta(k_ice, k_air)
    q = k_air * (1.0 + b) / (1.0 - b * f)
    return excess_constant(k_ice, k_air) * math.pi * R / (L * L) * q * q


if __name__ == "__main__":
    L, R = 1.0e-3, 2.5e-4
    f = area_fraction(R, L)
    mg, ks = k_maxwell_garnett(f), k_sharp(f)
    s = arith_slope(R, L)
    print(f"f = {f:.6f}   beta = {beta():.6f}")
    print(f"Maxwell-Garnett {mg:.7f}   Rayleigh/Perrins {ks:.7f}   "
          f"(rel {(ks - mg) / mg:.1e})")
    print(f"arith first-order slope dk/deps = {s:.2f} W/m^2/K\n")
    print(f"{'eps':>9}  {'arith pred':>11}  {'bias':>7}")
    for d in (100, 200, 400, 800):
        e = L / d
        print(f"   L/{d:<4d}  {ks + s * e:11.7f}  {s * e / ks:+7.2%}")
    # Self-check: the header of the geometry file quotes these two numbers.
    ok = abs(ks - 0.0295684) < 1e-7 and abs(mg - 0.0295632) < 1e-7
    print("\nself-check vs geometry-file header:", "OK" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)
