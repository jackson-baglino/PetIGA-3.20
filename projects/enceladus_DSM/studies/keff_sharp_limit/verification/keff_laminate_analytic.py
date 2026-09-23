"""Closed-form effective conductivity of a diffuse-interface laminate.

This module is the single source of truth for the analytic predictions used by
both the pass/fail gate (verify_keff_sharp_limit.sh) and the figures
(plot_keff_sharp_limit.py). Nothing here duplicates solver logic: it implements
the formulas derived in

    projects/effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex

sections 6 (the surface excess) and 7 (the laminate), which are what this study
exists to test.

THE ONE-PARAGRAPH VERSION
------------------------
A phase-field cell and a sharp cell use the SAME conductivity law
K(phi) = K_a + (K_i - K_a) phi; they differ only in the argument, the smoothed
phi^eps versus the sharp indicator. To first order in eps the diffuse band acts
as a sharp interface carrying two Gibbs surface excesses per unit area, a
tangential (conductance) one and a normal (resistance) one:

    Sigma_t = int [ K(phi^eps) - K_star ] ds        = 0        exactly
    Sigma_n = int [ 1/K(phi^eps) - 1/K_star ] ds    = -eps*C   C > 0

with C = (1/K_a - 1/K_i) ln(K_i/K_a). Sigma_t vanishes because the logistic
profile is antisymmetric, sigma(-u) = 1 - sigma(u), so the arithmetic law moves
no ice across the interface and the cell's ice fraction is exact at every eps.
Sigma_n does not vanish because 1/K is convex in phi: the band short-circuits
what should be a series resistance. Hence k_parallel is eps-exact and k_perp is
biased HIGH, first order in eps, with no free constants.

WHY RESISTIVITY, NOT CONDUCTIVITY
---------------------------------
The prediction is linear in eps in RESISTIVITY <1/K> and a hyperbola in
conductivity. Fit and plot 1/k_perp against eps. Done that way the test has two
independently predicted numbers and no fitted ones -- see resistivity_line().
Plotted as k_perp the same data is a curve whose agreement can only be eyeballed.

Verified against direct numerical quadrature of <1/K(phi^eps)> over the cell,
using the same nearest-interface signed-distance profile the solver builds in
FormInitialIceSlab2D: agreement is 5e-5 relative at eps = L/50 and machine
precision by eps = L/128 (the residual is tail overlap between neighbouring
interfaces, O(exp(-L/4eps))).

THE TENSOR LAW (-keff_interp tensor)
------------------------------------
Arithmetic along the interface, harmonic across it. On the laminate the normal
is exactly the layer normal, so k_perp = 1/<1/K_harm(phi)>, and 1/K_harm is
AFFINE in phi: <1/K_harm> = phi_bar/K_i + (1 - phi_bar)/K_a, the sharp value,
at every eps. The predicted ladder slope is therefore ZERO, and the whole
prediction is the flat line at the sharp intercept. Every function below takes
interp="arith" (default) or "tensor".
"""

import math

# Solver defaults, src/enceladus_main.c:58-59. Override per call if they change.
K_ICE_DEFAULT = 2.29
K_AIR_DEFAULT = 0.02


def excess_constant(k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """C = (1/K_a - 1/K_i) ln(K_i/K_a), the normal surface-excess coefficient.

    Sigma_n = -eps * C is the excess resistance per unit interface area, from
    Eq. (eq:sigman) of the note. Units m K / W. At K_i/K_a = 2.29/0.02 this is
    234.959, and it is the entire content of the predicted ladder slope.
    """
    return (1.0 / k_air - 1.0 / k_ice) * math.log(k_ice / k_air)


def k_parallel(phi_bar, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """Conductivity along the layers: the arithmetic (Voigt) mean.

    EXACT at every eps, because the tangential surface excess vanishes
    identically. A measured drift here is a real defect, not interface bias.

    Pass the MEASURED phi_bar from the k_eff CSV, not the nominal 0.5 -- see
    phi_bar_note().
    """
    return k_air + (k_ice - k_air) * phi_bar


def resistivity_sharp(phi_bar, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """<1/K> for the sharp laminate: the harmonic (Reuss) mean, inverted."""
    return phi_bar / k_ice + (1.0 - phi_bar) / k_air


INTERPS = ("arith", "tensor")


def _slope_factor(interp):
    """1 for the arithmetic law (Sigma_n = -eps*C), 0 for the tensor law."""
    if interp not in INTERPS:
        raise ValueError(f"interp must be one of {INTERPS}, got {interp!r}")
    return 1.0 if interp == "arith" else 0.0


def resistivity(eps, phi_bar, L, n_gamma=2,
                k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT, interp="arith"):
    """<1/K> for the DIFFUSE laminate at interface decay length eps.

    Eq. (eq:lamperp). n_gamma is the number of interfaces per period: 2 under
    -periodic 1, because the y=0 seam is a real interface and must be resolved
    like the one at y = Ly/2.

    Exact up to tail overlap between neighbouring interfaces, O(exp(-L/4eps)).
    """
    return (resistivity_sharp(phi_bar, k_ice, k_air)
            - _slope_factor(interp) * n_gamma * eps / L
            * excess_constant(k_ice, k_air))


def k_perp(eps, phi_bar, L, n_gamma=2,
           k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT, interp="arith"):
    """Conductivity across the layers. Biased HIGH, first order in eps, under
    the arithmetic law; eps-exact under the tensor law."""
    return 1.0 / resistivity(eps, phi_bar, L, n_gamma, k_ice, k_air, interp)


def k_perp_sharp(phi_bar, k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT):
    """The eps -> 0 limit of k_perp."""
    return 1.0 / resistivity_sharp(phi_bar, k_ice, k_air)


def resistivity_line(phi_bar, L, n_gamma=2,
                     k_ice=K_ICE_DEFAULT, k_air=K_AIR_DEFAULT, interp="arith"):
    """The two predicted numbers of the test, as (intercept, slope).

        1/k_perp(eps) = intercept - slope * eps

    intercept = <1/K>_sharp   and   slope = (n_gamma/L) * C.

    BOTH are predicted, neither is fitted. That is the whole point of working
    in resistivity: the gate compares a measured line against a known line
    rather than eyeballing a curve.
    """
    return (resistivity_sharp(phi_bar, k_ice, k_air),
            _slope_factor(interp) * n_gamma / L * excess_constant(k_ice, k_air))


def phi_bar_note():
    """Why the analytic forms are evaluated at the measured phi_bar."""
    return (
        "Evaluate every prediction at the phi_bar column of the k_eff CSV, not "
        "at the nominal 0.5. The discrete field is a spline fit to nodal values "
        "and its mean lands near but not on 0.5 (~4e-4 relative at Nx=Ny=256). "
        "Feeding 0.5 in would predict k_00 = 1.155000 where the correct "
        "discrete answer is 1.155454, and the gate would flag a 4e-4 'solver "
        "error' that is really the initial condition's quadrature error."
    )


if __name__ == "__main__":
    # Self-check: reproduce the two numbers stated independently in
    # inputs/geometry/iceslab/iceslab_2D_L1mm_eps20um_keff.opts:11-12.
    # Failing this means THIS module is wrong, not the solver.
    L = 1.0e-3
    expect = {"k_par": 1.155000, "k_perp": 0.0396538}
    got = {"k_par": k_parallel(0.5), "k_perp": k_perp_sharp(0.5)}
    # 1e-5, not tighter: the header quotes k_perp to seven places as 0.0396538
    # where the exact value is 0.03965372, so it is a last-digit rounding off by
    # 3e-6. The tolerance has to admit the header's own precision.
    TOL = 1e-5
    print("self-check against the opts-file header (phi_bar = 0.5, sharp):")
    ok = True
    for key, want in expect.items():
        rel = abs(got[key] - want) / want
        flag = "OK " if rel < TOL else "FAIL"
        ok &= rel < TOL
        print(f"  {flag} {key:8s} got {got[key]:.7f}  want {want:.7f}  "
              f"(rel {rel:.1e})")
    print(f"\nexcess constant C = {excess_constant():.4f} m K / W")
    icept, slope = resistivity_line(0.5, L)
    print(f"predicted line:  1/k_perp = {icept:.4f} - {slope:.4e} * eps\n")
    print(f"{'eps':>10}  {'<1/K>':>10}  {'k_perp':>10}  {'bias':>8}")
    for d in (50, 64, 128, 256, 512):
        e = L / d
        print(f"     L/{d:<4d}  {resistivity(e, 0.5, L):10.4f}  "
              f"{k_perp(e, 0.5, L):10.6f}  "
              f"{100*(k_perp(e, 0.5, L)/k_perp_sharp(0.5) - 1):+7.1f}%")
    raise SystemExit(0 if ok else 1)
