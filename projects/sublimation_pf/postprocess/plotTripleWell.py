#!/usr/bin/env python3
"""
3D surface plot of a triple-well energy potential:

    F^tri = 1/2 φ_i^2 (1 - φ_i)^2
          + 1/2 φ_s^2 (1 - φ_s)^2
          + 1/2 φ_a^2 (1 - φ_a)^2
          + Λ φ_i^2 φ_s^2 φ_a^2

Under the constraint φ_a = 1 - φ_i, we visualize
F as a function of (φ_i, φ_s).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# -----------------------------
# Triple-well potential
# -----------------------------
def triple_well(phi_i, phi_s, phi_a, Lambda=1.0):
    """
    Triple-well free energy density F^tri(φ_i, φ_s, φ_a).

    Parameters
    ----------
    phi_i, phi_s, phi_a : array-like or float
        Phase fields (0..1).
    Lambda : float
        Coupling coefficient Λ.

    Returns
    -------
    F : array-like or float
        Free energy density.
    """
    term_i = 0.5 * phi_i**2 * (1.0 - phi_i)**2
    term_s = 0.5 * phi_s**2 * (1.0 - phi_s)**2
    term_a = 0.5 * phi_a**2 * (1.0 - phi_a)**2
    coupling = Lambda * phi_i**2 * phi_s**2 * phi_a**2

    return term_i + term_s + term_a + coupling


def main():
    # -----------------------------
    # Parameters
    # -----------------------------
    Lambda = 1.0          # coupling strength Λ
    npts   = 200          # resolution in each direction
    phi_min, phi_max = 0.0, 1.0


    # -----------------------------
    # Build grid in (φ_i, φ_s)
    # -----------------------------
    phi_i_vals = np.linspace(phi_min, phi_max, npts)
    phi_s_vals = np.linspace(phi_min, phi_max, npts)

    PHI_I, PHI_S = np.meshgrid(phi_i_vals, phi_s_vals, indexing="xy")
    PHI_A = 1.0 - PHI_I

    F_vals = triple_well(PHI_I, PHI_S, PHI_A, Lambda=Lambda)

    # -----------------------------
    # 3D surface plot
    # -----------------------------
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    # Use a stride to keep rendering smooth but not too heavy
    surf = ax.plot_surface(
        PHI_I,
        PHI_S,
        F_vals,
        rstride=3,
        cstride=3,
        cmap=cm.inferno,
        linewidth=0,
        antialiased=True,
        alpha=0.95,
    )

    # Labels and title
    ax.set_xlabel(r"$\phi_i$", fontsize=14, labelpad=10)
    ax.set_ylabel(r"$\phi_s$", fontsize=14, labelpad=10)
    ax.set_zlabel(r"$F^{\mathrm{tri}}$", fontsize=14, labelpad=10)

    title = (
        rf"Triple-well energy $F^{{\mathrm{{tri}}}}(\phi_i,\phi_s,\phi_a)$" + "\n"
        + rf"$\Lambda = {Lambda:.2g},\ \phi_a = 1 - \phi_i$"
    )
    ax.set_title(title, fontsize=14, pad=10)

    # Colorbar for F
    cbar = fig.colorbar(surf, shrink=0.7, aspect=20, pad=0.08)
    cbar.set_label(r"$F^{\mathrm{tri}}$", fontsize=12)

    # Make it a bit prettier
    ax.view_init(elev=30, azim=-135)  # tweak viewing angle if you like
    ax.tick_params(axis="both", which="major", labelsize=10)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()