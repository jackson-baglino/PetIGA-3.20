#!/usr/bin/env python3
"""
Plot a double-well energy potential.

Default:  F(φ) = 1/4 * (φ^2 − 1)^2   (minima at φ = ±1, barrier at φ = 0)

Other common forms are included below—choose with `--form` at runtime.
Forms include: quartic_shifted, phi2_phi4, cahn_hilliard, double-well.

Double-well options:
- "quartic_shifted":   F = (1/4) * (phi**2 - 1)**2
- "phi2_phi4":         F = a*phi**2 + b*phi**4   (with a<0, b>0)
- "cahn_hilliard":     F = W * phi**2 * (1 - phi)**2  (minima at 0 and 1)
- "double-well":       Two-phase double-well with parameters sigma_i, sigma_a, Lambda
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

# --- Theming helper: light vs dark foreground on transparent background ---

def apply_theme(ax, theme: str = "light"):
    """Apply a light/dark foreground theme while keeping background transparent.
    Returns a tuple (fg_color, line_color).
    """
    fg = "black" if theme == "light" else "white"
    line = fg
    plt.rcParams.update({
        "text.color": fg,
        "axes.labelcolor": fg,
        "axes.edgecolor": fg,
        "xtick.color": fg,
        "ytick.color": fg,
    })
    # Update existing axes spines to match
    for sp in ax.spines.values():
        sp.set_color(fg)
    ax.title.set_color(fg)
    ax.xaxis.label.set_color(fg)
    ax.yaxis.label.set_color(fg)
    ax.tick_params(colors=fg)
    # Keep figure & axes backgrounds transparent explicitly
    ax.set_facecolor((0, 0, 0, 0))
    fig = ax.get_figure()
    fig.patch.set_alpha(0)
    return fg, line

def potential(phi, params, form="quartic_shifted"):
    """
    Double-well options:
    - "quartic_shifted":   F = (1/4) * (phi**2 - 1)**2
    - "phi2_phi4":         F = a*phi**2 + b*phi**4   (with a<0, b>0)
    - "cahn_hilliard":     F = W * phi**2 * (1 - phi)**2  (minima at 0 and 1)
    - "double-well":       Two-phase double-well with parameters sigma_i, sigma_a, Lambda
    """
    if form == "quartic_shifted":
        return 0.25 * (phi**2 - 1.0)**2
    elif form == "phi2_phi4":
        a = params.get("a", -1.0)
        b = params.get("b", 1.0)
        return a*phi**2 + b*phi**4
    elif form == "cahn_hilliard":
        W = params.get("W", 1.0)
        return W * phi**2 * (1.0 - phi)**2
    elif form == "double-well":
        # Two-phase double-well: 1/2*sigma_i*phi^2(1-phi)^2 + 1/2*sigma_a*(1-phi)^2*phi^2 + Lambda*phi^2(1-phi)^2
        phi_a = 1.0 - phi
        sigma_i = params.get("sigma_i", 10.0)
        sigma_a = params.get("sigma_a", 10.0)
        Lambda  = params.get("Lambda",  10.0)
        return 0.5 * sigma_i * phi**2 * (1.0 - phi)**2 \
             + 0.5 * sigma_a * phi_a**2 * (1.0 - phi_a)**2 \
             + Lambda * phi**2 * phi_a**2
    else:
        raise ValueError(f"Unknown form: {form}")

def main():
    p = argparse.ArgumentParser(description="Plot a double-well potential")
    p.add_argument("--form", type=str, default="quartic_shifted",
                   choices=["quartic_shifted", "phi2_phi4", "cahn_hilliard", "double-well"],
                   help="Functional form of the potential")
    p.add_argument("--a", type=float, default=-1.0, help="a for phi2_phi4")
    p.add_argument("--b", type=float, default= 1.0, help="b for phi2_phi4")
    p.add_argument("--W", type=float, default= 1.0, help="W for cahn_hilliard")
    p.add_argument("--sigma_i", type=float, default=1.0, help="sigma_i for double-well")
    p.add_argument("--sigma_a", type=float, default=1.0, help="sigma_a for double-well")
    p.add_argument("--Lambda",  type=float, default=0.0, help="Lambda coupling for double-well")
    p.add_argument("--phi_min", type=float, default=-1.5, help="x-range min")
    p.add_argument("--phi_max", type=float, default= 1.5, help="x-range max")
    p.add_argument("--n", type=int, default=1000, help="number of x samples")
    p.add_argument("--outfile", type=str, default="outputs/double_well.svg",
                   help="output image (e.g., .svg, .pdf, .png)")
    p.add_argument("--dpi", type=int, default=300, help="DPI (for raster formats)")
    p.add_argument("--theme", type=str, default="light", choices=["light", "dark"],
                   help="Foreground color theme for transparent output (light=black fg, dark=white fg)")
    args = p.parse_args()

    params = {"a": args.a, "b": args.b, "W": args.W,
              "sigma_i": args.sigma_i, "sigma_a": args.sigma_a, "Lambda": args.Lambda}

    # High-quality figure defaults
    plt.rcParams.update({
        "figure.figsize": (6.0, 4.0),
        "font.size": 14,
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    phi = np.linspace(args.phi_min, args.phi_max, args.n)
    F = potential(phi, params, form=args.form)

    fig, ax = plt.subplots()
    fg, line = apply_theme(ax, args.theme)
    ax.plot(phi, F, lw=2, color=line)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0.0, 2.0)

    # Labels/titles tailored to the selected form
    if args.form == "quartic_shifted":
        ax.set_title(r"$F(\phi) = \frac{1}{4}(\phi^2 - 1)^2$")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$F(\phi)$")
        # Annotate minima/barrier
        ax.scatter([-1, 1, 0], potential(np.array([-1, 1, 0.0]), {}, "quartic_shifted"),
                   zorder=3, s=25, color=fg)
        ax.annotate("min", (-1, 0), xytext=(-1.25, 0.1), color=fg,
                    arrowprops=dict(arrowstyle="->", color=fg), ha="right")
        ax.annotate("min", ( 1, 0), xytext=( 1.25, 0.1), color=fg,
                    arrowprops=dict(arrowstyle="->", color=fg), ha="left")
        ax.annotate("barrier", (0, 0.25), xytext=(0.3, 0.45), color=fg,
                    arrowprops=dict(arrowstyle="->", color=fg))

    elif args.form == "phi2_phi4":
        ax.set_title(r"$F(\phi) = a\,\phi^2 + b\,\phi^4$"
                     + f"  (a={args.a:g}, b={args.b:g})")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$F(\phi)$")

    elif args.form == "cahn_hilliard":
        ax.set_title(r"$F(\phi) = W\,\phi^2(1-\phi)^2$" + f"  (W={args.W:g})")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$F(\phi)$")

    elif args.form == "double-well":
        # ax.set_title(r"Double-Well")
        # ax.set_ylabel(r"$F(\phi)$")
        # Optional: mark wells at 0 and 1
        ax.scatter([0, 1], potential(np.array([0.0, 1.0]), params, "double-well"), s=20, zorder=3, color=fg)

    fig.tight_layout()
    fig.savefig(args.outfile, dpi=args.dpi, bbox_inches="tight", transparent=True)
    print(f"Saved: {args.outfile}")

if __name__ == "__main__":
    main()