# condensation_coefficient.py
# Calculates the condensation coefficient (alpha) using the formulation from
# Kaempfer and Plapp (2009), "Phase changes and related mass fluxes in dry snow metamorphism"
# Inputs: Temperature (T, in Celsius), vapor density at saturation (rho_vs, kg/m^3),
# ice density (rho_i, kg/m^3), and kinetic coefficient (beta, m/s).

def condensation_coefficient(T, rhoI_rhoVS, beta):
    """
    Calculate the condensation coefficient (alpha) as per Kaempfer and Plapp (2009).
    Parameters:
        T (float): Temperature in Celsius
        rhoI_rhoVS (float): Ratio of ice density to saturation vapor density (kg/m^3)
        beta (float): Kinetic coefficient (m/s)
    Returns:
        alpha (float): Condensation coefficient (dimensionless)
    """
    # Convert temperature to Kelvin
    T_K = T + 273.15

    # Constants
    kB = 1.380649e-23  # Boltzmann constant (J/K)
    m_water = 18.01528e-3 / 6.02214076e23  # Mass of a water molecule (kg)

    # Calculate the mean thermal speed of water vapor molecules (m/s)
    import math
    v_kin = math.sqrt(kB * T_K / (2 * math.pi * m_water))

    # Condensation coefficient (alpha)
    alpha = (1 / beta) * rhoI_rhoVS * v_kin

    return alpha

# Example usage:
if __name__ == "__main__":
    # Example values
    T = -20.0  # Celsius
    rhoI_rhoVS = 1.820794 * 10 ** 5
    beta = 7.688954e-1    # m/s (example value, replace with actual)

    alpha = condensation_coefficient(T, rhoI_rhoVS, beta)
    print(f"Condensation coefficient (alpha): {alpha:.4f}")