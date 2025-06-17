import numpy as np
import matplotlib.pyplot as plt
import warnings

# Constants
kB = 1.38e-23         # Boltzmann constant (J/K)
m_H2O = 3.0e-26       # Mass of water molecule (kg)
pi = np.pi
eps = 9.0e-7
rho_i = 917.0              # Density of ice (kg/m³)
rho_a = 1.341
k_a = 0.02
k_i = 2.29
cp_a = 1005.0          # Specific heat capacity of air (J/(kg·K))
cp_i = 1.96e3         # Specific heat capacity of ice (J/(kg·K))
diff_sub = 0.5 * (k_a / rho_a / cp_a + k_i / rho_i / cp_i)
dif_vap = 2.178e-5  


def sigma0(temp_c):
  sig = np.array([3.0e-3, 4.1e-3, 5.5e-3, 8.0e-3, 4.0e-3,
                  6.0e-3, 3.5e-2, 7.0e-2, 1.1e-1, 0.75])
  tem = np.array([-0.0001, -2.0, -4.0, -6.0, -7.0,
                  -10.0, -20.0, -30.0, -40.0, -100.0])

  if temp_c > tem[0] or temp_c < tem[-1]:
    warnings.warn("Temperature out of range in sigma0.")

  interv = -1
  for i in range(len(tem)):
    if temp_c <= tem[i]:
      interv = i

  if interv == -1:
    return sig[0]
  elif interv == len(tem) - 1:
    return sig[-1]
  else:
    t0, t1 = abs(tem[interv]), abs(tem[interv + 1])
    s0, s1 = sig[interv], sig[interv + 1]
    log_interp = np.log10(s0) + (np.log10(s1) - np.log10(s0)) / (np.log10(t1) - np.log10(t0)) * (np.log10(abs(temp_c)) - np.log10(t0))
    return 10 ** log_interp

import numpy as np

def rho_vs_I(temp_c):
  """
  Computes the saturation vapor density and optionally its derivative.
  
  Parameters
  ----------
  temp_c : float
      Temperature in Celsius.
  user : dict
      Must include:
        - 'rho_air': air density [kg/m³]
  
  Returns
  -------
  rho_vs : float
      Saturation vapor density [kg/m³]
  d_rhovs : float
      Derivative of rho_vs with respect to temperature [kg/m³/°C]
  """
  # Constants
  Patm = 101325.0       # Pa
  bb = 0.62             # Empirical constant
  temp_K = temp_c + 273.15

  # Empirical coefficients for saturation vapor pressure
  K0 = -0.5865e4
  K1 = 0.2224e2
  K2 = 0.1375e-1
  K3 = -0.3403e-4
  K4 = 0.2697e-7
  K5 = 0.6918

  # Saturation vapor pressure [Pa]
  Pvs = np.exp(K0 * temp_K**-1 + K1 + K2 * temp_K +
               K3 * temp_K**2 + K4 * temp_K**3 + K5 * np.log(temp_K))

  # Derivative of saturation vapor pressure w.r.t. temperature [Pa/K]
  Pvs_T = Pvs * (-K0 * temp_K**-2 + K2 + 2.0 * K3 * temp_K +
                 3.0 * K4 * temp_K**2 + K5 / temp_K)

  # Saturation vapor density [kg/m³]
  rho_vs = rho_a * bb * Pvs / (Patm - Pvs)

  # Derivative of vapor density w.r.t. temperature [kg/m³/°C]
  denom = (Patm - Pvs)**2
  d_rhovs = rho_a * bb * (Pvs_T * (Patm - Pvs) + Pvs * Pvs_T) / denom

  return rho_vs, d_rhovs

def compute_updated_params(temp_c, rhovs):
  T = temp_c + 273.15
  rho_rhovs = rho_i / rho_vs_I(temp_c)[0]

  alpha_k = 

  d0 = 2.548e-8 / temp_c
  beta0 = 1.0 / (alpha_k * v_kin)

  d0_sub = d0 / rho_rhovs
  beta_sub = beta0 / rho_rhovs

  a1 = 5.0
  a2 = 0.1581

  lambda_sub = 

  tau_sub = 

  mob_sub = 
  alph_sub = 

  return beta_sub, mob_sub, alph_sub

def main():
  temps = np.linspace(-80, -5, 300)
  betas, mobs, alphas = [], [], []

  for T in temps:
    beta, mob, alpha = compute_updated_params(T, rhovs=0.002, user=user_params)
    betas.append(beta)
    mobs.append(mob)
    alphas.append(alpha)

  fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)

  axes[0].plot(temps, betas, lw=2)
  axes[0].set_ylabel(r'$\beta_{\mathrm{sub}}$', fontsize=12)
  axes[0].grid(True)

  axes[1].plot(temps, mobs, lw=2)
  axes[1].set_ylabel('mob_sub', fontsize=12)
  axes[1].grid(True)

  axes[2].plot(temps, alphas, lw=2)
  axes[2].set_ylabel(r'$\alpha_{\mathrm{sub}}$', fontsize=12)
  axes[2].set_xlabel('Temperature (°C)', fontsize=12)
  axes[2].grid(True)

  fig.suptitle('Temperature Dependence of Updated Sublimation Parameters', fontsize=14)
  axes[2].invert_xaxis()
  fig.tight_layout(rect=[0, 0, 1, 0.96])
  plt.show()

if __name__ == "__main__":
  main()