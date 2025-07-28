import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

# --- Global Plot Style ---
plt.rcParams.update({
  "font.size": 16,
  "font.family": "sans-serif",
  "axes.titlesize": 20,
  "axes.titleweight": "bold",
  "axes.labelsize": 18,
  "axes.labelweight": "bold",
  "axes.edgecolor": "#333333",
  "xtick.labelsize": 15,
  "ytick.labelsize": 15,
  "legend.fontsize": 15,
  "figure.facecolor": "#f8fafc",
  "axes.facecolor": "#f3f4f6",
  "grid.color": "#d1d5db",
  "grid.linewidth": 1.2,
  "axes.grid": True,
  "axes.axisbelow": True,
  "figure.autolayout": True
})

# ------------------------
# Global Physical Constants
# ------------------------
GAMMA_IA = 0.109
RHO_I = 917.0
THCOND_ICE = 2.29
THCOND_AIR = 0.02
CP_ICE = 1.96e3
CP_AIR = 1.044e3
M_MOL = 3.0e-26
K_BOLTZ = 1.380649e-23
RHO_AIR = 1.341
PATM = 101325.0
BB = 0.62
SIGMA_SURF = 0.999
EPS = 5.00e-09
a1 = 5.0
a2 = 0.15812
K0 = -0.5865e4
K1 = 0.2224e2
K2 = 0.1375e-1
K3 = -0.3403e-4
K4 = 0.2697e-7
K5 = 0.6918
PLOT_DIR = "./outputs/plot_sigma0/"

def ensure_plot_dir():
  os.makedirs(PLOT_DIR, exist_ok=True)

def load_sigma0_data():
  temperature_C = -np.array([
    2.0662284513557885, 4.138849483297185, 5.323308309678149,
    6.203550141721592, 7.204643971446396, 8.575858332006952,
    10.339056175102106, 11.39213060745516, 12.477159925905502,
    13.504294635360388, 15.513046129354148, 20.806359097398456,
    31.253004787164496, 40.757870820038875
  ])
  sigma0 = np.array([
    0.004214474496253623, 0.0054980285416265295, 0.007521232999978934,
    0.007379749932627774, 0.004019064140763167, 0.004723107309565331,
    0.006220353575957829, 0.013045585034418642, 0.019438327725533045,
    0.015042465999957863, 0.020000000000000014, 0.0360329075195012,
    0.07004194652678876, 0.1298370424308925
  ])
  return temperature_C, sigma0

def interpolate_sigma0(temp_C, sigma0):
  temp_C = np.asarray(temp_C)
  sigma0 = np.asarray(sigma0)
  sorted_indices = np.argsort(temp_C)
  temp_sorted = temp_C[sorted_indices]
  sigma_sorted = sigma0[sorted_indices]
  return interp1d(temp_sorted, sigma_sorted, kind='linear', fill_value='extrapolate')

def compute_d0(temp_C):
  temp_C = np.asarray(temp_C)
  c_i = RHO_I / M_MOL
  d0 = GAMMA_IA / (c_i * K_BOLTZ * (temp_C + 273.15))
  print(f"Coefficient is: {GAMMA_IA / (c_i * K_BOLTZ)}")
  return d0

def compute_alpha_basal(temp_C, sigma0):
  temp_C = np.asarray(temp_C)
  sigma0 = np.asarray(sigma0)
  alpha_basal = np.exp(-sigma0 / SIGMA_SURF)
  return alpha_basal

def rhovs_I(temp_C):
  temp_C = np.asarray(temp_C)
  temp_K = temp_C + 273.15
  lnPvs = (K0 / temp_K + K1 + K2 * temp_K +
           K3 * temp_K**2 + K4 * temp_K**3 +
           K5 * np.log(temp_K))
  Pvs = np.exp(lnPvs)
  rho_vs = RHO_AIR * BB * Pvs / (PATM - Pvs)
  return rho_vs

def compute_beta_sub(temp_C, sigma0):
  temp_C = np.asarray(temp_C)
  sigma0 = np.asarray(sigma0)
  rhoI_vs = rhovs_I(temp_C)
  v_kin = (rhoI_vs / RHO_I) * np.sqrt(K_BOLTZ * (temp_C + 273.15) / (2 * np.pi * M_MOL))
  alpha_basal = compute_alpha_basal(temp_C, sigma0)
  beta_sub = 1 / (v_kin * alpha_basal)
  return beta_sub

def compute_lambda_sub(temp_C, sigma0):
  d0_sub = compute_d0(temp_C)
  rhoI_vs = rhovs_I(temp_C)
  rho_rhovs = RHO_I / rhoI_vs
  lambda_sub = a1 * (EPS / d0_sub) * rho_rhovs
  return lambda_sub

def compute_tau_sub(temp_C, sigma0):
  lambda_sub = compute_lambda_sub(temp_C, sigma0)
  beta_sub = compute_beta_sub(temp_C, sigma0)
  dif_vap = 2.178e-5
  diff_sub = 0.5 * (THCOND_AIR / RHO_AIR / CP_AIR + THCOND_ICE / RHO_I / CP_ICE)
  tau_sub = EPS * lambda_sub * (beta_sub / a1 + a2 * EPS / diff_sub + a2 * EPS / dif_vap)
  return tau_sub

def compute_mobility(temp_C, sigma0):
  tau_sub = compute_tau_sub(temp_C, sigma0)
  M0 = EPS / (3 * tau_sub)
  return M0

def compute_alpha_sub(temp_C, sigma0):
  lambda_sub = compute_lambda_sub(temp_C, sigma0)
  tau_sub = compute_tau_sub(temp_C, sigma0)
  alpha_sub = lambda_sub / tau_sub
  return alpha_sub

def pretty_axes(ax):
  """Apply final pretty touches to axes."""
  ax.spines['top'].set_visible(False)
  ax.spines['right'].set_visible(False)
  ax.spines['left'].set_linewidth(1.8)
  ax.spines['bottom'].set_linewidth(1.8)
  ax.xaxis.label.set_weight('bold')
  ax.yaxis.label.set_weight('bold')
  ax.title.set_weight('bold')
  ax.tick_params(axis='both', which='major', width=1.5, length=7)
  ax.tick_params(axis='both', which='minor', width=1, length=4)
  ax.grid(True, which='major', linestyle='-', linewidth=1.2, alpha=0.23)
  ax.grid(True, which='minor', linestyle=':', linewidth=0.7, alpha=0.18)
  ax.set_axisbelow(True)

def plot_sigma0(temp_C, sigma0_interp_func):
  ensure_plot_dir()
  temp_C = np.asarray(temp_C)
  temp_fine = np.linspace(min(temp_C), max(temp_C), 400)
  sigma_interp = sigma0_interp_func(temp_fine)
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(temp_C, sigma0_interp_func(temp_C), 'o', markersize=8, markerfacecolor='#2563eb',
          markeredgecolor='k', label='Data')
  ax.plot(temp_fine, sigma_interp, '-', linewidth=3.0, color='#0ea5e9', label='Linear fit')
  ax.set_xlabel("Temperature [°C]")
  ax.set_ylabel(r"$\sigma_0$")
  ax.set_title("Temperature Dependence of $\sigma_0$", pad=18)
  ax.legend(frameon=False)
  pretty_axes(ax)
  plt.tight_layout()
  plt.savefig(os.path.join(PLOT_DIR, "sigma0_vs_temperature.png"), dpi=300)
  plt.close()

def plot_d0(temp_C, d0_values):
  ensure_plot_dir()
  temp_C = np.asarray(temp_C)
  d0_values = np.asarray(d0_values)
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(temp_C, d0_values * 1e9, 's-', color='#059669', linewidth=2.5, markersize=7,
          markerfacecolor='#bbf7d0', markeredgewidth=1.4, markeredgecolor='#065f46')
  ax.set_xlabel("Temperature [°C]")
  ax.set_ylabel(r"$d_0$ [nm]")
  ax.set_title("Capillary Length $d_0$ vs Temperature", pad=18)
  pretty_axes(ax)
  plt.tight_layout()
  plt.savefig(os.path.join(PLOT_DIR, "d0_vs_temperature.png"), dpi=300)
  plt.close()

def plot_beta_sub(temp_C, beta_sub_values):
  ensure_plot_dir()
  temp_C = np.asarray(temp_C)
  beta_sub_values = np.asarray(beta_sub_values)
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(temp_C, beta_sub_values, '^-', color='#2563eb', linewidth=2.8, markersize=9,
          markerfacecolor='#dbeafe', markeredgewidth=1.2, markeredgecolor='#1e3a8a')
  ax.set_xlabel("Temperature [°C]")
  ax.set_ylabel(r"$\beta_{\rm sub}$")
  ax.set_title(r"Attachment Coefficient $\beta_{\rm sub}$ vs Temperature", pad=18)
  ax.set_yscale('log')
  pretty_axes(ax)
  plt.tight_layout()
  plt.savefig(os.path.join(PLOT_DIR, "beta_sub_vs_temperature.png"), dpi=300)
  plt.close()

def plot_mobility(temp_C, M0_values):
  ensure_plot_dir()
  temp_C = np.asarray(temp_C)
  M0_values = np.asarray(M0_values)
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(temp_C, M0_values, 'D-', color='#a21caf', linewidth=2.5, markersize=8,
          markerfacecolor='#f3e8ff', markeredgewidth=1.4, markeredgecolor='#6d28d9')
  ax.set_xlabel("Temperature [°C]")
  ax.set_ylabel(r"$M_0$")
  ax.set_title(r"Mobility $M_0$ vs Temperature", pad=18)
  ax.set_yscale('log')
  pretty_axes(ax)
  plt.tight_layout()
  plt.savefig(os.path.join(PLOT_DIR, "mobility_M0_vs_temperature.png"), dpi=300)
  plt.close()

def plot_alpha_sub(temp_C, alpha_sub_values):
  ensure_plot_dir()
  temp_C = np.asarray(temp_C)
  alpha_sub_values = np.asarray(alpha_sub_values)
  fig, ax = plt.subplots(figsize=(8, 5))
  ax.plot(temp_C, alpha_sub_values, 'o-', color='#eab308', linewidth=2.5, markersize=8,
          markerfacecolor='#fef9c3', markeredgewidth=1.4, markeredgecolor='#92400e')
  ax.set_xlabel("Temperature [°C]")
  ax.set_ylabel(r"$\alpha_{\rm sub}$")
  ax.set_title(r"Attachment Factor $\alpha_{\rm sub}$ vs Temperature", pad=18)
  ax.set_yscale('log')
  pretty_axes(ax)
  plt.tight_layout()
  plt.savefig(os.path.join(PLOT_DIR, "alpha_sub_vs_temperature.png"), dpi=300)
  plt.close()

def main():
  temp_C, sigma0 = load_sigma0_data()
  sigma0_func = interpolate_sigma0(temp_C, sigma0)
  d0_values = compute_d0(temp_C)
  beta_values = compute_beta_sub(temp_C, sigma0)
  M0_values = compute_mobility(temp_C, sigma0)
  alpha_sub_values = compute_alpha_sub(temp_C, sigma0)

  print("Sample values:")
  for T, s, d0, beta, M0, alpha_sub in zip(temp_C, sigma0, d0_values, beta_values, M0_values, alpha_sub_values):
    print(f"  T = {T:.2f}°C | σ₀ = {s:.3e} | d0 = {d0:.3e} | β = {beta:.3e} | M0 = {M0:.3e} | α_sub = {alpha_sub:.3e}")

  plot_sigma0(temp_C, sigma0_func)
  plot_d0(temp_C, d0_values)
  plot_beta_sub(temp_C, beta_values)
  plot_mobility(temp_C, M0_values)
  plot_alpha_sub(temp_C, alpha_sub_values)

  print(f"\nAll plots saved to: {os.path.abspath(PLOT_DIR)}")

if __name__ == "__main__":
  main()