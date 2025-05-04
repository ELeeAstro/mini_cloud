import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Full slip correction expression (Fuchs–Sutugin style) ---
def beta_full(Kn, A1=1.165, A2=0.483, A3=0.997):
    return 1 + Kn * (A1 + A2 * np.exp(-A3 / Kn))

# --- Log-linear model: log(β) ≈ log(1 + A·Kn) ---
def log_beta_linear(Kn, A):
    return np.log10(1 + A * Kn)

# --- Log-quadratic model: log(β) ≈ log(1 + A·Kn + B·Kn²) ---
def log_beta_quadratic(Kn, A, B):
    return np.log10(1 + A * Kn + B * Kn**2)

# --- Generate data over wide Kn range ---
Kn_vals = np.logspace(-3, 2, 1000)
beta_vals = beta_full(Kn_vals)
log_beta_vals = np.log10(beta_vals)

# --- Fit log-linear model ---
popt_lin, _ = curve_fit(log_beta_linear, Kn_vals, log_beta_vals)
A_lin = popt_lin[0]
log_fit_lin = log_beta_linear(Kn_vals, A_lin)
beta_fit_lin = 10.0**(log_fit_lin)
rel_error_lin = (beta_fit_lin - beta_vals) / beta_vals * 100

# --- Fit log-quadratic model ---
popt_quad, _ = curve_fit(log_beta_quadratic, Kn_vals, log_beta_vals)
A_quad, B_quad = popt_quad
log_fit_quad = log_beta_quadratic(Kn_vals, A_quad, B_quad)
beta_fit_quad = 10.0**(log_fit_quad)
rel_error_quad = (beta_fit_quad - beta_vals) / beta_vals * 100

# --- Plot the results ---
fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

axs[0].plot(Kn_vals, beta_vals, label='Full β(Kn)', linewidth=2)
axs[0].plot(Kn_vals, beta_fit_lin, '--', label=f'Log-Linear: 1 + {A_lin:.3f} Kn', linewidth=2)
axs[0].plot(Kn_vals, beta_fit_quad, '--', label=f'Log-Quadratic: 1 + {A_quad:.3f} Kn + {B_quad:.3f} Kn²', linewidth=2)
axs[0].set_ylabel("β(Kn)")
axs[0].set_title("Slip Correction: Full vs Log-Fitted Linear and Quadratic Models")
axs[0].legend()
axs[0].grid(True, which='both', linestyle=':')

axs[1].plot(Kn_vals, rel_error_lin, '--', label='Log-Linear Fit Error [%]', linewidth=2)
axs[1].plot(Kn_vals, rel_error_quad, '--', label='Log-Quadratic Fit Error [%]', linewidth=2)
axs[1].axhline(0, color='gray', linestyle='--')
axs[1].set_xscale('log')
axs[1].set_xlabel("Knudsen Number (Kn)")
axs[1].set_ylabel("Relative Error [%]")
axs[1].grid(True, which='both', linestyle=':')
axs[1].legend()

plt.tight_layout()
plt.show()

# --- Print fitted coefficients ---
print(f"Log-Linear fit coefficient: A = {A_lin:.6f}")
print(f"Log-Quadratic fit coefficients: A = {A_quad:.6f}, B = {B_quad:.6f}")
