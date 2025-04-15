import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# --- Full expression for β(Kn) (Fuchs–Sutugin type) ---
def beta_full(Kn, A1=1.165, A2=0.483, A3=0.997):
    """Full slip correction expression based on empirical data"""
    return 1 + Kn * (A1 + A2 * np.exp(-A3 / Kn))

# --- Linear model to fit: β(Kn) ≈ 1 + A * Kn ---
def beta_linear(Kn, A):
    return 1 + A * Kn

# --- Fit β(Kn) to a linear approximation ---
def fit_beta_linear(kn_min=1e-4, kn_max=100, num_points=1000):
    Kn_vals = np.logspace(np.log10(kn_min), np.log10(kn_max), num_points)
    beta_vals = beta_full(Kn_vals)
    
    # Perform curve fitting
    popt, _ = curve_fit(beta_linear, Kn_vals, beta_vals)
    A_fit = popt[0]
    
    # Evaluate fit
    beta_fit = beta_linear(Kn_vals, A_fit)
    rel_error = (beta_fit - beta_vals) / beta_vals * 100

    # Plot results
    fig, axs = plt.subplots(2, 1, figsize=(9, 8), sharex=True)

    axs[0].plot(Kn_vals, beta_vals, label='Full β(Kn)', linewidth=2)
    axs[0].plot(Kn_vals, beta_fit, '--', label=f'Linear fit: 1 + {A_fit:.3f} Kn', linewidth=2)
    axs[0].set_ylabel("β(Kn)")
    axs[0].set_title("Slip Correction: Full vs Linear Fit")
    axs[0].legend()
    axs[0].grid(True, which='both', linestyle=':')

    axs[1].plot(Kn_vals, rel_error, '--', color='tab:red', label='Relative Error (%)', linewidth=2)
    axs[1].axhline(0, color='gray', linestyle='--')
    axs[1].set_xscale('log')
    axs[1].set_xlabel("Knudsen Number (Kn)")
    axs[1].set_ylabel("Relative Error [%]")
    axs[1].grid(True, which='both', linestyle=':')
    axs[1].legend()

    plt.tight_layout()
    plt.show()

    return A_fit

# --- Run the fit and report result ---
A_opt = fit_beta_linear()
print(f"Optimal linear coefficient A = {A_opt:.6f}")

print(1.165 + 0.483)
