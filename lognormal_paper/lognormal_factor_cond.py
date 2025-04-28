import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln  # Use log-gamma instead of gamma
import seaborn as sns

# Define nu range
nu = np.logspace(-1, 1, 100)

# Compute log-Gamma values
log_gamma_nu = gammaln(nu)
log_gamma_nu_p1_3 = gammaln(nu + 1.0/3.0)
log_gamma_nu_p2_3 = gammaln(nu + 2.0/3.0)

# Compute the terms using exponentials of log-Gamma ratios
v_l = nu**(-1.0/3.0) * np.exp(log_gamma_nu_p1_3 - log_gamma_nu)
v_h = nu**(-2.0/3.0) * np.exp(log_gamma_nu_p2_3 - log_gamma_nu)

# Plot results
fig = plt.figure()
col = sns.color_palette('colorblind')

plt.plot(nu, v_l, label=r'Cond. (${\rm Kn}$ $\ll$ 1)', c=col[0])
plt.plot(nu, v_h, label=r'Cond. (${\rm Kn}$ $\gg$ 1)', c=col[1])

plt.vlines(1.0, 0.4, 1, colors='black', ls='dotted')
plt.text(1.1, 0.5, r'Exponential' "\n" r'distribution', c='black')

plt.legend()
plt.ylim(0.4, 1)
plt.xlim(1e-1, 10)
plt.xscale('log')

# Add minor ticks every 0.02 on y-axis
y_major_ticks = np.arange(0.4, 1.02, 0.1)  # Major ticks every 0.1
y_minor_ticks = np.arange(0.4, 1.025, 0.025)  # Minor ticks every 0.02

plt.yticks(y_major_ticks)  # Set major ticks
plt.minorticks_on()  # Enable minor ticks
plt.gca().set_yticks(y_minor_ticks, minor=True)  # Set minor ticks

plt.xlabel(r'$\nu$', fontsize=16)
plt.ylabel(r'$\frac{dm}{dt}$(gamma/mono)', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=10, length=4)  # Customize minor ticks

plt.savefig('gamma_fact_cond.pdf', dpi=144, bbox_inches='tight')

plt.show()
