import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln  # Use log-gamma instead of gamma
import seaborn as sns

# Define nu range
sig = np.linspace(1.0, 3.0, 100)
lnsig2 = np.log(sig)**2

# Compute the terms using log-normal expressions
v_l = np.exp(1/18 * lnsig2)*(np.exp(-0.5 * lnsig2))**1.0
v_h = np.exp(2/9 * lnsig2)*(np.exp(-0.5 * lnsig2))**2.0

# Plot results
fig = plt.figure()
col = sns.color_palette('colorblind')

plt.plot(sig, v_l, label=r'Cond. (${\rm Kn}$ $\ll$ 1)', c=col[0])
plt.plot(sig, v_h, label=r'Cond. (${\rm Kn}$ $\gg$ 1)', c=col[1])


plt.legend()
#plt.ylim(0.4, 1)
plt.xlim(1, 3)
#plt.xscale('log')

# Add minor ticks every 0.02 on y-axis
y_major_ticks = np.arange(0.4, 1.02, 0.1)  # Major ticks every 0.1
y_minor_ticks = np.arange(0.4, 1.025, 0.025)  # Minor ticks every 0.02

#plt.yticks(y_major_ticks)  # Set major ticks
#plt.minorticks_on()  # Enable minor ticks
plt.gca().set_yticks(y_minor_ticks, minor=True)  # Set minor ticks

plt.xlabel(r'$\sigma_{\rm g}$', fontsize=16)
plt.ylabel(r'$\frac{dm}{dt}$(log-normal/mono)', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=10, length=4)  # Customize minor ticks

plt.savefig('lognormal_fact_cond.pdf', dpi=144, bbox_inches='tight')

plt.show()
