import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln  # Use log-gamma instead of gamma
import seaborn as sns

# Define sigma range
sig = np.linspace(1.0, 3.0, 100)
lnsig2 = np.log(sig)**2

rfac = np.exp(-0.5 * lnsig2)

H = 1.0/np.sqrt(2.0)
Kn = 0.01
A = 1.165
B = 0.483
C = 0.997

# Compute the terms using log-normal resuls
B_l = (2.0 * (1.0 + np.exp(1.0/9.0 * lnsig2) \
  + Kn*1.639*(np.exp(1.0/18.0 * lnsig2) + np.exp(5.0/18.0 * lnsig2)))) / (4.0*(1.0 + Kn * (A + B*np.exp(-C/Kn))))

B_h = (H * np.sqrt(8.0)*rfac**2*(np.exp(25.0/72.0 * lnsig2) + 2.0*np.exp(5.0/72.0 * lnsig2) +  np.exp(1.0/72.0 * lnsig2))) / 8.0

grav = (rfac**2*(np.exp(2.0/9.0 * lnsig2) + np.exp(1.0/9.0 * lnsig2))) / 2.0

# Plot results
fig = plt.figure()
col = sns.color_palette('colorblind')

plt.plot(sig, B_l, label=r'Coag. (${\rm Kn}$ $\ll$ 1)', c=col[0])
plt.plot(sig, B_h, label=r'Coag. (${\rm Kn}$ $\gg$ 1)', c=col[1])
plt.plot(sig, grav, label=r'Coal.', c=col[2])

plt.legend()
plt.ylim(0.5, 1.5)
plt.xlim(1.0, 3)

# Add minor ticks every 0.02 on y-axis
y_major_ticks = np.arange(0.7, 1.32, 0.1)  # Major ticks every 0.1
y_minor_ticks = np.arange(0.7, 1.32, 0.025)  # Minor ticks every 0.02

#plt.yticks(y_major_ticks)  # Set major ticks
#plt.minorticks_on()  # Enable minor ticks
plt.gca().set_yticks(y_minor_ticks, minor=True)  # Set minor ticks

plt.xlabel(r'$\sigma_{\rm g}$', fontsize=16)
plt.ylabel(r'$\frac{dN_{\rm c}}{dt}$(log-normal/mono)', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('lognormal_fact_coll.pdf', dpi=144, bbox_inches='tight')

plt.show()
