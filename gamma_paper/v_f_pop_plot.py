import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln  # Use log-gamma instead of gamma
import seaborn as sns

# Define nu range
nu = np.logspace(-1, 1, 100)


Kn = 0.01
A = 1.165
B = 0.483
C = 0.997

beta = 1.0 + Kn * (A + B*np.exp(-C/Kn))

# Compute log-Gamma values
log_gamma_nu = gammaln(nu)
log_gamma_nu_p1_3 = gammaln(nu + 1.0/3.0)
log_gamma_nu_p2_3 = gammaln(nu + 2.0/3.0)

vf_l_0 = (nu**(-2.0/3.0) * np.exp(log_gamma_nu_p2_3 - log_gamma_nu) + A*Kn*nu**(-1.0/3.0) * np.exp(log_gamma_nu_p1_3 - log_gamma_nu))/beta
vf_h_0 = nu**(-1.0/3.0) * np.exp(log_gamma_nu_p1_3 - log_gamma_nu)


log_gamma_nu_p1 = gammaln(nu + 1.0)
log_gamma_nu_p4_3 = gammaln(nu + 4.0/3.0)
log_gamma_nu_p5_3 = gammaln(nu + 5.0/3.0)

vf_l_1 = (nu**(-2.0/3.0) * np.exp(log_gamma_nu_p5_3 - log_gamma_nu_p1) + A*Kn*nu**(-1.0/3.0) * np.exp(log_gamma_nu_p4_3 - log_gamma_nu_p1))/beta
vf_h_1 = nu**(-1.0/3.0) * np.exp(log_gamma_nu_p4_3 - log_gamma_nu_p1)

log_gamma_nu_p2 = gammaln(nu + 2.0)
log_gamma_nu_p7_3 = gammaln(nu + 7.0/3.0)
log_gamma_nu_p8_3 = gammaln(nu + 8.0/3.0)

vf_l_2 = (nu**(-2.0/3.0) * np.exp(log_gamma_nu_p8_3 - log_gamma_nu_p2) + A*Kn*nu**(-1.0/3.0) * np.exp(log_gamma_nu_p7_3 - log_gamma_nu_p2))/beta
vf_h_2 = nu**(-1.0/3.0) * np.exp(log_gamma_nu_p7_3 - log_gamma_nu_p2)

# Plot results
fig = plt.figure()
col = sns.color_palette('colorblind')

plt.plot(nu, vf_l_0, label=r'${\rm Kn}$ $\ll$ 1', c=col[0])
plt.plot(nu, vf_h_0, label=r'${\rm Kn}$ $\gg$ 1', c=col[1])

plt.plot(nu, vf_l_1, c=col[0],ls='dashed')
plt.plot(nu, vf_h_1, c=col[1],ls='dashed')

plt.plot(nu, vf_l_2, c=col[0],ls='dashdot')
plt.plot(nu, vf_h_2, c=col[1],ls='dashdot')

plt.vlines(1.0, 0.1, 10, colors='black', ls='dotted')
plt.text(1.1, 0.2, r'Exponential' "\n" r'distribution', c='black')

plt.hlines(1.0, 0.1, 10, colors='black', ls='dotted')


plt.legend()
plt.ylim(0.1, 10)
plt.xlim(1e-1, 10)
plt.xscale('log')
plt.yscale('log')

# Add minor ticks every 0.02 on y-axis
#y_major_ticks = np.arange(0.4, 1.02, 0.1)  # Major ticks every 0.1
#y_minor_ticks = np.arange(0.4, 1.025, 0.025)  # Minor ticks every 0.02

#plt.yticks(y_major_ticks)  # Set major ticks
#plt.minorticks_on()  # Enable minor ticks
#plt.gca().set_yticks(y_minor_ticks, minor=True)  # Set minor ticks

plt.xlabel(r'$\nu$', fontsize=16)
plt.ylabel(r'$\overline{v_{\rm f}}$ (gamma/mono)', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.tick_params(axis='both', which='minor', labelsize=10, length=4)  # Customize minor ticks

plt.savefig('gamma_fact_vf.pdf', dpi=144, bbox_inches='tight')

plt.show()