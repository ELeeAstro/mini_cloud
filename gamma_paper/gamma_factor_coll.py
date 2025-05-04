import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln  # Use log-gamma instead of gamma
import seaborn as sns

# Define nu range
#nu = np.logspace(np.log10(1.0/2.0), 1, 100)
nu = np.logspace(-1, 1, 1000)



# Compute the terms using exponentials of log-gamma for Kn << 1 with clamping
Kn = 10.0
A = 1.165
B = 0.483
C = 0.997

nu_1 = np.maximum(0.667,nu)
log_gamma_nu_p1_3 = gammaln(nu_1 + 1.0/3.0)
log_gamma_nu_m1_3 = gammaln(nu_1 - 1.0/3.0)
log_gamma_nu_m2_3 = gammaln(nu_1 - 2.0/3.0)
log_gamma_nu = gammaln(nu_1)
B_l = (2.0 * (1.0 + np.exp(log_gamma_nu_p1_3 + log_gamma_nu_m1_3 - 2.0 * log_gamma_nu) + 
      1.639*Kn*nu_1**(1.0/3.0) 
      * (np.exp(log_gamma_nu_m1_3 - log_gamma_nu) 
      + np.exp(log_gamma_nu_p1_3 + log_gamma_nu_m2_3 - 2.0*log_gamma_nu)))) / (4.0*(1.0 + Kn * (A + B*np.exp(-C/Kn))))

nu_1 = np.maximum(0.334,nu)
log_gamma_nu_p1_3 = gammaln(nu_1 + 1.0/3.0)
log_gamma_nu_m1_3 = gammaln(nu_1 - 1.0/3.0)
log_gamma_nu_m2_3 = gammaln(nu_1 - 2.0/3.0)
log_gamma_nu = gammaln(nu_1)
B_l_2 = (2.0 * (1.0 + np.exp(log_gamma_nu_p1_3 + log_gamma_nu_m1_3 - 2.0 * log_gamma_nu)))/4.0


H = 1.0/np.sqrt(2.0)

nu_2 = np.maximum(0.501,nu)
log_gamma_nu = gammaln(nu_2)
log_gamma_nu_p1_3 = gammaln(nu_2 + 1.0/3.0)
log_gamma_nu_m1_3 = gammaln(nu_2 - 1.0/3.0)
log_gamma_nu_p2_3 = gammaln(nu_2 + 2.0/3.0)
log_gamma_nu_p1_6 = gammaln(nu_2 + 1.0/6.0)
log_gamma_nu_m1_2 = gammaln(nu_2 - 1.0/2.0)
log_gamma_nu_m1_6 = gammaln(nu_2 - 1.0/6.0)

B_h = (0.85 * np.sqrt(8.0) * nu_2**(-1.0/6.0) * 
       (np.exp(log_gamma_nu_p2_3 + log_gamma_nu_m1_2 - 2.0 * log_gamma_nu) + 
        2.0*np.exp(log_gamma_nu_p1_3 + log_gamma_nu_m1_6 - 2.0 * log_gamma_nu) + 
        np.exp(log_gamma_nu_p1_6 - log_gamma_nu))) / 8.0
B_h_2 = (H * np.sqrt(8.0) * nu_2**(-1.0/6.0) * 
       (np.exp(log_gamma_nu_p2_3 + log_gamma_nu_m1_2 - 2 * log_gamma_nu) + 
        2.0*np.exp(log_gamma_nu_p1_3 + log_gamma_nu_m1_6 - 2 * log_gamma_nu) + 
        np.exp(log_gamma_nu_p1_6 - log_gamma_nu))) / 8.0



log_gamma_nu = gammaln(nu)
log_gamma_nu_p2_3 = gammaln(nu + 2.0/3.0)  
log_gamma_nu_p1_3 = gammaln(nu + 1.0/3.0)        
grav = (nu**(-2.0/3.0) * 
        (np.exp(log_gamma_nu_p2_3 - log_gamma_nu) + 
         np.exp(2 * log_gamma_nu_p1_3 - 2.0 * log_gamma_nu))) / 2.0

# Plot results
fig = plt.figure()
col = sns.color_palette('colorblind')

#plt.plot(nu, B_l, label=r'Coag. (${\rm Kn}$ $\ll$ 1)', c=col[0])
plt.plot(nu, B_l_2, label=r'Coag. (${\rm Kn}$ $\ll$ 1)', c=col[0])
plt.plot(nu, B_h, label=r'Coag. (${\rm Kn}$ $\gg$ 1, $H$ = 0.85)', c=col[1])
plt.plot(nu, B_h_2, label=r'Coag. (${\rm Kn}$ $\gg$ 1, $H$ = 1/$\sqrt{2}$)', c=col[1], ls = 'dashed')
plt.plot(nu, grav, label=r'Coal.', c=col[2])

plt.hlines(1.0, np.log10(1.0/2.0), 10, colors='black', ls='dashed')

plt.vlines(1.0, 0.1, 200.0, colors='black', ls='dotted')
plt.text(1.1, 0.2, r'Exponential' "\n" r'distribution', c='black')

plt.legend()
plt.ylim(0.1, 200)
plt.xlim(0.1, 10)
plt.xscale('log')
plt.yscale('log')

# Add minor ticks every 0.02 on y-axis
#y_major_ticks = np.arange(0.7, 1.8, 0.1)  # Major ticks every 0.1
#y_minor_ticks = np.arange(0.7, 1.7, 0.025)  # Minor ticks every 0.02

#plt.yticks(y_major_ticks)  # Set major ticks
#plt.minorticks_on()  # Enable minor ticks
#plt.gca().set_yticks(y_minor_ticks, minor=True)  # Set minor ticks

plt.xlabel(r'$\nu$', fontsize=16)
plt.ylabel(r'$\frac{dN_{\rm c}}{dt}$(gamma/mono)', fontsize=16)

plt.tick_params(axis='both', which='major', labelsize=14)
plt.savefig('gamma_fact_coll.pdf', dpi=144, bbox_inches='tight')

plt.show()
