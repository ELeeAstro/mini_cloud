import numpy as np
import matplotlib.pylab as plt
from scipy.special import gammaln, gammaincc, gamma
import seaborn as sns

nnu = 1000
nu = np.logspace(-1,1,nnu)

rho = 1.99
rmin = 1e-7
xmin = 4.0/3.0 * np.pi * rmin**3 * rho

gi_m1_3 = np.zeros(nnu)
for i in range(nnu):
  if (nu[i] > 1.0/3.0):
    gi_m1_3[i] = gammaincc(nu[i] - 1.0/3.0, xmin) * gamma(nu[i] - 1.0/3.0)
  else:
    gi_m1_3[i] = ((gammaincc(nu[i]-1.0/3.0+1.0,xmin)*gamma(nu[i]-1.0/3.0+1.0) - xmin**(nu[i]-1.0/3.0)*np.exp(-xmin))) \
    /(nu[i]-1.0/3.0)

Kn_av_n_1 = nu**(1.0/3.0) * gamma(nu - 1.0/3.0) / gamma(nu)
Kn_av_n_2 = nu**(1.0/3.0) * gi_m1_3 / gamma(nu)
# Mass-weighted and mass^2-weighted Knudsen numbers are already safe
Kn_av_m  = nu**(1.0/3.0) * gamma(nu + 2.0/3.0) / gamma(nu + 1.0)
Kn_av_m2 = nu**(1.0/3.0) * gamma(nu + 5.0/3.0) / gamma(nu + 2.0)

fig = plt.figure()


col = sns.color_palette('colorblind')

plt.plot(nu, Kn_av_n_1,c=col[0],label=r'$\overline{{\rm Kn}_{N}}$ (gamma)')
plt.plot(nu, Kn_av_n_2,c=col[0],label=r'$\overline{{\rm Kn}_{N}}$ (inc. gamma)',ls='dashed')

plt.plot(nu, Kn_av_m,c=col[1],label=r'$\overline{{\rm Kn}_{m}}$')
plt.plot(nu, Kn_av_m2,c=col[2],label=r'$\overline{{\rm Kn}_{m^{2}}}$')

plt.legend()

plt.hlines(1, 1e-2, 100, colors='black', linestyles='dotted')

plt.ylabel(r'$\overline{\rm Kn}$(gamma)/$\overline{\rm Kn}$(mono)',fontsize=16)
plt.xlabel(r'$\nu$',fontsize=16)

plt.ylim(1e-1,1e4)
plt.xlim(1e-1,10)

plt.yscale('log')
plt.xscale('log')

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Kn_pop_average.pdf',dpi=144,bbox_inches='tight')

plt.show()



