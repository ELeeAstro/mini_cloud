import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns

nu = np.logspace(-2,2,1000)

print(nu)

nu_b = np.maximum(nu,0.334)

Kn_av_2 = nu_b**(1.0/3.0) * gamma(nu_b - 1.0/3.0) / gamma(nu_b)
# Mass-weighted and mass^2-weighted Knudsen numbers are already safe
Kn_av_m  = nu**(1.0/3.0) * gamma(nu + 2.0/3.0) / gamma(nu + 1.0)
Kn_av_m2 = nu**(1.0/3.0) * gamma(nu + 5.0/3.0) / gamma(nu + 2.0)

fig = plt.figure()


col = sns.color_palette('colorblind')

plt.plot(nu, Kn_av_2,c=col[0],label=r'$\overline{{\rm Kn}_{N}}$')

plt.plot(nu, Kn_av_m,c=col[1],label=r'$\overline{{\rm Kn}_{m}}$')
plt.plot(nu, Kn_av_m2,c=col[2],label=r'$\overline{{\rm Kn}_{m^{2}}}$')

plt.legend()

plt.hlines(1, 1e-2, 100, colors='black', linestyles='dotted')

plt.ylabel(r'Kn(gamma)/Kn(mono)',fontsize=16)
plt.xlabel(r'$\nu$',fontsize=16)

plt.ylim(1e-2,1000)
plt.xlim(1e-2,100)

plt.yscale('log')
plt.xscale('log')

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Kn_pop_average.pdf',dpi=144,bbox_inches='tight')

plt.show()



