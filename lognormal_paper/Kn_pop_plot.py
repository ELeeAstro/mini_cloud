import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns

sig = np.linspace(1,3,1000)
lnsig2 = np.log(sig)**2

print(sig)

Kn_av_N = np.exp(1.0/18.0 * lnsig2)
Kn_av_m = np.exp(-5.0/18.0 * lnsig2)
Kn_av_m2 = np.exp(-11.0/18.0 * lnsig2)

fig = plt.figure()


col = sns.color_palette('colorblind')

plt.plot(sig, Kn_av_N,c=col[0],label=r'$\overline{{\rm Kn}_{N}}$')
plt.plot(sig, Kn_av_m,c=col[1],label=r'$\overline{{\rm Kn}_{m}}$')
plt.plot(sig, Kn_av_m2,c=col[2],label=r'$\overline{{\rm Kn}_{m^{2}}}$')

plt.legend()

#plt.hlines(1, 1e-2, 100, colors='black', linestyles='dotted')

plt.ylabel(r'Kn(lognormal)/Kn(mono)',fontsize=16)
plt.xlabel(r'$\sigma_{\rm g}$',fontsize=16)

plt.ylim(0.4,1.1)
#plt.xlim(1e-2,100)

#plt.yscale('log')
#plt.xscale('log')

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Kn_lognormal_pop_average.pdf',dpi=144,bbox_inches='tight')

plt.show()



