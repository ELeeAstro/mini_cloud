import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns


Nc = 1.0
rho_d = 1.0


nr = 1000
r = np.logspace(-3,1,nr) * 1e-4

m = r**3*(4.0*np.pi*rho_d/3.0)

rc = 1.0 * 1e-4
mc = rc**3*(4.0*np.pi*rho_d/3.0)

nu = [0.5,1,2,3,10]
nu = np.array(nu)
nnu = len(nu)

lam = mc/nu

fx = np.zeros((nr,nnu))

for i in range(nnu):
  fx[:,i] = Nc/(lam[i]**nu[i]*gamma(nu[i])) * m[:]**(nu[i]- 1) * np.exp(-m[:]/lam[i])

fx2 = np.zeros((nr))



# for i in range(nnu):
#   fx[:,i] = Nc/(lam[i]) * m[:]**(nu[i]- 1) * np.exp(-m[:]/lam[i])

fig = plt.figure()

col = sns.color_palette('colorblind')

for i in range(nnu):
  plt.plot(r[:]*1e4,fx[:,i]*m,label=str(nu[i]),c=col[i])

plt.vlines(1.0, 1e-9, 2, colors='black',ls='dotted')

plt.legend(title=r'$\nu$')

plt.ylim(1e-9,2)
plt.xlim(1e-3,10)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)
plt.tick_params(axis='both',which='major',labelsize=14)

plt.savefig('gamma_comp.pdf',dpi=144,bbox_inches='tight')

plt.show()