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

nu = 1.0

lam = mc/nu

fx = np.zeros((nr))

fx[:] = Nc/(lam**nu*gamma(nu) * m[:]**(nu- 1.0)) * np.exp(-m[:]/lam)

r_c = rc
r_n = r_c * nu**(-1.0/3.0) * gamma(nu+1.0/3.0)/gamma(nu)
r_p = r_c * ((nu + 1.0)/nu)**(1.0/3.0)

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(r[:]*1e4,fx[:]*m,label=str(nu),c=col[0])

plt.vlines(r_n*1e4,0,1,color=col[1],ls='dashed')
plt.vlines(r_c*1e4,0,1,color=col[2],ls='dashed')
plt.vlines(r_p*1e4,0,1,color=col[3],ls='dashed')

plt.legend(title=r'$\nu$')

plt.ylim(1e-3,1)
plt.xlim(1e-2,10)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)
plt.tick_params(axis='both',which='major',labelsize=14)

#plt.savefig('gamma_comp.pdf',dpi=144,bbox_inches='tight')

plt.show()