import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns


Nc = 1.0
rho_d = 1.0


nr = 1000
r = np.logspace(-3,1,nr) * 1e-4

m = r**3*(4.0*np.pi*rho_d/3.0)

r_med = 1.0 * 1e-4
m_med = r_med**3*(4.0*np.pi*rho_d/3.0)

sig = [1.1,1.5,2.0,3.0]
sig = np.array(sig)
nsig = len(sig)

lnsig = np.log(sig)
lnsig2 = np.log(sig)**2

fx_m = np.zeros((nr,nsig))
fx_r = np.zeros((nr,nsig))

for i in range(nsig):
  fx_m[:,i] = Nc/(m*np.sqrt(2.0*np.pi)*lnsig[i]) * np.exp(-(np.log(m/m_med)**2)/(2.0*lnsig2[i]))
  fx_r[:,i] = Nc/(r*np.sqrt(2.0*np.pi)*lnsig[i]) * np.exp(-(np.log(r/r_med)**2)/(2.0*lnsig2[i]))


fig = plt.figure()

col = sns.color_palette('colorblind')

for i in range(nsig):
  plt.plot(r[:]*1e4,fx_m[:,i]*m,label=str(sig[i]),c=col[i])


plt.vlines(r_med*1e4, 1e-9, 1, colors='black',ls='dotted',label='Monodisperse')

plt.legend(title=r'$\sigma_{\rm g}$')

plt.ylim(1e-3,10)
plt.xlim(1e-2,10)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)
plt.tick_params(axis='both',which='major',labelsize=14)

plt.savefig('lognormal_comp_1.pdf',dpi=144,bbox_inches='tight')

fig = plt.figure()

col = sns.color_palette('colorblind')

for i in range(nsig):
  plt.plot(m[:],fx_m[:,i]*m,label=str(sig[i]),c=col[i])


plt.vlines(m_med, 1e-9, 1, colors='black',ls='dotted',label='Monodisperse')

plt.legend(title=r'$\sigma_{\rm g}$')

plt.ylim(1e-3,10)
#plt.xlim(1e-2,10)

plt.xscale('log')
plt.yscale('log')

plt.xlabel(r'$m$ [g]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)
plt.tick_params(axis='both',which='major',labelsize=14)

plt.savefig('lognormal_comp_2.pdf',dpi=144,bbox_inches='tight')

plt.show()