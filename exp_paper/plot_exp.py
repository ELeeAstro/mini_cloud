import numpy as np
import matplotlib.pylab as plt
import seaborn as sns


nr = 1000
r = np.logspace(-3,1,nr)*1e-4
rho_d = 1.99

m = r**3 * ((4.0*np.pi*rho_d)/3.0)

r_c = [1e-3,1e-2,1e-1,1]
n = 4
r_c = np.array(r_c) * 1e-4

m_c = r_c**3 * ((4.0*np.pi*rho_d)/3.0)

Nc = 1.0

fr = np.zeros((n,nr))
for i in range(n):
  fr[i,:] = Nc/r_c[i] * np.exp(-r[:]/r_c[i])

fr_m = np.zeros((n,nr))
for i in range(n):
  fr_m[i,:] = Nc/m_c[i] * np.exp(-m[:]/m_c[i])

col = sns.color_palette('colorblind')

fig = plt.figure()

for i in range(n):
  plt.plot(r[:]*1e4,fr[i,:]*r[:],c=col[i])
  plt.vlines(r_c[i]*1e4,0,1.0,ls='dashed',color=col[i])

plt.xscale('log')

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$r$ $\cdot$ $f(r)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('exp_v_mono_1.pdf',dpi=144,bbox_inches='tight')

fig = plt.figure()

for i in range(n):
  plt.plot(m[:],fr_m[i,:]*m[:],c=col[i])
  plt.vlines(m_c[i],0,1.0,ls='dashed',color=col[i])

plt.xscale('log')

plt.xlim(7e-21,1e-10)

plt.xlabel(r'$m$ [g]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('exp_v_mono_2.pdf',dpi=144,bbox_inches='tight')

plt.show()

