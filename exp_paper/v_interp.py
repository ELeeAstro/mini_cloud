import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
R_gas = 8.31446261815324e7
amu = 1.66053906892e-24

T = 500.0
p = 1.0e-3 * 1e6

mu = 2.33
rho_d = 1.99
grav = 10.0**3.25

rho = (p*mu*amu)/(kb * T)

molg_g = 2.01588
LJ_g = 59.7
d_g = 2.827e-8

eta = (5.0/16.0) * (np.sqrt(np.pi*(molg_g*amu)*kb*T)/(np.pi*d_g**2)) \
  * ((((kb*T)/LJ_g)**(0.16))/1.22)

cT = np.sqrt((2.0 * kb * T) / (mu * amu))

mfp = (2.0*eta/rho) * np.sqrt((np.pi * mu)/(8.0*R_gas*T))


nr = 1000
r_c = np.logspace(-4,2,nr) * 1e-4

Kn = np.zeros(nr)
Kn[:] = mfp/r_c[:]

beta = np.zeros(nr)
beta[:] = 1.0 + Kn*(1.165 + 0.483 * np.exp(-0.997/Kn[:]))

vf_s = np.zeros(nr)
vf_s[:] = (2.0 * beta[:] * grav * r_c[:]**2 * (rho_d - rho))/(9.0 * eta) \
* (1.0 + ((0.45*grav*r_c[:]**3*rho*rho_d)/(54.0*eta**2))**(0.4))**(-1.25)

vf_e = np.zeros(nr)
vf_e[:] = (np.sqrt(np.pi)*grav*rho_d*r_c[:])/(2.0*cT*rho)

fx = np.zeros(nr)
fx[:] = 0.5 * (1.0 - np.tanh(2.0*np.log10(Kn[:])))

vf = np.zeros(nr)
vf[:] = fx[:]*vf_s[:] + (1.0 - fx[:])*vf_e[:]

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(Kn, vf_s, label='Stokes',c=col[0])
plt.plot(Kn, vf_e, label='Epstein',c=col[1])
plt.plot(Kn, vf, ls='dashed',c='black', label='interpolation')

plt.vlines(1, 1e-1, 1e5, colors='black', linestyles='dotted')

plt.legend()

plt.xlim(1e-2,1e2)
plt.ylim(1e-1,1e5)

plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$\overline{{\rm Kn}}$',fontsize=16)
plt.ylabel(r'$v_{\rm f}$ [cm s$^{-1}$]',fontsize=16)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('v_f_interp.pdf',dpi=144,bbox_inches='tight')


plt.show()

