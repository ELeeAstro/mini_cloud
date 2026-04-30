import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
V_seed = 4.0/3.0 * np.pi * r_seed**3

rho_d = [1.99, 3.9]
m_seed = V_seed * rho_d[0]

fname = 'tracers.txt'

data = np.loadtxt(fname,skiprows=3)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:8]
q_0 = data[:,10]
q_1 = data[:,11:13]
vf = data[:,15:20]

nlay = len(pl)

nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

N_c = np.zeros(nlay)
N_c[:] = q_0[:]*nd_atm[:]

rho = np.zeros(nlay)
rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

rho_c = np.zeros((nlay,2))
rho_c_t = np.zeros(nlay)
for j in range(2):
  rho_c[:,j] = q_1[:,j]*rho[:]
  rho_c_t[:] = rho_c_t[:] + rho_c[:,j]

m_c = np.zeros(nlay)
m_c[:] = np.maximum(rho_c_t[:]/N_c[:],m_seed)

rho_d_m = np.zeros(nlay)
for j in range(2):
  rho_d_m[:] = rho_d_m[:] + (rho_c[:,j])/rho_c_t[:] * rho_d[j]

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0),r_seed) * 1e4

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

col = sns.color_palette('colorblind')

p_rc = ax1.plot(r_c,pl,c=col[0],label=r'$r_{\rm c}$')
p_nc = ax2.plot(N_c,pl,c=col[1],label=r'$N_{\rm c}$',ls='dashed')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')

plt.gca().invert_yaxis()
yticks = [100,10,1,0.1,0.01,1e-3,]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$']
ax1.set_yticks(yticks,yticks_lab)

plt.ylim(300,3e-3)

ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)

ax1.set_xlabel(r'$r_{\rm c}$ [$\mu$m]',fontsize=16)
ax2.set_xlabel(r'$N_{\rm c}$ [cm$^{-3}$]',fontsize=16)
ax1.set_ylabel(r'$p_{\rm gas}$ [bar]',fontsize=16)

ax2.set_zorder(1)
lns = p_rc + p_nc
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='upper left')

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.show()
