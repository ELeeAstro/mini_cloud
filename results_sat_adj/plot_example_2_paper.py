import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99
sig = 2.0

fname = 'tracers.txt'

data = np.loadtxt(fname)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:7]
q_v = data[:,8]
q_c = data[:,9]
r_med = data[:,10] 
vf = data[:,11]

nlay = len(pl)

nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

rho = np.zeros(nlay)
rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

N_c = np.zeros(nlay)
N_c[:] = (3.0 * q_c[:] * rho[:])/(4.0*np.pi*rho_d*r_med[:]**3) \
  * np.exp(-9.0/2.0 * np.log(sig)**2)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

col = sns.color_palette('colorblind')

p_rc = ax1.plot(r_med*1e4,pl,c=col[0],label=r'$r_{\rm med}$')
p_nc = ax2.plot(N_c,pl,c=col[1],label=r'$N_{\rm c}$',ls='dashed')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')

plt.gca().invert_yaxis()
yticks = [100,10,1,0.1,0.01,1e-3,]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$']
ax1.set_yticks(yticks,yticks_lab)

#ax1.set_xlim(1e-3,1e1)
#ax2.set_xlim(1e-5,1e6)

plt.ylim(300,3e-3)

ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)

ax1.set_xlabel(r'$r_{\rm med}$ [$\mu$m]',fontsize=16)
ax2.set_xlabel(r'$N_{\rm c}$ [cm$^{-3}$]',fontsize=16)
ax1.set_ylabel(r'$p_{\rm gas}$ [bar]',fontsize=16)

# added these three lines
ax2.set_zorder(1)
lns = p_rc + p_nc
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='upper left')


plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

#plt.savefig('Y_425_mono.pdf',dpi=144,bbox_inches='tight')

plt.show()



