import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99

dirs = ['../results_2_mono/','../results_3_lognormal/']
ndir = len(dirs)

fname = 'tracers_425.txt'

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

col = sns.color_palette('colorblind')

lss = ['dashed','solid']

for i in range(ndir):

  data = np.loadtxt(dirs[i]+fname)

  ni = data[:,1]

  Tl = data[:,2]
  pl = data[:,3]/1e5
  grav = data[:,4]
  mu = data[:,5]
  VMR = data[:,6:7]
  q_v = data[:,8]
  q_0 = data[:,9]
  q_1 = data[:,10]
  vf = data[:,11]

  nlay = len(pl)

  if (i == 0):
    nd = np.zeros((nlay,2))
    r_c = np.zeros((nlay,2))
    

  nd_atm = np.zeros(nlay)
  nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

  nd[:,i] = q_0[:]*nd_atm[:]

  rho = np.zeros(nlay)
  rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

  m_c = np.zeros(nlay)
  m_c[:] = (q_1[:]*rho[:])/(q_0[:]*nd_atm[:])

  r_c[:,i] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d))**(1.0/3.0),r_seed) * 1e4


nd_rel = np.zeros(nlay)
nd_rel[:] = (nd[:,1]-nd[:,0])/nd[:,0] * 100.0
r_c_rel = np.zeros(nlay)
r_c_rel[:] = (r_c[:,1]-r_c[:,0])/r_c[:,0] * 100.0

p_rc = ax1.plot(r_c_rel,pl,c=col[0],label=r'relative $r_{\rm c}$',ls=lss[i])
p_nc = ax2.plot(nd_rel,pl,c=col[1],label=r'relative $N_{\rm c}$',ls=lss[i])

plt.vlines(0, 3e-3, 300, colors='black', linestyles='dotted')

ax1.set_yscale('log')
#ax1.set_xscale('log')
#ax2.set_xscale('log')

plt.gca().invert_yaxis()
yticks = [100,10,1,0.1,0.01,1e-3,]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$']
ax1.set_yticks(yticks,yticks_lab)

#ax1.set_xlim(-10,10)
#ax2.set_xlim(-100,100)
#ax1.set_xticks(np.arange(-10, 12, step=2))
#ax2.set_xticks(np.arange(-100, 120, step=20))
#ax1.set_xticks(np.arange(-10, 11, step=1),minor=True)
#ax2.set_xticks(np.arange(-100, 110, step=10),minor=True)

ax1.set_xlim(-100,100)
ax2.set_xlim(-500,500)
ax1.set_xticks(np.arange(-100, 120, step=20))
ax2.set_xticks(np.arange(-500, 600, step=100))
ax1.set_xticks(np.arange(-100, 110, step=10),minor=True)
ax2.set_xticks(np.arange(-500, 510, step=10),minor=True)


plt.ylim(300,3e-3)

ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)

ax1.set_xlabel(r'$\delta$$r_{\rm c}$ [$\%$]',fontsize=16)
ax2.set_xlabel(r'$\delta$$N_{\rm c}$ [$\%$]',fontsize=16)
ax1.set_ylabel(r'$p_{\rm gas}$ [bar]',fontsize=16)

# added these three lines
ax2.set_zorder(1)
lns = p_rc + p_nc
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='lower right')


plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Y_425_mono_lognormal_rel.pdf',dpi=144,bbox_inches='tight')

plt.show()



