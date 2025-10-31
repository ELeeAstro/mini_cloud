import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99
V_seed = 4.0/3.0 * np.pi * r_seed**3
m_seed = V_seed * rho_d

dirs = ['../results_2_mono/','../results_2_exp/','../results_3_gamma/']
ndir = len(dirs)

fname = 'tracers_325.txt'

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()

col = sns.color_palette('colorblind')

lss = ['dashed','dotted','solid']

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
  q_2 = data[:,11]

  nlay = len(pl)

  nd_atm = np.zeros(nlay)
  nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

  nd = np.zeros(nlay)
  nd[:] = q_0[:]*nd_atm[:]

  rho = np.zeros(nlay)
  rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

  m_c = np.zeros(nlay)
  m_c[:] = np.maximum((q_1[:]*rho[:])/(q_0[:]*nd_atm[:]),m_seed)

  r_c = np.zeros(nlay)
  r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d))**(1.0/3.0),r_seed) * 1e4

  sig2 = np.zeros(nlay)
  sig2[:] = np.maximum((q_2[:]*rho[:]**2)/(q_0[:]*nd_atm[:]) - ((q_1[:]*rho[:])/(q_0[:]*nd_atm[:]))**2,m_seed**2)

  nu = np.zeros(nlay)
  nu[:] = m_c[:]**2/sig2[:]
  nu[:] = np.maximum(m_c[:]**2/sig2[:],0.001)
  nu[:] = np.minimum(nu[:],100.0)

  lam = np.zeros(nlay)
  lam[:] = m_c[:]/nu[:]
  m_p = np.zeros(nlay)
  #m_p[:] = np.maximum((q_2[:]*rho[:]**2)/(q_1[:]*rho[:]),m_seed)
  m_p[:] = np.maximum(lam[:] * (nu[:] + 1.0),m_seed)

  r_p = np.zeros(nlay)
  r_p[:] = np.maximum(((3.0*m_p[:])/(4.0*np.pi*rho_d))**(1.0/3.0),r_seed) * 1e4
  #r_p[:] = np.maximum(r_c[:]/1e4 * ((nu[:] + 1.0)/nu[:])**(1.0/3.0),r_seed) * 1e4
  
  p_rc = ax1.plot(r_c,pl,c=col[0],label=r'$r_{\rm c}$',ls=lss[i])
  p_nc = ax2.plot(nd,pl,c=col[1],label=r'$N_{\rm c}$',ls=lss[i])

  if (i == 2):
    p_rp = ax1.plot(r_p,pl,c=col[2],label=r'$r_{\rm p}$',ls=lss[i])


# Read CARMA values
fname = 'Carma_325.txt'
data = np.loadtxt(fname,delimiter=',')
r_Ca = data[:,0]
p_Ca = data[:,1]

p_Ca = ax1.plot(r_Ca, p_Ca, c=col[3],ls='dotted',label=r'$r_{\rm eff}$ (Gao+ 2018)')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax2.set_xscale('log')



plt.gca().invert_yaxis()
yticks = [100,10,1,0.1,0.01,1e-3,]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$']
ax1.set_yticks(yticks,yticks_lab)

ax1.set_xlim(1e-3,1e2)


ax2.set_xlim(1e-4,1e-1)
#ax2.set_xlim(1e-1,1e2)

plt.ylim(300,3e-3)

ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)

ax1.set_xlabel(r'$r$ [$\mu$m]',fontsize=16)
ax2.set_xlabel(r'$N_{\rm c}$ [cm$^{-3}$]',fontsize=16)
ax1.set_ylabel(r'$p_{\rm gas}$ [bar]',fontsize=16)

# added these three lines
ax2.set_zorder(1)
lns = p_rc + p_rp + p_nc + p_Ca
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='lower right',ncol=2)


plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Y_325_mono_gamma.pdf',dpi=144,bbox_inches='tight')

plt.show()



