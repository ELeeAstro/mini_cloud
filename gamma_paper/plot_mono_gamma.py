import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

R = 8.31446261815324e7
kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99

mol_w_sp = 74.5513
Rd_v = R/mol_w_sp

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
#  q_2 = data[:,11]
#  vf = data[:,12]

  nlay = len(pl)

  nd_atm = np.zeros(nlay)
  nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

  nd = np.zeros(nlay)
  nd[:] = q_0[:]*nd_atm[:]

  rho = np.zeros(nlay)
  rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

  m_c = np.zeros(nlay)
  m_c[:] = (q_1[:]*rho[:])/(q_0[:]*nd_atm[:])

  r_c = np.zeros(nlay)
  r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d))**(1.0/3.0),r_seed) * 1e4

  q_s = np.zeros(nlay)
  q_s[:] = (np.exp(-2.69250e4/Tl[:] + 3.39574e+1 - 2.04903e-3*Tl[:]  -2.83957e-7*Tl[:]**2 + 1.82974e-10*Tl[:]**3))/(Rd_v * Tl[:])
  q_s[:] = q_s[:]/rho[:]

  p_T = ax1.plot(Tl,pl,c=col[0],label=r'T',ls=lss[i])
  p_qc = ax2.plot(q_1,pl,c=col[1],label=r'$q_{\rm 1}$',ls=lss[i])
  p_qv = ax2.plot(q_v,pl,c=col[2],label=r'$q_{\rm v}$',ls=lss[i])
  p_qs = ax2.plot(q_s,pl,c=col[4],label=r'$q_{\rm s}$',ls='dashdot')
  #p_0 = ax2.plot(q_0,pl,c=col[5],label=r'$q_{\rm 0}$',ls=lss[i])
  #p_2 = ax2.plot(q_2,pl,c=col[6],label=r'$q_{\rm 2}$',ls=lss[i])

ax1.set_yscale('log')
#ax1.set_xscale('log')
ax2.set_xscale('log')

plt.gca().invert_yaxis()
yticks = [100,10,1,0.1,0.01,1e-3,]
yticks_lab = ['100','10','1','0.1','0.01','10$^{-3}$']
ax1.set_yticks(yticks,yticks_lab)

ax1.set_xlim(0,1500)
ax2.set_xlim(1e-8,1e-5)
#ax2.set_xlim(1e-7,1e-3)

plt.ylim(300,3e-3)

ax1.tick_params(axis='both',which='major',labelsize=14)
ax2.tick_params(axis='both',which='major',labelsize=14)

ax1.set_xlabel(r'$T_{\rm gas}$ [K]',fontsize=16)
ax2.set_xlabel(r'$q$ [g g$^{-1}$]',fontsize=16)
ax1.set_ylabel(r'$p_{\rm gas}$ [bar]',fontsize=16)

# added these three lines
ax2.set_zorder(1)
lns = p_T + p_qc + p_qv + p_qs
labs = [l.get_label() for l in lns]
ax2.legend(lns, labs,fontsize=10,loc='lower left')


plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Y_325_mono_gamma_Tq.pdf',dpi=144,bbox_inches='tight')

plt.show()



