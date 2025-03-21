import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.colors as mcolors
from scipy.special import gamma


kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99
V_seed = 4.0/3.0 * np.pi * r_seed**3
m_seed = V_seed * rho_d

dirs = ['../results_3_gamma/']
ndir = len(dirs)

fname = 'tracers.txt'

fig, ax =  plt.subplots()

col = sns.color_palette('colorblind')

lss = ['solid']

na = 1000
r = np.logspace(-3,2,na) * 1e-4
m = 4.0/3.0 * np.pi * r**3 * rho_d

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
  vf = data[:,12]

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

  sig2 = np.zeros(nlay)  
  sig2[:] = np.maximum((q_2[:]*rho[:]**2)/(q_0[:]*nd_atm[:]) - ((q_1[:]*rho[:])/(q_0[:]*nd_atm[:]))**2,m_seed**2)


  lam = np.zeros(nlay)
  lam[:] = sig2[:]/m_c[:]

  nu = np.zeros(nlay)
  nu[:] = m_c[:]**2/sig2[:]
  nu[:] = np.minimum(nu[:],100.0)

  fx = np.zeros((na,nlay))
  for n in range(na):
    fx[n,:] = nd[:]/(lam[:]**nu[:]*gamma(nu[:])) * m[n]**(nu[:] - 1.0) * np.exp(-m[n]/lam[:]) * m[n]


# setup the normalization and the colormap
normalize = mcolors.Normalize(vmin=np.log10(pl[0]), vmax=np.log10(pl[-1]))
cmap = sns.color_palette("crest", as_cmap=True)


for i in range(0,nlay,5):
  if nu[i] > 0.1: 
    print(pl[i],nu[i])
    plt.plot(r[:]*1e4,fx[:,i],c=cmap(normalize(np.log10(pl[i]))))

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=cmap)
scalarmappaple.set_array(np.log10(pl))
cbar = plt.colorbar(scalarmappaple,ax=ax)

cticks = np.linspace(np.log10(pl[0]),np.log10(pl[-1]),10)
cbar.set_ticks(cticks)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel('$\log_{10}$ $p$ [bar]',fontsize=14)

plt.yscale('log')
plt.xscale('log')

plt.tick_params(axis='both',which='major',labelsize=14)

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$m$ $\cdot$ $f(m)$ [cm$^{-3}$]',fontsize=16)

plt.ylim(1e-1,1e4)
#plt.ylim(1e-3,1e6)
plt.xlim(1e-3,1e2)


plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('Y_425_gamma_dist.pdf',dpi=144,bbox_inches='tight')

plt.show()



