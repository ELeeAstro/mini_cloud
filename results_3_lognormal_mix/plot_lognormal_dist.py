import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors


kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
V_seed = 4.0/3.0 * np.pi * r_seed**3

fname = 'tracers.txt'

ndust = 4
rho_d = [4.23, 3.986, 7.874, 3.21]
m_seed = V_seed * rho_d[0]

fig, ax = plt.subplots()

na = 1000

r_min = r_seed
r_max = 100.0 * 1e-4

m_min = 4.0/3.0 * np.pi * r_min**3 * rho_d[0]
m_max = 4.0/3.0 * np.pi * r_max**3 * rho_d[0]

m = np.logspace(np.log10(m_min),np.log10(m_max),na)

data = np.loadtxt(fname,skiprows=3)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:8]
q_v = data[:,8:12]
q_0 = data[:,12]
q_1 = data[:,13:17]
q_2 = data[:,17:21]
vf = data[:,21:30]
dTdt = data[:,30]

nlay = len(pl)

nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

N_c = np.zeros(nlay)
N_c[:] = q_0[:]*nd_atm[:]

rho = np.zeros(nlay)
rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

rho_c = np.zeros((nlay,ndust))
rho_c_t = np.zeros(nlay)
for j in range(ndust):
  rho_c[:,j] = q_1[:,j]*rho[:]
  rho_c_t[:] = rho_c_t[:] + rho_c[:,j]

m_c = np.zeros(nlay)
m_c[:] = np.maximum(rho_c_t[:]/N_c[:],m_seed)

rho_d_m = np.zeros(nlay)
for j in range(ndust):
  rho_d_m[:] = rho_d_m[:] + (rho_c[:,j])/rho_c_t[:] * rho_d[j]

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0),r_seed) * 1e4

Z_c_t = np.sum(q_2[:,:],axis=1)*rho[:]**2
lnsig_raw = np.sqrt(np.maximum(np.log(N_c[:] * Z_c_t[:] / rho_c_t[:]**2), 0.0))
sig_g = np.clip(np.exp(lnsig_raw), 1.01, 3.0)
lnsig2 = np.log(sig_g)**2
lnsig = np.sqrt(lnsig2)

m_med = m_c[:] * np.exp(-0.5 * lnsig2[:])

fx = np.zeros((na,nlay))
rr = np.zeros((na,nlay))
for n in range(na):
  fx[n,:] = N_c[:] / (np.sqrt(2.0*np.pi) * lnsig[:]) * \
    np.exp(-0.5 * (np.log(m[n]/m_med[:]))**2 / lnsig2[:])
  rr[n,:] = ((3.0*m[n])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0)

normalize = mcolors.Normalize(vmin=np.log10(pl[0]), vmax=np.log10(pl[-1]))
cmap = sns.color_palette("crest", as_cmap=True)

for i in range(0,nlay,5):
  if lnsig[i] >= 0.001:
    plt.plot(rr[:,i]*1e4,fx[:,i],c=cmap(normalize(np.log10(pl[i]))))

scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=cmap)
scalarmappaple.set_array(np.log10(pl))
cbar = plt.colorbar(scalarmappaple,ax=ax)

cticks = np.linspace(np.log10(pl[0]),np.log10(pl[-1]),10)
cbar.set_ticks(cticks)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel(r'$\log_{10}$ $p$ [bar]',fontsize=14)

plt.yscale('log')
plt.xscale('log')

plt.tick_params(axis='both',which='major',labelsize=14)

plt.xlabel(r'$r$ [$\mu$m]',fontsize=16)
plt.ylabel(r'$f(m)$ [cm$^{-3}$]',fontsize=16)

plt.ylim(1e-3,1e1)
plt.xlim(1e-3,1e2)

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.show()
