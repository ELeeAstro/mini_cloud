import numpy as np
import matplotlib.pylab as plt
import seaborn as sns


R = 8.31446261815324e7
kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
V_seed = 4.0/3.0 * np.pi * r_seed**3

fname = 'tracers.txt'

ndust = int(np.loadtxt(fname, max_rows=1))
rho_d = np.loadtxt(fname, skiprows=1, max_rows=1)
mol_w_sp = np.loadtxt(fname, skiprows=2, max_rows=1)
data = np.loadtxt(fname, skiprows=3)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:8]
q_v = data[:,8:8+ndust]
q_0 = data[:,8+ndust]
q_1 = data[:,9+ndust:9+2*ndust]
q_2 = data[:,9+2*ndust:9+3*ndust]
vf = data[:,9+3*ndust:10+5*ndust]
dTdt = data[:,10+5*ndust]

nlay = len(pl)

nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

rho = np.zeros(nlay)
rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

N_c = np.zeros(nlay)
N_c[:] = q_0[:]*nd_atm[:]

rho_c = np.zeros((nlay,ndust))
rho_c_t = np.zeros(nlay)
for j in range(ndust):
  rho_c[:,j] = q_1[:,j]*rho[:]
  rho_c_t[:] = rho_c_t[:] + rho_c[:,j]

m_seed = V_seed * rho_d[0]

m_c = np.zeros(nlay)
m_c[:] = np.maximum(rho_c_t[:]/N_c[:],m_seed)

rho_d_m = np.zeros(nlay)
for j in range(ndust):
  rho_d_m[:] = rho_d_m[:] + (rho_c[:,j])/rho_c_t[:] * rho_d[j]

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0),r_seed) * 1e4

V_mix = np.zeros((nlay,ndust))
for j in range(nlay):
  V_mix[j,:] = rho_c[j,:]/rho_d[:]/(sum(rho_c[j,:]/rho_d[:]))

sig2 = np.zeros(nlay)
sig2[:] = np.maximum(np.sum(q_2[:,:],axis=1)*rho[:]**2/N_c[:] - m_c[:]**2,m_seed**2)

nu = np.zeros(nlay)
nu[:] = m_c[:]**2/sig2[:]

fig = plt.figure()

col = sns.color_palette('colorblind')

for j in range(ndust):
  plt.plot(q_v[:,j],pl,c=col[j])
  plt.plot(q_1[:,j],pl,c=col[j],ls='dashed')
  plt.plot(q_2[:,j],pl,c=col[j],ls='dotted')
plt.plot(q_0,pl,c=col[ndust],ls='dashdot')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(N_c,pl,c=col[0],label=r'N_{\rm c}')
plt.plot(r_c,pl,c=col[1],label=r'r_{\rm c}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

for j in range(ndust):
  plt.plot(V_mix[:,j],pl,c=col[j])

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(nu,pl,c=col[0],label=r'$\nu$')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(dTdt,pl,c=col[0],label=r'$dTdt$')

plt.yscale('log')
plt.gca().invert_yaxis()

plt.show()
