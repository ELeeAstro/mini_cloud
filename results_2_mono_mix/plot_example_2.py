import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

bar = 1e6
kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99

fname = 'tracers.txt'

ndust = 2
rho_d = [1.99, 3.9]
mol_w_sp = [74.5513, 97.4450]

data = np.loadtxt(fname,skiprows=3)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:8]
q_v = data[:,8:10]
q_0 = data[:,10]
q_1 = data[:,11:13]
vf = data[:,13]

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

m_c = np.zeros(nlay)
m_c[:] = rho_c_t[:]/N_c[:]

rho_d_m = np.zeros(nlay)
for j in range(ndust):
  rho_d_m[:] = rho_d_m[:] + (rho_c[:,j])/rho_c_t[:] * rho_d[j]
print(rho_d_m[:])

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0),r_seed) * 1e4

q_s = np.zeros((nlay,ndust))
q_s[:,0] = np.exp(-2.69250e4/Tl[:] + 3.39574e+1 - 2.04903e-3*Tl[:]  -2.83957e-7*Tl[:]**2 + 1.82974e-10*Tl[:]**3)
q_s[:,0] = np.maximum(q_s[:,0]/1e6/pl[:],1e-30)

q_s[:,1] = np.exp(-4.75507888e4/Tl[:] + 3.66993865e1 - 2.49490016e-3*Tl[:] + 7.29116854e-7*Tl[:]**2 - 1.12734453e-10*Tl[:]**3)
q_s[:,1] = np.maximum(q_s[:,1]/1e6/pl[:],1e-30)

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(vf,pl,c=col[3],label=r'v_{\rm f}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(N_c,pl,c=col[0],label=r'N_{\rm c}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(r_c,pl,c=col[1],label=r'r_{\rm c}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(q_v[:,0],pl,c=col[0],label=r'q_{\rm v}')
plt.plot(q_v[:,1],pl,c=col[1],label=r'q_{\rm v}')
plt.plot(q_0,pl,c=col[3],label=r'q_{0}')
plt.plot(q_1[:,0],pl,c=col[0],label=r'q_{1}',ls='dashed')
plt.plot(q_1[:,1],pl,c=col[1],label=r'q_{1}',ls='dashed')
plt.plot(q_s[:,0],pl,c=col[0],label=r'q_{s}',ls='dotted')
plt.plot(q_s[:,1],pl,c=col[1],label=r'q_{s}',ls='dotted')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

#plt.xlim(1e-20,1)




plt.show()



