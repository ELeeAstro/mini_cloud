import numpy as np
import matplotlib.pylab as plt
import seaborn as sns


R = 8.31446261815324e7
bar = 1e6
atm = 1.01325e6
kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7

fname = 'tracers.txt'

ndust = 4
rho_d = [4.23, 3.21, 7.874, 3.986]
mol_w_sp = [79.8658, 140.693, 55.8450, 101.961]

Rd_v = [R/mol_w_sp[0],R/mol_w_sp[1],R/mol_w_sp[2],R/mol_w_sp[3]]

data = np.loadtxt(fname,skiprows=3)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:8]
q_v = data[:,8:12]
q_0 = data[:,12]
q_1 = data[:,13:17]
vf = data[:,17]

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

r_c = np.zeros(nlay)
r_c[:] = np.maximum(((3.0*m_c[:])/(4.0*np.pi*rho_d_m[:]))**(1.0/3.0),r_seed) * 1e4

q_s = np.zeros((nlay,ndust))
q_s[:,0] = (np.exp(-7.70443e4/Tl[:] +  4.03144e1 - 2.59140e-3*Tl[:] + \
  6.02422e-7*Tl[:]**2 - 6.86899e-11*Tl[:]**3))/(Rd_v[0] * Tl[:])
q_s[:,0] = q_s[:,0]/rho[:]

q_s[:,1] = (np.exp(-62279.0/Tl[:] + 20.944) * atm)/(Rd_v[1] * Tl[:])
q_s[:,1] = q_s[:,1]/rho[:]

q_s[:,2] = (np.exp(15.71 - 47664.0/Tl[:]) * bar)/(Rd_v[2] * Tl[:])
q_s[:,2] = q_s[:,2]/rho[:]

q_s[:,3] = (np.exp(-73503.0/Tl[:] + 22.005) * atm)/(Rd_v[3] * Tl[:])
q_s[:,3] = q_s[:,3]/rho[:]

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

plt.plot(q_v[:,0],pl,c=col[0],label=r'$q_{\rm v}$')
plt.plot(q_v[:,1],pl,c=col[1])
plt.plot(q_v[:,2],pl,c=col[2])
plt.plot(q_v[:,3],pl,c=col[3])

plt.plot(q_1[:,0],pl,c=col[0],label=r'$q_{1}$',ls='dashed')
plt.plot(q_1[:,1],pl,c=col[1],ls='dashed')
plt.plot(q_1[:,2],pl,c=col[2],ls='dashed')
plt.plot(q_1[:,3],pl,c=col[3],ls='dashed')

plt.plot(q_s[:,0],pl,c=col[0],label=r'$q_{\rm s}$',ls='dotted')
plt.plot(q_s[:,1],pl,c=col[1],ls='dotted')
plt.plot(q_s[:,2],pl,c=col[2],ls='dotted')
plt.plot(q_s[:,3],pl,c=col[3],ls='dotted')

plt.plot(q_0,pl,c=col[4],label=r'$q_{0}$',ls='dashdot')

plt.legend()

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

#plt.xlim(1e-20,1)


plt.show()



