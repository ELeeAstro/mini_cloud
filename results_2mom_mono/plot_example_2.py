import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906660e-24
r_seed = 1e-7
rho_d = 1.99

fname = 'tracers.txt'

data = np.loadtxt(fname)

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

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(vf,pl,c=col[3],label=r'v_{\rm f}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(nd,pl,c=col[0],label=r'N_{\rm c}')

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

plt.plot(q_v,pl,c=col[0],label=r'q_{\rm v}')
plt.plot(q_0,pl,c=col[1],label=r'q_{0}')
plt.plot(q_1,pl,c=col[2],label=r'q_{1}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

#plt.xlim(1e-20,1)




plt.show()



