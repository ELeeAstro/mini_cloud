import numpy as np
import matplotlib.pylab as plt
import seaborn as sns


n_dust = 4


kb = 1.380649e-16
amu = 1.66053906660e-24

fname = 'tracers.txt'
data = np.loadtxt(fname)

#t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_0, q_1, v_f
time = data[:,1]
T = data[:,2]
P = data[:,3]
grav = data[:,4]
mu = data[:,5]
rho_d = data[:,6]
VMR_bg = data[:,7:7+3]
i1 = 7+3
q_v = data[:,i1:i1+n_dust]
i2 = i1+n_dust
q_0 = data[:,i2]
i3 = i2 + 1
q_1s = data[:,i3:i3+n_dust]
i4 = i3+n_dust
v_f = data[:,i4]

nd_atm = (P*10.0)/(kb*T)  
rho = (P*10.0*mu*amu)/(kb * T) 

m_c = np.zeros(len(time))
for i in range(len(time)):
  m_c[i] = (sum(q_1s[i,:])*rho[i])/(q_0[i]*nd_atm[i])

r_c = ((3.0*m_c)/(4.0*np.pi*rho_d*1000.0))**(1.0/3.0) * 1e4

col = sns.color_palette('colorblind')

fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='black', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
for i in range(n_dust):
  ax2.plot(time, q_v[:,i], label='q_v'+str(i),ls='dotted',c=col[i])
ax2.plot(time, q_0, c='black', label='q_0')
for i in range(n_dust):
  ax2.plot(time, q_1s[:,i], label='q_1'+str(i),c=col[i])
ax2.set_ylabel('q')
ax2.legend(loc=2)
plt.yscale('log')

fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='red', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
ax2.plot(time, v_f, c='blue', label='v_f')
ax2.set_ylabel('v_f [m s$^{-1}$]')
ax2.legend(loc=2)
plt.yscale('log')

fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='red', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
ax2.plot(time, m_c, c='blue', label='m_c')
ax2.plot(time, r_c, c='black', label='r_c')
ax2.set_ylabel('m,r')
ax2.legend(loc=2)
plt.yscale('log')

fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='red', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
ax2.plot(time, rho_d, c='blue', label='rho_d')
ax2.set_ylabel(r'$\rho_{\rm d}$ [g cm$^{-3}$]')
ax2.legend(loc=2)
#plt.yscale('log')

plt.show()
