import numpy as np
import matplotlib.pylab as plt



kb = 1.380649e-23
amu = 1.66053906892e-27
R_gas = 8.31446261815324

fname = 'tracers.txt'
data = np.loadtxt(fname)

#t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_0, q_1, q_2, v_f
time = data[:,1]
T = data[:,2]
P = data[:,3]
grav = data[:,4]
mu = data[:,5]
q_v = data[:,-5]
q_0 = data[:,-4]
q_1 = data[:,-3]
q_2 = data[:,-2]
v_f = data[:,-1]

rho_d = 4.23
nd_atm = P/(kb*T)  
rho = (P*mu*amu)/(kb * T) 

m_c = (q_1*rho)/(q_0*nd_atm)
r_c = ((3.0*m_c)/(4.0*np.pi*rho_d*1000.0))**(1.0/3.0) * 1e6

fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='red', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
ax2.plot(time, q_v, c='blue', label='q_v')
ax2.plot(time, q_0, c='black', label='q_0')
ax2.plot(time, q_1, c='grey', label='q_1')
ax2.plot(time, q_2, c='orange', label='q_2')
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


fname = 'opac.txt'
data = np.loadtxt(fname)

plt.show()