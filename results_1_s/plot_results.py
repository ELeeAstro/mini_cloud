import numpy as np
import matplotlib.pylab as plt


fname = 'tracers.txt'
data = np.loadtxt(fname)

#t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_c, v_f
time = data[:,1]
T = data[:,2]
q_v = data[:,-3]
q_c = data[:,-2]
v_f = data[:,-1]


fig, ax1 = plt.subplots()
ax1.plot(time, T, ls='dashed', c='red', label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
ax2 = ax1.twinx()
ax2.plot(time, q_v, c='blue', label='q_v')
ax2.plot(time, q_c, c='grey', label='q_c')
ax2.set_ylabel('q [kg kg$^{-1}$]')
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


fname = 'opac.txt'
data = np.loadtxt(fname)

plt.show()