import numpy as np
import matplotlib.pylab as plt

## Edit this code to plot the example output for testing
## TODO: add dimensionality to fortran output to avoid annoying editing of this script
example = 1

nmom = 4
ndust = 3
neps = 3

data1 = np.loadtxt('tracers.txt')
time = data1[:,1]
nt = len(time)
Tg = data1[:,2]
Pg = data1[:,3]
ntracer = nmom + ndust + neps
tracer = np.zeros((nt,ntracer))
tracer[:,:] = data1[:,4:]

amean = tracer[:,1]/tracer[:,0]  * 1e4
amean[np.where(tracer[:,0] < 1e-20)] = 0.0

Vsp = np.zeros((nt,ndust))
for i in range(ndust):
  Vsp[:,i] = tracer[:,nmom+i]/tracer[:,3]
Vsp[np.where(tracer[:,3] < 1e-20)] = 0.0

spname = ['TiO2','Fe','MgSiO3']
epsname = ['TiO2','Fe','MgSiO3']
c = []


if example == 1 :
  yval = Tg
  print('Example 1, pressure = ', Pg * 1e-5, 'bar')
elif example == 2:
  yval = Pg * 1e-5
  print('Example 2, temperature = ',Tg , 'K')


fig, ax1 = plt.subplots()
ax1.plot(time,yval, ls='dashed',c='red',label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
if example == 2:
    plt.yscale('log')
ax2 = ax1.twinx()
ax2.plot(time,tracer[:,0],label='nd')
ax2.set_ylabel('nd [cm-3]')
plt.yscale('log')
ax2.legend(loc=2)

fig, ax1 = plt.subplots()
ax1.plot(time,yval, ls='dashed',c='red',label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
if example == 2:
    ax1.set_ylabel('P [bar]')
    plt.yscale('log')
ax2 = ax1.twinx()
ax2.plot(time,amean,label='<a>')
ax2.set_ylabel('<a> [um]')
ax2.legend(loc=2)
plt.yscale('log')

fig, ax1 = plt.subplots()
ax1.plot(time,yval, ls='dashed',c='red',label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
if example == 2:
    ax1.set_ylabel('P [bar]')
    plt.yscale('log')
ax2 = ax1.twinx()
for i in range(ndust):
  ax2.plot(time,Vsp[:,i],label=spname[i])
ax2.set_ylabel('Vs')
ax2.legend(loc=2)
plt.yscale('log')

fig, ax1 = plt.subplots()
ax1.plot(time,yval, ls='dashed',c='red',label='T')
ax1.set_ylabel('T [K]')
ax1.set_xlabel('Time [s]')
ax1.legend(loc=1)
if example == 2:
    ax1.set_ylabel('P [bar]')
    plt.yscale('log')
ax2 = ax1.twinx()
for i in range(neps):
  ax2.plot(time,tracer[:,nmom+ndust+i],label=epsname[i])
for i in range(neps):
  ax2.axhline(tracer[0,nmom+ndust+i],ls='dotted')
ax2.set_ylabel('eps')
ax2.legend(loc=2)
plt.yscale('log')

plt.show()
