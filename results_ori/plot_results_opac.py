import numpy as np
import matplotlib.pylab as plt

## Edit this code to plot the example output for testing
## TODO: add dimensionality to fortran output to avoid annoying editing of this script

wl = np.loadtxt('opac.txt',max_rows=1)
nwl = len(wl)

print(wl)

data1 = np.loadtxt('opac.txt',skiprows=1)
it = data1[:,0]
time = data1[:,1]
nt = len(time)

k_ext = np.zeros((nt,nwl))
k_ext[:,:] = data1[:,2:2+nwl]

a = np.zeros((nt,nwl))
a[:,:] = data1[:,2+nwl:2+nwl+nwl]

g = np.zeros((nt,nwl))
g[:,:] = data1[:,2+nwl+nwl:]


fig, ax1 = plt.subplots()
for i in range(nwl):
  plt.plot(time,k_ext[:,i],label='{:.2f}'.format(wl[i]))
plt.ylabel('$\kappa_{ext}$ [cm$^{2}$ g$^{-1}$]')
plt.yscale('log')
plt.legend()

fig, ax1 = plt.subplots()
for i in range(nwl):
  plt.plot(time,a[:,i],label='{:.2f}'.format(wl[i]))
plt.ylabel('a')
plt.yscale('log')
plt.legend()

fig, ax1 = plt.subplots()
for i in range(nwl):
  plt.plot(time,g[:,i],label='{:.2f}'.format(wl[i]))
plt.ylabel('g')
plt.yscale('log')
plt.legend()

plt.show()