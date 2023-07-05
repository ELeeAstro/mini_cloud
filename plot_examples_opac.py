import numpy as np
import matplotlib.pylab as plt

## Edit this code to plot the example output for testing
## TODO: add dimensionality to fortran output to avoid annoying editing of this script
example = 1

nwl = 11

data1 = np.loadtxt('opac.txt')
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
plt.plot(time,k_ext[:,0])
plt.plot(time,k_ext[:,-1])
plt.ylabel('k_ext')
plt.yscale('log')
plt.show()

fig, ax1 = plt.subplots()
plt.plot(time,a[:,0])
plt.plot(time,a[:,-1])
plt.ylabel('a')
plt.yscale('log')
plt.show()

fig, ax1 = plt.subplots()
plt.plot(time,g[:,0])
plt.plot(time,g[:,-1])
plt.ylabel('g')
plt.yscale('log')
plt.show()
