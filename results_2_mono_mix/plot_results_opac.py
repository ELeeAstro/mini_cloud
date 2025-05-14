import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

## Edit this code to plot the example output for testing
## TODO: add dimensionality to fortran output to avoid annoying editing of this script

wl = np.loadtxt('opac_k.txt',max_rows=1)
nwl = len(wl)

print(wl)

data_k = np.loadtxt('opac_k.txt',skiprows=1)

pl_k = data_k[:,2]/1e6
k_ext = data_k[:,3:]

data_a = np.loadtxt('opac_a.txt',skiprows=1)

pl_a = data_a[:,2]/1e6
ssa = data_a[:,3:]


data_g = np.loadtxt('opac_g.txt',skiprows=1)

pl_g = data_g[:,2]/1e6
g = data_g[:,3:]

fig = plt.figure()

col = sns.color_palette("husl", nwl)

for i in range(nwl):
  plt.plot(k_ext[:,i],pl_k[:],c=col[i],label='{:.2f}'.format(wl[i]))

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()
plt.legend()

fig = plt.figure()

col = sns.color_palette("husl", nwl)

for i in range(nwl):
  plt.plot(ssa[:,i],pl_a[:],c=col[i],label='{:.2f}'.format(wl[i]))

plt.yscale('log')
#plt.xscale('log')

plt.gca().invert_yaxis()
plt.legend()

fig = plt.figure()

col = sns.color_palette("husl", nwl)

for i in range(nwl):
  plt.plot(g[:,i],pl_g[:],c=col[i],label='{:.2f}'.format(wl[i]))

plt.yscale('log')
#plt.xscale('log')

plt.gca().invert_yaxis()
plt.legend()

plt.show()