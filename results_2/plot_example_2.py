import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

n_sp = 2

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


fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(q_v,pl,c=col[0])
plt.plot(q_0,pl,c=col[1])
plt.plot(q_1,pl,c=col[2])


plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

plt.xlim(1e-20,1)



plt.show()



