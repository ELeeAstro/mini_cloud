import numpy as np
import matplotlib.pylab as plt


Kn = np.logspace(-4,4,100)

A = 1.165
B = 0.483
C = 0.997

beta1 = 1.0 + Kn * (A + B*np.exp(-C/Kn))
beta2 = 1.0 + Kn * (A + B)


print((A + B))

plt.figure()


plt.plot(Kn, beta1)
plt.plot(Kn, beta2)

plt.xscale('log')

plt.figure()


plt.plot(Kn, beta1/beta2)

plt.xscale('log')

plt.show()
