import numpy as np
import matplotlib.pylab as plt

Kn = np.logspace(-3,3,100)

print(Kn)

beta = 1.0 + Kn*(1.165 + 0.483*np.exp(-0.997/Kn))

figure = plt.figure()

plt.plot(Kn, beta)
plt.xscale('log')
plt.show()
