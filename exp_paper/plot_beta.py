import numpy as np
import matplotlib.pylab as plt

Kn = np.logspace(-3,3,100)

print(Kn)

beta_1 = 1.0 + Kn*(1.165 + 0.483*np.exp(-0.997/Kn))
beta_2 = 1.0 + Kn * 1.639

figure = plt.figure()

plt.plot(Kn, beta_1)
plt.plot(Kn, beta_2)

plt.xscale('log')

plt.show()
