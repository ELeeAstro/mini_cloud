import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma

nu = np.logspace(-2,2,1000)

print(nu)

Kn_av = nu**(1.0/3.0) * gamma(nu - 1.0/3.0)/gamma(nu)
Kn_av_m = nu**(1.0/3.0) * gamma(nu + 2.0/3.0)/gamma(nu + 1.0)
Kn_av_m2 = nu**(1.0/3.0) * gamma(nu + 5.0/3.0)/gamma(nu + 2.0)

fig = plt.figure()


plt.plot(nu, Kn_av)
plt.plot(nu, Kn_av_m)
plt.plot(nu, Kn_av_m2)

plt.yscale('log')
plt.xscale('log')

plt.show()



