import numpy as np
from scipy.special import gamma


cond1 = gamma(4.0/3.0)
cond2 = gamma(5.0/3.0)

print(cond1,cond2)

Kn = 100.0
beta = 1.0 + Kn*(1.165 + 0.483*np.exp(-0.997/Kn))

coag1 = (2.0*(1.0 + gamma(4.0/3.0)*gamma(2.0/3.0)))/4.0
coag1_b = (2.0*(1.0 + gamma(4.0/3.0)*gamma(2.0/3.0) + Kn*1.639 * (gamma(2.0/3.0) + gamma(4.0/3.0)*gamma(1.0/3.0))))/(4.0*beta)
coag2 = (np.sqrt(8.0)*(gamma(5.0/3.0)*gamma(1.0/2.0) + 2.0*gamma(4.0/3.0)*gamma(5.0/6.0) + gamma(7.0/6.0)))/8.0
coal = (gamma(5.0/3.0) + gamma(4.0/3.0)**2)/2.0

print(coag1,coag1_b,coag2,coal)

