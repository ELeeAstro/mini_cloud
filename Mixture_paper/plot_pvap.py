import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

def func(T, b, a):
    return b - a/T

bar = 1e6

nT = 10000
T = np.linspace(100,3000,nT)

a = [-2.69250E+04,  3.39574E+01, -2.04903E-03, -2.83957E-07,  1.82974E-10]

p_vap_sp = np.exp(a[0]/T +  a[1] + a[2]*T \
    + a[3]*T**2 + a[4]*T**3)

p_vap_sp_2 = 10.0**(7.611 - 11382.0/T) * bar

ln_p_vap_sp = np.log(p_vap_sp)

# popt, pcov = curve_fit(func, T, ln_p_vap_sp)
# print(popt)
# b = popt[0]
# a = popt[1]
# p_vap_fit = np.exp(b - a/T)

fig = plt.figure()

plt.plot(T,p_vap_sp)
plt.plot(T,p_vap_sp_2)
#plt.plot(T,p_vap_fit * bar,ls='dashed',c='black')

plt.yscale('log')

plt.show()