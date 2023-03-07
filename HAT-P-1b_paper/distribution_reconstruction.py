import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns

#Results

k0 = 215.19428555942420
k1 =  2.1588715875281238E-002
k2 =  2.3624659435413312E-006
k3 =    3.0427467964386399E-010

#k0 = 166.61488812049220
#k1 = 1.9451216474316985E-002
#k2 = 2.4047090067741973E-006
#k3 = 3.2214254633654980E-010

Ev = k1/k0
Var = k2/k0 - Ev**2

print(Ev, Var)

#Var = 1.9e-10

na = 1000
a = np.logspace(-6,-3,na)

# log normal
mu = np.log(Ev/np.sqrt(Var/Ev**2 + 1))
sig = np.sqrt(np.log(Var/Ev**2 + 1))
fln = k0/(a*sig*np.sqrt(2.0*np.pi)) * np.exp(-(np.log(a) - mu)**2/(2.0*sig**2))

# normal
#mu = Ev
#sig = np.sqrt(var)
#fnorm = k0/(a*sig*np.sqrt(2.0*np.pi)) * np.exp(-(np.log(a) - mu)**2/(2.0*sig**2))

# Gamma
alpha = Ev**2/Var
beta = Ev/Var
fgam = k0*beta**(alpha)/gamma(alpha) * a**(alpha-1) * np.exp(-beta*a)

# Inv Gamma
alpha = np.minimum(Ev**2/Var + 2.0,50)
beta = Ev*(alpha - 1.0)
figam = k0*beta**(alpha)/gamma(alpha) * (1.0/a)**(alpha+1) * np.exp(-beta/a)

# Rayleigh distribution
sig = Ev/np.sqrt(np.pi/2.0)
#sig = np.sqrt(Var/(2.0 - np.pi/2.0))
fRay = k0*a/sig**2*np.exp(-a**2/(2.0*sig**2))

#Potential Exponential
BB = (2.0*k1*k3 - 3.0*k2**2)/(k2**2 - k1*k3)
CC = (2.0 + BB)*k1/k2
AA = np.log(k1) + (2.0 + BB)*np.log(CC) - np.log(gamma(2.0 + BB))

fpot = a**BB*np.exp(AA - CC*a)

# Hansen
aeff = k3/k2
veff = (k1*k3)/k2**2 - 1.0
print(aeff,veff)
fHan = k0/gamma((1.0 - 2.0*veff)/veff) * (aeff*veff)**((2.0*veff - 1.0)/veff) * a**((1.0 - 3*veff)/veff) * np.exp(-a/(aeff*veff))

#Exponential
lam = 1.0/Ev
fexp = k0*lam*np.exp(-lam*a)

# Bell shaped
av = Ev
sig2 = (k2*k0 - k1**2)/k0**2

fbell = np.exp(-(a - av)**2/(2.0*sig))

# beta 
nu = (Ev*(1.0 - Ev))/Var - 1.0
alpha = Ev*nu
beta = (1.0 - Ev)*nu

Bbeta = (gamma(alpha)*gamma(beta))/gamma(alpha + beta)
print(gamma(alpha), gamma(beta), gamma(alpha + beta))
print(Bbeta,nu,alpha,beta)

fbeta = k0 * (a**(alpha - 1.0) * (1.0 - a)**(beta - 1.0))/Bbeta

# Initial guess for the parameters
#k_init = 1
#lambda_init = Ev / np.e
#xbar = Ev
#s2 = Var
# Define the function to optimize
#def weibull_moments(params, xbar, s2):
#    k, lambda_ = params
#    mean = lambda_ * gamma(1 + 1/k)
#    var = lambda_**2 * gamma(1 + 2/k) - mean**2
#    return (mean - xbar)**2 + (var - s2)**2
# Optimize the parameters
#from scipy.optimize import minimize
#params_est = minimize(weibull_moments, x0=[k_init, lambda_init], args=(xbar, s2))
# Estimated parameters
#k, lam = params_est.x
#fweb = k0 * k/lam * (a/lam)**(k-1.0)*np.exp(-(a/lam)**k)

fig = plt.figure()

c = sns.color_palette('colorblind') 

plt.plot(a*1e4,fln,label='log normal',c=c[0])
plt.plot(a*1e4,figam,label='Inv. Gamma',c=c[1])
plt.plot(a*1e4,fgam,label='Gamma',c=c[2])
plt.plot(a*1e4,fexp,label='Exponential',c=c[3])
plt.plot(a*1e4,fRay,label='Rayleigh',c=c[4])
plt.plot(a*1e4,fpot,label='Pot. Exp.',c=c[5])
plt.plot(a*1e4,fHan,label='Hansen',c=c[6])
#plt.plot(a*1e4,fbell,label='Bell',c=c[7])
plt.plot(a*1e4,fbeta,label='Beta',c=c[7])

#plt.plot(a*1e4,fweb,label='Web',c=c[7])


plt.ylabel('f(a) [cm$^{-3}$ cm$^{-1}$]',fontsize=16)
plt.xlabel('a [$\mu$m]',fontsize=16)

plt.xscale('log')

plt.xlim(1e-2,1e1)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.legend()

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('dist.pdf',dpi=300,bbox_inches='tight')

plt.show()
