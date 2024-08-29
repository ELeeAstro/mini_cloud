import numpy as np
import matplotlib.pylab as plt
from scipy.special import gamma
import seaborn as sns

#Results

rho_d = 1.0

Nc = 1319124335.8608873/1e6
rhoc = 4.2398316889244740E-006/1e6*1000.0
Zc = 2.9241692112278429E-014/1e6*(1000.0)**2

m_c = rhoc/Nc

r_c = ((3.0*m_c)/(4.0*np.pi*rho_d))**(1.0/3.0) * 1e4

Ev_m = rhoc/Nc
Var_m = Zc/Nc - Ev_m**2

Ev_m = 1e-11
Var_m = 1e-11

m_c = Ev_m

r_c = ((3.0*m_c)/(4.0*np.pi*rho_d))**(1.0/3.0) * 1e4

print(m_c, r_c)
print(Ev_m, Var_m)

na = 1000
a = np.logspace(-5,-3,na)
m = 4.0/3.0 * np.pi * a**3 * rho_d

Nc = 1.0

# log normal
S = np.exp(np.sqrt(np.log(1.0 + Var_m**2/Ev_m**2)))
r_m = Ev_m/(np.exp(np.log(S)**2))
print(r_m, S)
fln_m = Nc/(m*S*np.sqrt(2.0*np.pi)) * np.exp(-(np.log(m) - np.log(r_m))**2/(2.0*S**2))

fln = (fln_m*m) / (a/3.0) 

# # Inv Gamma
alpha_m = Ev_m**2/Var_m + 2.0
print(alpha_m)
beta_m = Ev_m*(alpha_m - 1.0)
figam_m = Nc*beta_m**(alpha_m)/gamma(alpha_m) * (1.0/m)**(alpha_m+1.0) * np.exp(-beta_m/m)

figam = (figam_m*m) / (a/3.0)

# # Gamma
alpha_m = Ev_m**2/Var_m
beta_m = Ev_m/Var_m
print(alpha_m, beta_m)
fgam_m = Nc*beta_m**(alpha_m)/gamma(alpha_m) * m**(alpha_m-1.0) * np.exp(-beta_m*m)

fgam = (fgam_m*m) / (a/3.0)

# #Exponential
lam_m = 1.0/Ev_m
fexp_m = Nc*lam_m*np.exp(-lam_m*m)

fexp = (fexp_m*m) / (a/3.0)

# # Rayleigh distribution
sig_m = Ev_m/np.sqrt(np.pi/2.0)
fRay_m = Nc*m/sig_m**2*np.exp(-m**2/(2.0*sig_m**2))

fRay = (fRay_m*m) / (a/3.0)


fln_int = np.trapezoid(fln_m,m)
fln_int2 = np.trapezoid(fln,a)

figam_int = np.trapezoid(figam_m,m)
figam_int2 = np.trapezoid(figam,a)

fgam_int = np.trapezoid(fgam_m,m)
fgam_int2 = np.trapezoid(fgam,a)

fexp_int = np.trapezoid(fexp_m,m)
fexp_int2 = np.trapezoid(fexp,a)

fRay_int = np.trapezoid(fRay_m,m)
fRay_int2 = np.trapezoid(fRay,a)

print(fln_int,  fln_int2)
print(figam_int,  figam_int2)
print(fgam_int,  fgam_int2)
print(fexp_int,  fexp_int2)
print(fRay_int,  fRay_int2)


fig = plt.figure()

c = sns.color_palette('colorblind') 

plt.plot(a*1e4,fln*a,label='log normal',c=c[0])
plt.plot(a*1e4,figam*a,label='Inv. Gamma',c=c[1])
plt.plot(a*1e4,fgam*a,label='Gamma',c=c[2])
plt.plot(a*1e4,fexp*a,label='Exponential',c=c[3])
plt.plot(a*1e4,fRay*a,label='Rayleigh',c=c[4])
plt.ylabel('f(a) [cm$^{-3}$ cm$^{-1}$]',fontsize=16)
plt.xlabel('a [$\mu$m]',fontsize=16)
plt.xscale('log')
#plt.xlim(1e-2,1e1)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('dist.pdf',dpi=300,bbox_inches='tight')

fig = plt.figure()

c = sns.color_palette('colorblind') 

plt.plot(m,fln_m*m,label='log normal',c=c[0])
plt.plot(m,figam_m*m,label='Inv. Gamma',c=c[1])
plt.plot(m,fgam_m*m,label='Gamma',c=c[2])
plt.plot(m,fexp_m*m,label='Exponential',c=c[3])
plt.plot(m,fRay_m*m,label='Rayleigh',c=c[4])
plt.ylabel('f(m) [cm$^{-3}$ g$^{-1}$]',fontsize=16)
plt.xlabel('m [g]',fontsize=16)
plt.xscale('log')
#plt.xlim(1e-2,1e1)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend()
plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('dist.pdf',dpi=300,bbox_inches='tight')

plt.show()
