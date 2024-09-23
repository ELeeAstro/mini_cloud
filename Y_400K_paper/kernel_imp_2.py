import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

kb = 1.380649e-16
amu = 1.66053906892e-24
R_gas = 8.31446261815324e7

T = 750.0

p = 1e-3 * 1e6

mu = 2.33

d_H2 = 2.827e-8
LJ_H2 = 59.7 * kb
molg_H2 = 2.01588

grav = 1e4

nd_atm = p/(kb*T)  

rho = (p*mu*amu)/(kb * T)

nr = 1000

rho_d = 2.0

eta = (5.0/16.0) * (np.sqrt(np.pi*(molg_H2*amu)*kb*T)/(np.pi*d_H2**2)) \
  * ((((kb*T)/LJ_H2)**(0.16))/1.22)

mfp = (2.0*eta/rho) * np.sqrt((np.pi * mu)/(8.0*R_gas*T))

l_k = 1.0
nu = eta/rho
eps_d = nu**3/l_k**4

#Begin r dependent calculations


r2 = np.logspace(-7,-2,nr)
r1 = np.logspace(-7,-2,nr)

Kn1 = mfp/r1
Kn2 = mfp/r2

beta1 = 1.0 + Kn1*(1.257 + 0.4 * np.exp(-1.1/Kn1))
beta2 = 1.0 + Kn2*(1.257 + 0.4 * np.exp(-1.1/Kn2))


cont = np.sqrt((8.0*np.pi/3.0)) * (r1 + r2)**2


# Turbulent inertial collision rate 2 [cm3 s-1]
gam = 10.0
vf2 = (gam * np.sqrt(eps_d*nu))/0.183
b = (3.0*rho)/(2.0*rho_d + rho)

tau1 = (beta1 * (2.0 * rho_d + rho)*r1**2)/(9.0*eta)
tau2 = (beta2 * (2.0 * rho_d + rho)*r2**2)/(9.0*eta)

dvf2dt = 1.16 * eps_d**(3.0/2.0)/nu**0.5

# Stauffman kernel
f_s = cont * np.sqrt((3.0*(1.0 - rho/rho_d)**2 * (tau1 - tau2)**2 * dvf2dt + 1.0/5.0*(r1 + r2)**2*eps_d/nu))


T_L = (0.4*vf2)/eps_d

th1 = tau1/T_L
th2 = tau2/T_L

c1 = np.sqrt((1.0 + th1 + th2)/((1.0 + th1)*(1.0 + th2)))
c2 = (1.0/(((1.0 + th1)*(1.0 + th2))) - 1.0/(((1.0 + gam*th1)*(1.0 + gam*th2))))
wa2 = 3.0*(1.0-b)**2*vf2*(gam/(gam-1.0)) * (((th1 + th2)**2 - 4.0*th1*th2*c1)/(th1 + th2)) * c2


#Turbulent shear collision rate 2 [cm3 s-1]
vi2vf2_1 = gam/(gam - 1.0) * ((1.0 + b**2*th1)/(1.0 + th1)) - ((1.0 + b**2*gam*th1)/(gam*(1.0 + gam*th1)))
vi2vf2_2 = gam/(gam - 1.0) * ((1.0 + b**2*th2)/(1.0 + th2)) - ((1.0 + b**2*gam*th2)/(gam*(1.0 + gam*th2)))


top = (th1 + th2 + 2.0*th1*th2)*b*(th1**2 + th2**2 - 2.0*th1*th2) + b**2*(th1**2*th2 + th1*th2**2 + 2.0*th1*th2)
bot = (th1 + th2)*(1.0 + th1)*(1.0 + th2)
t1 = top/bot
top = (th1 + th2 + 2.0*gam*th1*th2)*b*gam*(th1**2 + th2**2 - 2.0*th1*th2) + b**2*(gam**2*th1**2*th2 + gam**2*th1*th2**2 + 2.0*gam*th1*th2)
bot = gam*(th1 + th2)*(1.0 + gam*th1)*(1.0 + gam*th2)
t2 = top/bot
v1v2vf2 =  gam/(gam - 1.0) * (t1 - t2)

ws2 = 0.283*b*vf2*(vi2vf2_1*(th1/beta1) + vi2vf2_2*(th2/beta2) + 2.0*v1v2vf2*np.sqrt((th1*th2)/(beta1*beta2)))

f_p = cont * np.sqrt(ws2 + wa2)

print(f_s,f_p)



fig = plt.figure()

col = sns.color_palette('colorblind')

r1 = r1 * 1e4

plt.plot(r1,f_s,c='black',label=r'fs',ls='dashed',zorder=4)
plt.plot(r1,f_p,label=r'fd',zorder=4,c=col[0])

plt.xscale('log')
plt.yscale('log')


plt.xlim(1e-3,1e2)
plt.ylim(1e-13,1e-4)

plt.ylabel(r'$K$ [cm$^{3}$ s$^{-1}$]', fontsize=16)
plt.xlabel(r'$r_{\rm c}$ [$\mu$m]', fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend(title=r'p = 1 mbar, T = 750 K')

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('K_rate_2.pdf',dpi=144,bbox_inches='tight')

plt.show()


