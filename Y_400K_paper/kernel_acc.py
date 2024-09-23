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

r = np.logspace(-7,-2,nr)

rho_d = 2.0
m = 4.0/3.0 * np.pi * r**3 * rho_d

eta = (5.0/16.0) * (np.sqrt(np.pi*(molg_H2*amu)*kb*T)/(np.pi*d_H2**2)) \
  * ((((kb*T)/LJ_H2)**(0.16))/1.22)

mfp = (2.0*eta/rho) * np.sqrt((np.pi * mu)/(8.0*R_gas*T))

#Begin r dependent calculations
Kn = mfp/r

beta = 1.0 + Kn*(1.257 + 0.4 * np.exp(-1.1/Kn))

# Coagulation
D_r = (kb*T*beta)/(6.0*np.pi*eta*r)
V_r = np.sqrt((8.0*kb*T)/(np.pi*m))
lam_r = (8.0*D_r)/(np.pi*V_r)
del_r = ((2.0*r + lam_r)**3 - (4.0*r**2 + lam_r**2)**(1.5))/(6.0*r*lam_r) \
  - 2.0*r
phi = 2.0*r/(2.0*r+ np.sqrt(2.0)*del_r) + (4.0*D_r)/(r*np.sqrt(2.0)*V_r) 

# Coagulation flux (Zeroth moment) [cm3 s-1]
f_coag = (4.0 * kb * T * beta)/(3.0 * eta * phi)

# Coalesence 
vf = (2.0 * beta * grav * r**2 * rho_d)/(9.0 * eta) \
   * (1.0 \
   + ((0.45*grav*r**3*rho*rho_d)/(54.0*eta**2))**(0.4))**(-1.25)

eps = 0.5
d_vf = eps * vf

E = np.zeros(nr)
for i in range(nr):
  if (Kn[i] >= 1.0):
    E[i] = 1.0
  else:
    Stk = (vf[i] * d_vf[i])/(grav * r[i])
    E[i] = np.maximum(0.0,1.0 - 0.42*Stk**(-0.75))

# Coalesence flux (Zeroth moment) [cm3 s-1]
f_coal = 2.0*np.pi*r**2*d_vf*E


l_k = 1.0
nu = eta/rho
eps_d = nu**3/l_k**4

# Acceleration
tau = (2.0*beta * (rho_d - rho)*r**2)/(9.0*eta)
dtau = 0.5 * tau
dudt2 = 1.16 * eps_d**(3.0/2.0)*nu**(-0.5)
w_acc2 = (1.0 - rho/rho_d)**2 * (eps*dtau)**2 * dudt2
f_acc = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_acc2)

# Turbulent inertial collision rate 2 [cm3 s-1]
cont = np.sqrt(8.0*np.pi/3.0) * (2.0*r)**2
gam = 10.0
vf2 = (gam * np.sqrt(eps_d*nu))/0.183
tau_i = (2.0*beta * (rho_d - rho)*r**2)/(9.0*eta)
b = (3.0*rho)/(2.0*rho_d + rho)
T_L = (0.4*vf2)/eps_d

th_i = tau_i/T_L

c1 = np.sqrt((1.0 + th_i + th_i)/((1.0 + th_i)*(1.0 + th_i)))
c2 = (1.0/(((1.0 + th_i)*(1.0 + th_i))) - 1.0/(((1.0 + gam*th_i)*(1.0 + gam*th_i))))
wa2 = 3.0*(1.0-b)**2*vf2*(gam/(gam-1.0)) * (((th_i + th_i)**2 - 4.0*th_i*th_i*c1)/(th_i + th_i)) * c2

f_acc_2 = cont * np.sqrt(wa2)


cont = np.sqrt(2.0/np.pi) * 2.0*np.pi*(2.0*r)**2
C = 1.0 + 0.6*np.exp(-1.0)
wa2_2 = 1.0/3.0 * C * vf2 * (wa2/vf2)
f_acc_3 = cont * np.sqrt(wa2_2)


fig = plt.figure()

c = sns.color_palette('colorblind')

r = r * 1e4

plt.plot(r,f_coag,c=c[0],label=r'Brownian coagulation')
plt.plot(r,f_coal,c=c[1],label=r'Gravitational coalescence')
plt.plot(r,f_acc,c=c[4],label=r'Turbulent acceleration 1')
plt.plot(r,f_acc_2,c=c[4],label=r'Turbulent acceleration 2',ls='dashed')
plt.plot(r,f_acc_3,c=c[4],label=r'Turbulent acceleration 3',ls='dotted')



plt.xscale('log')
plt.yscale('log')


plt.xlim(1e-3,1e2)
plt.ylim(1e-13,1e-4)

plt.ylabel(r'$K$ [cm$^{3}$ s$^{-1}$]', fontsize=16)
plt.xlabel(r'$r_{\rm c}$ [$\mu$m]', fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend(title=r'p = 1 mbar, T = 750 K')

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
plt.savefig('K_rate_2.pdf',dpi=144,bbox_inches='tight')

plt.show()


# Turbulent inertial collision rate [cm3 s-1]
# f_ti = 2.0*((np.pi*eps_d**0.75)/(grav*nu**0.25)) * d_vf * r**2
# #f_ti = cont * np.sqrt(1.0/5.0 * 2.0*r**2 * eps_d/nu)


# #Turbulent shear collision rate [cm3 s-1]
# f_tsh = 4.0*np.sqrt((8.0*np.pi*eps_d)/(15.0*nu)) * r**3