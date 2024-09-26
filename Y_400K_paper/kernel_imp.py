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

l_k = 1.0
nu = eta/rho
eps_d = nu**3/l_k**4

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
    E[i] = np.maximum(1e-6,1.0 - 0.42*Stk**(-0.75))

Re = 2.0*r*vf/nu
Stk = (vf * d_vf)/(grav * r)

Ev = np.zeros(nr)
Ea = np.zeros(nr)
for i in range(nr):
  if (Stk[i] > 1.214):
    Ev[i] = (1.0 + (0.75*np.log(2.0*Stk[i]))/(Stk[i] - 1.214))**-2
    Ea[i] = Stk[i]**2/(Stk[i] + 0.5)**2
  else:
    Ev[i] = 0.0
    Ea[i] = Stk[i]**2/(Stk[i] + 0.5)**2

#E = (60.0*Ev + Ea*Re)/(60.0 + Re)

# Coalesence flux (Zeroth moment) [cm3 s-1]
f_coal = 2.0*np.pi*r**2*d_vf*E

# Shear kernel
w_shear2 = 1.0/15.0 *  (2.0*r)**2 * eps_d/nu
f_shear = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_shear2)

# Gravity kernel
v_f = (2.0 * beta * grav * r**2 * rho_d)/(9.0 * eta) \
   * (1.0 + \
   ((0.45*grav*r**3*rho*rho_d)/(54.0*eta**2))**(0.4))**(-1.25)
tau = v_f / grav

#tau = (2.0*beta * (rho_d - rho)*r**2)/(9.0*eta)
dtau = 0.5 * tau
w_grav2 = np.pi/8.0 * (eps*dtau)**2 * (1.0 - rho/rho_d)**2 * grav**2
f_grav = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_grav2)


# Coupling
lam_d = 10.0
dudt2 = 1.16 * eps_d**(3.0/2.0)*nu**(-0.5)
w_coup2 = 2.0*(1.0 - rho/rho_d)**2*tau**2*dudt2*(2.0*r)**2/lam_d**2
f_coup = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_coup2)

# Turbulent inertial collision rate 2 [cm3 s-1]
cont = np.sqrt(2.0) * (2.0*r)**2
gam = 10.0
vf2 = (gam * np.sqrt(eps_d*nu))/0.183
tau_i = v_f / grav
b = (3.0*rho)/(2.0*rho_d + rho)
T_L = (0.4*vf2)/eps_d

th_i = tau_i/T_L

c1 = np.sqrt((1.0 + th_i + th_i)/((1.0 + th_i)*(1.0 + th_i)))
c2 = (1.0/(((1.0 + th_i)*(1.0 + th_i))) - 1.0/(((1.0 + gam*th_i)*(1.0 + gam*th_i))))
w_acc = 3.0*(1.0-b)**2*vf2*(gam/(gam-1.0)) * (((th_i + th_i)**2 - 4.0*th_i*th_i*c1)/(th_i + th_i)) * c2

f_acc = cont * np.sqrt(w_acc)

f_tot = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_shear2 + w_acc + w_coup2) + f_coag + f_coal

fig = plt.figure()

c = sns.color_palette('colorblind')

r = r * 1e4

plt.plot(r,f_tot,c='black',label=r'Total',ls='dotted',zorder=4)
plt.plot(r,f_coag,c=c[0],label=r'Brownian coagulation',ls='dashed')
plt.plot(r,f_coal,c=c[1],label=r'Gravitational coalescence',ls='dashed')
plt.plot(r,f_shear,c=c[2],label=r'Turbulent shear')
plt.plot(r,f_acc,c=c[3],label=r'Turbulent acceleration')
#plt.plot(r,f_grav,c=c[4],label=r'Turbulent gravity')
plt.plot(r,f_coup,c=c[4],label=r'Turbulent coupling')


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

quit()


fig = plt.figure()

c = sns.color_palette('colorblind')

plt.plot(r,Kn,c='black',label=r'Total',ls='dotted',zorder=4)


plt.xscale('log')
plt.yscale('log')


plt.xlim(1e-3,1e2)
#plt.ylim(1e-13,1e-4)

plt.ylabel(r'Kn', fontsize=16)
plt.xlabel(r'$r_{\rm c}$ [$\mu$m]', fontsize=16)
plt.tick_params(axis='both',which='major',labelsize=14)
plt.legend(title=r'p = 1 mbar, T = 750 K')

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)
#plt.savefig('K_rate_1.pdf',dpi=144,bbox_inches='tight')

plt.show()


# Turbulent inertial collision rate [cm3 s-1]
# f_ti = 2.0*((np.pi*eps_d**0.75)/(grav*nu**0.25)) * d_vf * r**2
# #f_ti = cont * np.sqrt(1.0/5.0 * 2.0*r**2 * eps_d/nu)


# #Turbulent shear collision rate [cm3 s-1]
# f_tsh = 4.0*np.sqrt((8.0*np.pi*eps_d)/(15.0*nu)) * r**3

# #Turbulent shear collision rate 2 [cm3 s-1]
# vi2vf2 = gam/(gam - 1.0) * ((1.0 + b**2*th_i)/(1.0 + th_i)) - ((1.0 + b**2*gam*th_i)/(gam*(1.0 + gam*th_i)))

# top = (th_i + th_i + 2.0*th_i*th_i)*b*(th_i**2 + th_i**2 - 2.0*th_i*th_i) + b**2*(th_i**2*th_i + th_i*th_i**2 + 2.0*th_i*th_i)
# bot = (th_i + th_i)*(1.0 + th_i)*(1.0 + th_i)
# t1 = top/bot
# top = (th_i + th_i + 2.0*gam*th_i*th_i)*b*gam*(th_i**2 + th_i**2 - 2.0*th_i*th_i) + b**2*(gam**2*th_i**2*th_i + gam**2*th_i*th_i**2 + 2.0*gam*th_i*th_i)
# bot = gam*(th_i + th_i)*(1.0 + gam*th_i)*(1.0 + gam*th_i)
# t2 = top/bot
# v1v2vf2 =  gam/(gam - 1.0) * (t1 - t2)

# ws2 = 0.283*b*vf2*(vi2vf2*(th_i/beta) + vi2vf2*(th_i/beta) + 2.0*v1v2vf2*np.sqrt((th_i*th_i)/(beta*beta)))

# Acceleration
# w_acc2 = (1.0 - rho/rho_d)**2 * (eps*dtau)**2 * dudt2
# f_acc = np.sqrt(2.0)*(2.0*r)**2*np.sqrt(w_acc2)
