import numpy as np

kb = 1.380649e-16
amu = 1.66053906660e-23

d_H2 = 2.827e-8 
LJ_H2 = 59.7 * kb
molg_H2 = 2.01588
d_He = 2.511e-8
LJ_He = 10.22 * kb
molg_He = 4.002602

mixr_g = [0.85, 0.15]
d_g = [d_H2, d_He]
LJ_g = [LJ_H2, LJ_He]
molg_g = [molg_H2, molg_He]

ngmix = 2

T = 1000.0
mu = 2.33
P = 100.0 * 1e6
rho = (P * mu * amu)/(kb * T)
grav = 15.0 * 100.0

a_av = 1.0 * 1e-4

rho_mix = 3.0

nu_g = np.zeros(ngmix)

#Find the dynamical viscosity of each gas
for g in range(ngmix):
    # Dynamical viscosity formula - Rosner (2000/2012) using Ackerman & Marley (2001) constants
    nu_g[g] = (5.0/16.0) * (np.sqrt(np.pi*(molg_g[g]*amu)*kb*T)/(np.pi*d_g[g]**2)) \
       * (((kb*T)/LJ_g[g])**(0.16)/1.22)

# !! Now we loop twice over to find the mixture value using the Wilke (1950) mixing rule
nu_mix = 0.0
for i in range(ngmix):
    nu_sum = 0.0
    for j in range(ngmix):
      phi_ij_top = (1.0 + np.sqrt(nu_g[i]/nu_g[j]) * (molg_g[j]/molg_g[i])**(0.25))**2
      phi_ij_bot = (4.0/np.sqrt(2.0)) * np.sqrt(1.0 + (molg_g[i]/molg_g[j]))
      phi_ij = phi_ij_top  / phi_ij_bot
      nu_sum = nu_sum + mixr_g[j] * phi_ij
    nu_mix = nu_mix + (mixr_g[i] * nu_g[i]) / nu_sum

#nu_mix = nu_g[0]

#For consistency, we now use the dynamical viscosity to find the mean free path of the layer
l_scale = (nu_mix/P) * np.sqrt((np.pi * kb * T) / (2.0 * mu * amu))

#Knudsen number and Cunningham slip factor at mean grain size
Kn = l_scale / a_av
beta = 1.0 + Kn * (1.256 + 0.4 * np.exp(-1.1/Kn))

# Final v_f is negative (downward)
vf = -(2.0 * beta * a_av**2 * grav * (rho_mix - rho)) / (9.0 * nu_mix)
print(vf, vf*-rho*grav, nu_mix, l_scale, Kn, beta)