import numpy as np
import matplotlib.pylab as plt


# Scheme to try plot can calculate latent heat of various 
# Material from the saturation vapour pressure
# We need to convert the expression to ln p against A/T, where A/T is the slope
# Then latent heat is the slope * R / mw

R = 8.31446261815324e7

sp_list = []
L_list = []

# C
# p_vap_sp = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
mw = 12.01070
L = 8.65139e4 * R / mw
print('C',L)

# TiC

# SiC

# CaTiO3
# Kozasa et al. (1987)
#p_vap_sp = exp(-79568.2_dp/T + 42.0204_dp) * atm
mw = 135.9432
L = 79568.2 * R / mw
print('CaTiO3',L)

# TiO2
# GGChem 5 polynomial NIST fit
# p_vap_sp = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3*T &
  # + 6.02422e-7*T**2 - 6.86899e-11*T**3) * bar
mw = 79.8658
L = 7.70443e4 * R / mw # erg g-1
print('TiO2',L)

# Al2O3
# Kozasa et al. (1989)
# p_vap_sp = exp(-73503.0_dp/T + 22.01_dp) * atm
mw = 101.96128
L = 73503.0 * R / mw
print('Al2O3',L)

# Fe
# Ackerman & Marley et al. (2001)
# p_vap_sp = exp(9.86_dp - 37120.0_dp/T) * bar
mw = 55.8450
L = 37120.0 * R / mw
print('Fe',L)

# Mg2SiO4
# Kozasa et al. (1989)
# p_vap_sp = exp(-62279.0_dp/T + 20.944_dp) * atm
mw = 140.693
L = 62279.0 * R / mw
print('Mg2SiO4',L)

# MgSiO3
# Ackerman & Marley (2001)
# p_vap_sp = exp(-58663.0_dp/T + 25.37_dp) * bar
mw = 100.3887
L = 58663.0 * R / mw
print('MgSiO3',L)

# SiO2
# GGChem
# p_vap_sp = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
# & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) * bar
mw = 60.08430
L = 7.28086e4 * R / mw
print('SiO2',L)

# SiO
# Gail et al. (2013)
# p_vap_sp = exp(-49520.0_dp/T + 32.52_dp)
mw = 44.08490
L = 49520.0 * R / mw
print('SiO',L)

# Cr
# Morley et al. (2012)
# p_vap_sp = 10.0_dp**(7.490_dp - 20592.0_dp/T) * bar
mw = 51.99610
L = (20592.0*np.log(10.0)) * R /mw
print('Cr',L)

# MnS
# Morley et al. (2012)
# p_vap_sp = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
mw = 87.0030
L = (23810.0*np.log(10.0)) * R / mw
print('MnS',L)

# Na2S
# Morley et al. (2012)
# p_vap_sp =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
mw = 878.0445
L = (13889.0*np.log(10.0)) * R / mw
print('Na2S',L)

# ZnS
# Morley et al. (2012)
# p_vap_sp = 10.0_dp**(12.812_dp - 15873.0_dp/T) * bar
mw = 97.4450
L = (15873.0*np.log(10.0)) * R / mw
print('ZnS',L)

# KCl
# Morley et al. (2012)
# p_vap_sp = 10.0_dp**(7.611_dp - 11382.0_dp/T) * bar
mw = 74.5513
L = (11382.0*np.log(10.0)) * R / mw
print('KCl',L)

# NaCl
# Stull (1947)
# p_vap_sp = 10.0_dp**(5.07184_dp - 8388.497_dp / (T - 82.638_dp)) * bar
mw = 58.4428
L = (8388.497*np.log(10.0)) * R / mw
print('NaCl',L)

# NH4Cl
# p_vap_sp = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
mw = 53.4915
L = (4302.0*np.log(10.0)) * R / mw
print('NH4Cl',L)

# H2O

# NH3
# Ackerman & Marley (2001)
# p_vap_sp = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
# Dominant term approximation
mw = 17.03052
L = 2161.0 * R / mw
print('NH3',L)

# CH4
# Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1
# p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
mw = 16.0425
L = (443.028*np.log(10.0)) * R / mw
print('CH4',L)

# NH4SH
#  E.Lee's fit to Walker & Lumsden (1897)
# p_vap_sp = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
mw = 51.1114
L = (2409.4*np.log10(10.0)) * R / mw
print('NH4SH',L)

# H2S
# Stull (1947)
# p_vap_sp = 10.0_dp**(4.52887_dp - 958.587_dp/(T-0.539_dp)) * bar
mw = 34.0809
L = (958.587*np.log10(10.0)) * R / mw
print('H2S',L)

  