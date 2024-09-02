import numpy as np
import matplotlib.pylab as plt


# Scheme to try plot can calculate latent heat of various 
# Material from the saturation vapour pressure

sp = 'KCl'
R = 8.31446261815324


match sp:
  case 'KCl':
    # p_vap_sp = 10.0**(12.812 - 15873.0/T) * bar Neet to invert this to lnp vs 1/T
    mw = 74.5513
    L = (15873.0*np.log(10.0)) * R # J mol-1
    L = L / mw # J g mol -1
    #2923 kJ kg-1
    print(L, 15873.0*(np.log(10.0)))
  case _:
    print('Species not available!')