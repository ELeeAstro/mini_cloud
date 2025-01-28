import numpy as np

sb = 5.670374419e-8

T = 548.0
p = 4e4
g = 288.0
cp = 12045.0 
Rd = 3568.0

om = (2.0*np.pi)/(10.0*60.0*60.0)

R = 67202480.0

print(om)

tau_rad = p/g * (cp/(4.0*sb*T**3))

print('tau rad: ', tau_rad)

print('Tau: ', 2.0*om*tau_rad)


gam = 1.0
phi = gam * (Rd/cp)**2 * cp * T

print(phi)

RoT = phi/(2.0*om*R)**2

print(RoT) 

