import numpy as np
import matplotlib.pylab as plt

def minmod(nlay,q,dp):
    sig = np.zeros(nlay)
    sig[1] = 0.0
    for i in range(1,nlay-2):
      de_minus = (q[i] - q[i-1]) / dp[i]
      de_plus = (q[i+1] - q[i]) / dp[i]
      if ((de_minus > 0.0) and (de_plus > 0.0)):
        sig[i] = np.minimum(de_minus, de_plus)
      elif ((de_minus < 0.0) and (de_plus < 0.0)):
        sig[i] = np.maximum(de_minus, de_plus)
      else:
        sig[i] = 0.0
    sig[-1] = 0.0
    return sig

nlay = 54
nlev = nlay + 1
rdgas = 3577.0
grav = 7.46
cfl = 0.9
mu = 2.33 /1000.0
kb = 1.380649e-23
amu =  1.66053906660e-27

p = np.logspace(3,-6,nlev) * 1e5
p = p[::-1]
pmid = (p[0:nlev-1] +  p[1:nlev])/2.0
T = np.zeros(nlay)
vf = np.zeros(nlay)
q = np.zeros(nlay)
rho = np.zeros(nlay)

T[:] = 1000.0
q[2] = 100.0
rho[:]  = (pmid * mu * amu) / (kb * T[:])
vf[:] = 100.0 * rho * grav 

dp = np.zeros(nlay)
for i in range(nlay-2,-1,-1):
  dp[i] =  pmid[i+1] - pmid[i]

dp[-1] = (p[-1] - pmid[-1])

print(dp)

dt = 9.0e99
for i in range(0,nlay):
  D  = abs(vf[i])
  if (D <= 0.0):
    continue
  dt = np.minimum(dt,cfl*(dp[i])/D)
dt = dt * (1.0 + 1.e-12)

nit = 1

c = np.zeros(nlay)
sig = np.zeros(nlay)
qc = np.zeros(nlay)

print(q[:])

for n in range(nit):

    tnow = 0.0
    tend = 30.0
    iit = 1

    while tnow < tend and iit < 1e6:
        if (tnow + dt > tend):
            dt = tend - tnow
        c[:] = abs(vf[:]) * dt / dp[:]
        sig = minmod(nlay,q,dp)
        qc[:] = q[:]
        qc[:nlay-1] = q[:nlay-1] - sig[:nlay-1]*c[:nlay-1]*(q[1:nlay] - q[:nlay-1])
        q[1:nlay] = 0.5 * (q[1:nlay] + qc[1:nlay] - c[1:nlay]*(qc[1:nlay] - qc[:nlay-1]))

        q[:] = np.maximum(q[:],1e-30)

        tnow = tnow + dt
        iit = iit + 1
        q[-1] = 1e-30

    print(q[:])