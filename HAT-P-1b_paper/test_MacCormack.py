import numpy as np
import matplotlib.pylab as plt

def minmod(nlay,q,dh):
    sig = np.zeros(nlay)
    sig[1] = 0.0
    for i in range(1,nlay-2):
      de_minus = (q[i] - q[i-1]) / dh[i]
      de_plus = (q[i+1] - q[i]) / dh[i]
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

h = np.zeros(nlev)
p = np.logspace(3,-6,nlev)
p = p[::-1]
T = np.zeros(nlay)
vf = np.zeros(nlay)
q = np.zeros(nlay)

T[:] = 1000.0
vf[:] = -100.0
q[2] = 100.0


h[nlev-1] = 0.0
for k in range(nlay-1, -1, -1):
  h[k] = h[k+1] + (rdgas*T[k])/grav * np.log(p[k+1]/p[k])

print(h[:])

hmid = np.zeros(nlay)
hmid[:] = (h[0:nlev-1] +  h[1:nlev])/2.0
print(hmid[:])
dh = np.zeros(nlay)
for i in range(nlay-2,-1,-1):
  dh[i] =  hmid[i] - hmid[i+1]

dh[-1] = (hmid[-1] - h[-1])

print(dh[:])

dt = 9.0e99
for i in range(0,nlay):
  D  = abs(vf[i])
  if (D <= 0.0):
    continue
  dt = np.minimum(dt,cfl*(dh[i])/D)
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
        c[:] = abs(vf[:]) * dt / dh[:]
        sig = minmod(nlay,q,dh)
        qc[:] = q[:]
        qc[:nlay-1] = q[:nlay-1] - sig[:nlay-1]*c[:nlay-1]*(q[1:nlay] - q[:nlay-1])
        q[1:nlay] = 0.5 * (q[1:nlay] + qc[1:nlay] - c[1:nlay]*(qc[1:nlay] - qc[:nlay-1]))

        q[:] = np.maximum(q[:],1e-30)

        tnow = tnow + dt
        iit = iit + 1
        q[-1] = 1e-30

    print(q[:])