import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

R = 8.31446261815324e7
kb = 1.380649e-16
amu = 1.66053906660e-24
fname = 'tracers.txt'

with open(fname, "r", encoding="utf-8") as f:
    nsp = int(f.readline().split()[0])
    species = f.readline().split()
    rho_d = float(f.readline().split()[0])
    mol_w_sp = float(f.readline().split()[0])
    r_med = float(f.readline().split()[0])
    sig, dist = f.readline().split()
    sig = float(sig)
    dist = int(dist)

if nsp != 1:
    raise ValueError(f"plot_example_2.py expects one species, got {nsp}: {species}")

Rd_v = R/mol_w_sp

data = np.loadtxt(fname, skiprows=6)

Tl = data[:,2]
pl = data[:,3]/1e5
grav = data[:,4]
mu = data[:,5]
VMR = data[:,6:7]
q_v = data[:,8]
q_c = data[:,9]
vf = data[:,10]

nlay = len(pl)

nd_atm = np.zeros(nlay)
nd_atm[:] = (pl[:]*1e6)/(kb*Tl[:])

rho = np.zeros(nlay)
rho[:] = (pl[:]*1e6*mu[:]*amu)/(kb * Tl[:])

q_s = np.zeros(nlay)
q_s[:] = (np.exp(-2.69250e4/Tl[:] + 3.39574e+1 - 2.04903e-3*Tl[:]  -2.83957e-7*Tl[:]**2 + 1.82974e-10*Tl[:]**3))/(Rd_v * Tl[:])
q_s[:] = q_s[:]/rho[:]

N_c = np.zeros(nlay)
N_c[:] = (3.0 * q_c[:] * rho[:])/(4.0*np.pi*rho_d*r_med**3)
if dist == 1:
  N_c[:] = N_c[:] * np.exp(-9.0/2.0 * np.log(sig)**2)

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(vf,pl,c=col[3],label=r'v_{\rm f}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(N_c,pl,c=col[0],label=r'N_{\rm c}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.axvline(r_med*1e4,c=col[1],label=r'r_{\rm med}')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

fig = plt.figure()

col = sns.color_palette('colorblind')

plt.plot(q_v,pl,c=col[0],label=r'q_{\rm v}')
plt.plot(q_c,pl,c=col[2],label=r'q_{\rm c}')
plt.plot(q_s,pl,c=col[3],label=r'q_{\rm s}',ls='dashed')

plt.yscale('log')
plt.xscale('log')

plt.gca().invert_yaxis()

plt.show()


