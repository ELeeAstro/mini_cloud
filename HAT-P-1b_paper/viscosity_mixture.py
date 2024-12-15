import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

# Do H2 and He mixture
nsp = 2
kb = 1.380649e-16
amu = 1.66053906892e-24
sp = ['H2','He']
g_VMR = [0.85,0.15]
g_mu = [2.01588,4.002602]
g_LJ = [59.7 * kb, 10.22* kb]
g_d = [2.827e-8,2.511e-8]


nT = 1000
T = np.linspace(1,6000,nT)

eta_g = np.zeros((nsp,nT)) 
for i in range(nsp):
  eta_g[i,:] = (5.0/16.0) * (np.sqrt(np.pi*(g_mu[i]*amu)*kb*T[:])/(np.pi*g_d[i]**2)) \
    * ((((kb*T[:])/g_LJ[i])**(0.16))/1.22)

# Square root mixing law
eta_mix_1 = np.zeros(nT)
for n in range(nT):
  top = 0.0
  bot = 0.0
  for i in range(nsp):
     top = top + eta_g[i,n] * g_VMR[i] * np.sqrt(g_mu[i])
     bot = bot + g_VMR[i] * np.sqrt(g_mu[i])
  eta_mix_1[n] = top/bot

# Wilke mixing rule
eta_mix_2 = np.zeros(nT)
for n in range(nT):
  for i in range(nsp):
    eta_sum = 0.0
    for j in range(nsp):
      phi_ij_top = (1.0 + np.sqrt(eta_g[i,n]/eta_g[j,n]) * (g_mu[j]/g_mu[i])**(0.25))**2
      phi_ij_bot = np.sqrt(8.0*(1.0 + (g_mu[i]/g_mu[j])))
      phi_ij = phi_ij_top  / phi_ij_bot
      eta_sum = eta_sum + g_VMR[j] * phi_ij
    eta_mix_2[n] = eta_mix_2[n] + (g_VMR[i] * eta_g[i,n]) / eta_sum

# Davidson mixing rule
y = np.zeros(nsp)
bot = 0.0
for i in range(nsp):
  bot = bot + g_VMR[i] * np.sqrt(g_mu[i])

y[:] = (g_VMR[:] * np.sqrt(g_mu[:]))/bot

eta_mix_3 = np.zeros(nT)
for n in range(nT):
  for i in range(nsp):
    for j in range(nsp):
      Eij = ((2.0*np.sqrt(g_mu[i]*g_mu[j]))/(g_mu[i] + g_mu[j]))**0.375
      part = (y[i]*y[j])/(np.sqrt(eta_g[i,n]*eta_g[j,n])) * Eij
      eta_mix_3[n] = eta_mix_3[n] + part

eta_mix_3[:] = 1.0/eta_mix_3[:]

fig = plt.figure()

col = sns.color_palette('colorblind')


for i in range(nsp):
   plt.plot(T[:],eta_g[i,:],label=sp[i],c=col[i],ls='solid')

plt.plot(T[:],eta_mix_1[:],label='Sqrt mixing',c='black',ls='dashed')
plt.plot(T[:],eta_mix_2[:],label='Wilke',c='black',ls='dashdot')
plt.plot(T[:],eta_mix_3[:],label='Davidson',c='black',ls='dotted')



plt.legend()

plt.yscale('log')

plt.show()