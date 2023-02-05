import numpy as np
import matplotlib.pylab as plt
import matplotlib.ticker
import seaborn as sns

fname = 'HAT-P-1b_HELIOS.txt'
data = np.loadtxt(fname,skiprows=3)


T = data[:,1]
p = data[:,2]/1e6

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

## Total atmospheric pressure grid in bar
# Number of points, lowest pressure [bar], P_low, highest pressure, P_up [bar]
npoints = 100
P_low = 1e-6
P_up = 1000.0
P = np.logspace(np.log10(P_low),np.log10(P_up),npoints)

# Metallicity [M/H] (or [Fe/H]) in dex solar
met = 0.0

# Expressions for species supersaturation line (S = 1) from the references:
# Lodders & Fegley (2002),  Visscher et al. (2006, 2010), Morley et al. (2012) Wakeford et al. (2017)
# Tc_ [K] is the temperature (of condensation) where S = 1
Tc_CaTiO3 =  1.0e4 / (5.125 - 0.277 * np.log10(P) - 0.554 * met)
Tc_Al2O3 =  1.0e4 / (5.014 - 0.2179 * np.log10(P) + 2.264e-3 * np.log10(P)**2 - 0.585 * met)
Tc_Fe = 1.0e4 / (5.44 - 0.48 * np.log10(P) - 0.48 * met)
Tc_Mg2SiO4 = 1.0e4 / (5.89 - 0.37 * np.log10(P) - 0.73 * met)
Tc_MgSiO3 = 1.0e4 / (6.26 - 0.35 * np.log10(P) - 0.70 *  met)
Tc_Cr =  1.0e4 / (6.576 - 0.486 * np.log10(P) - 0.486 * met)
Tc_MnS =  1.0e4 / (7.447 - 0.42 * np.log10(P) - 0.84 * met)
Tc_Na2S = 1.0e4 / (10.045 - 0.72 * np.log10(P) - 1.08 * met)
Tc_ZnS = 1.0e4 / (12.527 - 0.63 * np.log10(P) - 1.26 * met)
Tc_KCl = 1.0e4 / (12.479 - 0.879 * np.log10(P) - 0.879 * met)
#Tc_NH4H2PO4 = 1.0e4 / ((29.99 - 0.2 * 11.0 * np.log10(P))) # Non trival metallicity dependence
Tc_H2O = 1.0e4 / (38.84 - 3.93 * met - 3.83  * np.log10(P) - 0.20 * met * np.log10(P))
Tc_NH4SH = 1.0e4 / (48.91 - 4.15 * np.log10(P) - 4.15 * met)
Tc_NH3 = 1.0e4 / (68.02 - 6.31 * np.log10(P) - 6.19 * met)
#Tc_H2S = 1.0e4 / ((86.49 - 8.54 * np.log10(P))) # No trival metallicity dependence

# Plot all species Tc_ against total atmospheric pressure
fig, ax = plt.subplots()

c = sns.color_palette('colorblind') 

plt.plot(T,p,label='HAT-P-1b',c='black',ls='dashed')

#plt.plot(Tc_CaTiO3,P,label=r'CaTiO$_{3}$',c='red')
plt.plot(Tc_Al2O3,P,label=r'Al$_{2}$O$_{3}$',c=c[0])
plt.plot(Tc_Fe,P,label=r'Fe',c=c[1])
plt.plot(Tc_Mg2SiO4,P,label=r'Mg$_{2}$SiO$_{4}$',c=c[2])
plt.plot(Tc_MgSiO3,P,label=r'MgSiO$_{3}$',c=c[3])
plt.plot(Tc_Cr,P,label=r'Cr',c=c[4])
plt.plot(Tc_MnS,P,label=r'MnS',c=c[5])
plt.plot(Tc_Na2S,P,label=r'Na$_{2}$S',c=c[6])
plt.plot(Tc_ZnS,P,label=r'ZnS',c=c[7])
plt.plot(Tc_KCl,P,label=r'KCl',c=c[8])
#plt.plot(Tc_H2O,P,label=r'H2O',c='cyan')
#plt.plot(Tc_NH4SH,P,label=r'NH$_{4}$SH',c='teal')
#plt.plot(Tc_NH3,P,label=r'NH$_{3}$',c='olive')

plt.yscale('log')
plt.ylim(1e-6,100)
yticks = [100,10,1,0.1,0.01,1e-3,1e-4,1e-5,1e-6]
yticks_lab = [r'100',r'10',r'1',r'0.1',r'0.01',r'10$^{-3}$',r'10$^{-4}$',r'10$^{-5}$',r'10$^{-6}$']
plt.yticks(yticks,yticks_lab)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.yaxis.set_ticks_position('both')

plt.xlim(500,2250)
xticks = [500,750,1000,1250,1500,1750,2000,2250]
xticks_lab = [r'500',r'750',r'1000',r'1250',r'1500',r'1750',r'2000',r'2250']
plt.xticks(xticks,xticks_lab)

plt.tick_params(axis='both',which='major',labelsize=14)

plt.gca().invert_yaxis()

plt.ylabel(r'p$_{\rm gas}$ [bar]',fontsize=14)
plt.xlabel(r'T$_{\rm gas}$ [K]',fontsize=14)


plt.legend()

plt.tight_layout(pad=1.05, h_pad=None, w_pad=None, rect=None)

plt.savefig('HAT-P-1b_HELIOS.pdf',dpi=300,bbox_inches='tight')

plt.show()
