import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

R = 8.31446261815324e7
kb = 1.380649e-16
amu = 1.66053906660e-24

with open("tracers.txt", "r", encoding="utf-8") as f:
    nsp = int(f.readline().split()[0])
    species = f.readline().split()
    rho_d = np.array([float(x) for x in f.readline().split()])
    mol_w_sp = np.array([float(x) for x in f.readline().split()])
    r_med = np.array([float(x) for x in f.readline().split()])
    sigma, dist = f.readline().split()
    sigma = float(sigma)
    dist = int(dist)

data = np.loadtxt("tracers.txt", skiprows=6)

Tl = data[:, 2]
pl = data[:, 3] / 1e5
grav = data[:, 4]
mu = data[:, 5]
q_v = data[:, 8:8 + nsp]
q_c = data[:, 8 + nsp:8 + 2 * nsp]
vf = data[:, 8 + 2 * nsp:8 + 3 * nsp]

rho = (pl * 1e6 * mu * amu) / (kb * Tl)

def p_vap_species(name, temp, pressure_dyn):
    if name == "TiO2":
        return np.exp(
            -7.70443e4 / temp
            + 4.03144e1
            - 2.59140e-3 * temp
            + 6.02422e-7 * temp**2
            - 6.86899e-11 * temp**3
        )
    if name == "Al2O3":
        return 10.0 ** (17.7 - 45892.6 / temp) * 1.0e6
    if name == "Fe":
        return 10.0 ** (7.23 - 20995.0 / temp) * 1.0e6
    if name == "Mg2SiO4":
        return 10.0 ** (14.88 - 32488.0 / temp - 0.2 * np.log10(pressure_dyn / 1e6)) * 1.0e6
    raise ValueError(f"unsupported species {name}")

q_s = np.zeros((len(pl), nsp))
for j, name in enumerate(species):
    rd_v = R / mol_w_sp[j]
    q_s[:, j] = (p_vap_species(name, Tl, pl * 1e6) / (rd_v * Tl)) / rho

N_c = np.zeros((len(pl), nsp))
for j in range(nsp):
    N_c[:, j] = (3.0 * q_c[:, j] * rho) / (4.0 * np.pi * rho_d[j] * r_med[j] ** 3)
    if dist == 1:
        N_c[:, j] *= np.exp(-9.0 / 2.0 * np.log(sigma) ** 2)

wl = np.loadtxt("opac.txt", max_rows=1)
opac = np.loadtxt("opac.txt", skiprows=1)
nwl = len(wl)

k_ext = opac[:, 2:2 + nwl]
ssa = opac[:, 2 + nwl:2 + 2 * nwl]
g = opac[:, 2 + 2 * nwl:2 + 3 * nwl]

col = sns.color_palette("colorblind")

fig, ax = plt.subplots()
for j, name in enumerate(species):
    ax.plot(q_v[:, j], pl, c=col[j % len(col)], ls="solid", label=fr"{name} $q_v$")
    ax.plot(q_c[:, j], pl, c=col[j % len(col)], ls="dashed", label=fr"{name} $q_c$")
    ax.plot(q_s[:, j], pl, c=col[j % len(col)], ls="dotted", label=fr"{name} $q_s$")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("mass mixing ratio")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()
ax.legend()

fig, ax = plt.subplots()
for j, name in enumerate(species):
    ax.plot(N_c[:, j], pl, c=col[j % len(col)], label=name)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$N_c$ [cm$^{-3}$]")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()

fig, ax = plt.subplots()
for j, name in enumerate(species):
    ax.axvline(r_med[j] * 1e4, c=col[j % len(col)], label=name)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$r_\mathrm{med}$ [$\mu$m]")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()

fig, ax = plt.subplots()
for j, name in enumerate(species):
    ax.plot(vf[:, j], pl, c=col[j % len(col)], label=name)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$v_f$ [cm s$^{-1}$]")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()

fig, ax = plt.subplots()
for i in range(nwl):
    ax.plot(k_ext[:, i], pl, label=f"{wl[i]:.2f}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\kappa_\mathrm{ext}$ [cm$^2$ g$^{-1}$]")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()
ax.legend(ncol=2, fontsize=8)

fig, ax = plt.subplots()
for i in range(nwl):
    ax.plot(ssa[:, i], pl, label=f"{wl[i]:.2f}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("single-scattering albedo")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()

fig, ax = plt.subplots()
for i in range(nwl):
    ax.plot(g[:, i], pl, label=f"{wl[i]:.2f}")
ax.set_xscale("linear")
ax.set_yscale("log")
ax.set_xlabel("asymmetry factor")
ax.set_ylabel("pressure [bar]")
ax.invert_yaxis()

plt.show()
