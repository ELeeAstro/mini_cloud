import argparse
import os

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from common import load_tracers, compute_bulk, vapour_mol_weights, equilibrium_mass_fraction


def main():
    parser = argparse.ArgumentParser(description="Plot scalar-q2 gamma mixed-cloud results.")
    parser.add_argument("file", nargs="?", default=None)
    parser.add_argument("--dir", default=".", help="Directory containing tracers.txt")
    parser.add_argument("--met", type=float, default=0.0, help="Metallicity used in some saturation fits")
    args = parser.parse_args()

    path = args.file or os.path.join(args.dir, "tracers.txt")
    tr = load_tracers(path)
    bulk = compute_bulk(tr)

    T = tr["T"]
    pl = tr["P_bar"]
    q_v = tr["q_v"]
    q_0 = tr["q_0"]
    q_1 = tr["q_1"]
    q_2 = tr["q_2"]
    dTdt = tr["dTdt"]
    vf = tr["v_f"]
    ndust = tr["ndust"]
    sp_tex = tr["species_tex"]

    mol_w_v = vapour_mol_weights(tr["species"], tr["mol_w_sp"])
    q_s = equilibrium_mass_fraction(tr["species"], T, pl, bulk["rho"], mol_w_v, met=args.met)

    col = sns.color_palette("colorblind", max(8, ndust + 3))

    # --- Mixing ratios ---
    plt.figure()
    for j in range(ndust):
        plt.plot(q_v[:, j], pl, c=col[j], label=sp_tex[j] + r" $q_{\rm v}$")
        plt.plot(q_1[:, j], pl, c=col[j], ls="dashed", label=sp_tex[j] + r" $q_1$")
        if np.any(np.isfinite(q_s[:, j])):
            plt.plot(q_s[:, j], pl, c=col[j], ls="dotted", label=sp_tex[j] + r" $q_{\rm s}$")
    plt.plot(q_0, pl, c=col[ndust], ls="dashdot", label=r"$q_0$")
    plt.plot(q_2, pl, c=col[ndust + 1], ls=":", label=r"$q_2$ total")
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$q$")
    plt.ylabel(r"$p$ [bar]")
    plt.legend(fontsize=8)

    # --- Number density and radius ---
    plt.figure()
    plt.plot(bulk["N_c"], pl, c=col[0], label=r"$N_{\rm c}$")
    plt.plot(bulk["r_c"], pl, c=col[1], label=r"$r_{\rm c}$ [$\mu$m]")
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$N_{\rm c}$ [cm$^{-3}$] / $r_{\rm c}$ [$\mu$m]")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    # --- Volume fractions ---
    plt.figure()
    for j in range(ndust):
        plt.plot(bulk["V_mix"][:, j], pl, c=col[j], label=sp_tex[j])
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$V_{\rm mix}$")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    # --- Settling velocities ---
    plt.figure()
    plt.plot(vf[:, 0], pl, c=col[0], label=r"$v_{{\rm f},0}$")
    for j in range(ndust):
        plt.plot(vf[:, 1 + j], pl, c=col[j + 1], label=sp_tex[j] + r" $q_1$")
    plt.plot(vf[:, ndust + 1], pl, c=col[ndust + 1], ls="dashed", label=r"$q_2$ total")
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$v_{\rm f}$ [cm s$^{-1}$]")
    plt.ylabel(r"$p$ [bar]")
    plt.legend(fontsize=8)

    # --- Gamma shape parameter ---
    plt.figure()
    plt.plot(bulk["nu"], pl, c=col[0], label=r"$\nu$")
    plt.plot(bulk["lam"], pl, c=col[1], label=r"$\lambda$ [g]")
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$\nu$ / $\lambda$")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    if dTdt is not None:
        plt.figure()
        plt.plot(dTdt, pl, c=col[0], label=r"$dT/dt$")
        plt.yscale("log")
        plt.gca().invert_yaxis()
        plt.xlabel(r"$dT/dt$ [K s$^{-1}$]")
        plt.ylabel(r"$p$ [bar]")
        plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
