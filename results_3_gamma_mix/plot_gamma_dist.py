import argparse
import os

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.special import gammaln

from common import load_tracers, compute_bulk


def gamma_distribution(tr, bulk, n_radius=1000, r_min_um=1.0e-3, r_max_um=1.0e2):
    r_cm = np.logspace(np.log10(r_min_um), np.log10(r_max_um), n_radius) * 1.0e-4
    r_um = r_cm * 1.0e4

    rho_d_m = bulk["rho_d_m"]
    N_c = bulk["N_c"]
    nu = bulk["nu"]
    lam = bulk["lam"]
    rho_c_t = bulk["rho_c_t"]

    m = 4.0 / 3.0 * np.pi * r_cm[:, None] ** 3 * rho_d_m[None, :]

    valid = (N_c > 0.0) & (rho_c_t > 0.0) & (nu > 0.0) & (lam > 0.0)
    mf = np.full_like(m, np.nan)

    if np.any(valid):
        log_m = np.log(np.maximum(m[:, valid], 1.0e-300))
        nu_v = nu[valid][None, :]
        lam_v = lam[valid][None, :]
        log_fx = (
            np.log(np.maximum(N_c[valid], 1.0e-300))[None, :]
            - nu_v * np.log(lam_v)
            - gammaln(nu_v)
            + (nu_v - 1.0) * log_m
            - m[:, valid] / lam_v
            + log_m
        )
        mf[:, valid] = np.exp(log_fx)

    return r_um, mf, valid


def main():
    parser = argparse.ArgumentParser(
        description="Plot gamma particle mass distributions from scalar-q2 mixed-cloud output."
    )
    parser.add_argument("file", nargs="?", default=None)
    parser.add_argument("--dir", default=".", help="Directory containing tracers.txt")
    parser.add_argument("--stride", type=int, default=5, help="Plot every Nth pressure layer")
    parser.add_argument("--output", "-o", default=None)
    parser.add_argument("--rmin", type=float, default=1.0e-3, help="Minimum radius [micron]")
    parser.add_argument("--rmax", type=float, default=1.0e2, help="Maximum radius [micron]")
    parser.add_argument("--ymin", type=float, default=1.0e-3)
    parser.add_argument("--ymax", type=float, default=1.0e1)
    args = parser.parse_args()

    path = args.file or os.path.join(args.dir, "tracers.txt")
    tr = load_tracers(path)
    bulk = compute_bulk(tr, nu_max=10.0)
    pl = tr["P_bar"]

    r_um, mf, valid = gamma_distribution(
        tr,
        bulk,
        n_radius=1000,
        r_min_um=args.rmin,
        r_max_um=args.rmax,
    )

    fig, ax = plt.subplots()
    logp = np.log10(pl)
    normalize = mcolors.Normalize(vmin=np.nanmin(logp), vmax=np.nanmax(logp))
    cmap = sns.color_palette("crest", as_cmap=True)

    for i in range(0, len(pl), args.stride):
        if valid[i] and np.any(np.isfinite(mf[:, i])):
            ax.plot(r_um, mf[:, i], c=cmap(normalize(np.log10(pl[i]))))

    scalarmappable = cm.ScalarMappable(norm=normalize, cmap=cmap)
    cbar = fig.colorbar(scalarmappable, ax=ax)
    cticks = np.linspace(np.nanmin(logp), np.nanmax(logp), 10)
    cbar.set_ticks(cticks)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r"$\log_{10}$ $p$ [bar]", fontsize=14)

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlabel(r"$r$ [$\mu$m]", fontsize=16)
    ax.set_ylabel(r"$m \cdot f(m)$ [cm$^{-3}$]", fontsize=16)
    ax.set_ylim(args.ymin, args.ymax)
    ax.set_xlim(args.rmin, args.rmax)

    fig.tight_layout(pad=1.05)
    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
