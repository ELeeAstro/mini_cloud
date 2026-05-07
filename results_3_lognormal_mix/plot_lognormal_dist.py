#!/usr/bin/env python3
"""Plot the reconstructed lognormal mass/radius distribution.

Updated for scalar total q_2 and additive-volume mixed density. Works for both
case 2 and case 3 because the species count is inferred from the tracer header.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

from common import R_SEED, compute_bulk, finish_figure, load_tracers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="tracers.txt", help="Path to tracers.txt")
    parser.add_argument("--output", default=None, help="Optional output PNG/PDF path")
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--stride", type=int, default=5, help="Plot every Nth pressure level")
    parser.add_argument("--r-max-micron", type=float, default=100.0)
    parser.add_argument("--n-grid", type=int, default=1000)
    args = parser.parse_args()

    tr = load_tracers(args.file)
    bulk = compute_bulk(tr)
    if "sig_g" not in bulk:
        raise ValueError("No q_2 column was found, so the lognormal width cannot be reconstructed.")

    p = np.asarray(tr["p_bar"])
    rho_d_m = np.asarray(bulk["rho_d_m"])
    N_c = np.asarray(bulk["N_c"])
    m_c = np.asarray(bulk["m_c"])
    lnsig2 = np.asarray(bulk["lnsig2"])
    lnsig = np.sqrt(np.maximum(lnsig2, 1.0e-30))
    m_med = np.asarray(bulk["m_med"])

    # Mass grid based on the layer-wise mixed density range. This avoids assuming
    # a single material density for all plotted radii.
    r_min = R_SEED
    r_max = args.r_max_micron * 1.0e-4
    rho_min = np.nanmin(rho_d_m[np.isfinite(rho_d_m) & (rho_d_m > 0.0)])
    rho_max = np.nanmax(rho_d_m[np.isfinite(rho_d_m) & (rho_d_m > 0.0)])
    m_min = 4.0 / 3.0 * np.pi * r_min**3 * rho_min
    m_max = 4.0 / 3.0 * np.pi * r_max**3 * rho_max
    m = np.logspace(np.log10(m_min), np.log10(m_max), args.n_grid)

    fig, ax = plt.subplots()
    norm = mcolors.Normalize(vmin=np.log10(np.nanmin(p)), vmax=np.log10(np.nanmax(p)))
    cmap = plt.get_cmap("viridis")

    for i in range(0, p.size, args.stride):
        if not np.isfinite(lnsig[i]) or lnsig[i] <= 1.0e-3:
            continue
        f_m = N_c[i] / (np.sqrt(2.0 * np.pi) * lnsig[i]) * np.exp(
            -0.5 * (np.log(m / m_med[i])) ** 2 / lnsig2[i]
        )
        r = ((3.0 * m) / (4.0 * np.pi * rho_d_m[i])) ** (1.0 / 3.0) * 1.0e4
        ax.plot(r, f_m, color=cmap(norm(np.log10(p[i]))))

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array(np.log10(p))
    cbar = fig.colorbar(sm, ax=ax)
    ticks = np.linspace(np.log10(np.nanmin(p)), np.log10(np.nanmax(p)), 8)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{10**t:.1e}" for t in ticks])
    cbar.ax.set_ylabel(r"$p_{\rm gas}$ [bar]")

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"$r$ [$\mu$m]")
    ax.set_ylabel(r"$f(m)$ [cm$^{-3}$]")
    ax.set_xlim(1.0e-3, args.r_max_micron)
    finite_positive = []
    # Preserve the old broad default while avoiding empty-looking plots for other cases.
    ax.set_ylim(1.0e-3, 1.0e1)
    finish_figure(fig, args.output, show=not args.no_show)


if __name__ == "__main__":
    main()
