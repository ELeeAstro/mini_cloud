#!/usr/bin/env python3
"""Plot opacity outputs for the mixed-grain examples.

This script is layout-independent because it reads opac_k/a/g files, not the
tracer moment columns. Included here for completeness with the updated plotting
set.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from common import finish_figure, load_opacity, setup_pressure_axis


def _save_path(save_dir: Path | None, name: str) -> Path | None:
    return None if save_dir is None else save_dir / name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="opac", help="Prefix for opac_k/a/g.txt")
    parser.add_argument("--save-dir", default=None)
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    op = load_opacity(args.prefix)
    p = op["pl"]
    wl = op["wl"]
    save_dir = Path(args.save_dir) if args.save_dir else None
    show = not args.no_show
    cmap = plt.get_cmap("viridis")

    for key, xlabel, fname, logx in [
        ("k_ext", r"$\kappa_{\rm ext}$ [cm$^2$ g$^{-1}$]", "opacity_extinction.png", True),
        ("ssa", r"single-scattering albedo", "opacity_ssa.png", False),
        ("g", r"asymmetry parameter $g$", "opacity_g.png", False),
    ]:
        fig, ax = plt.subplots()
        arr = op[key]
        for i in range(wl.size):
            colour = cmap(i / max(wl.size - 1, 1))
            ax.plot(arr[:, i], p, color=colour, label=f"{wl[i]:.2f}")
        setup_pressure_axis(ax, p)
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.legend(title=r"$\lambda$ [$\mu$m]", fontsize=8)
        finish_figure(fig, _save_path(save_dir, fname), show)


if __name__ == "__main__":
    main()
