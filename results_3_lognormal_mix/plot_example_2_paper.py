#!/usr/bin/env python3
"""Paper-style radius/number-density plot for the KCl/ZnS mixed-grain case.

Updated for scalar total q_2 and additive-volume mixed density.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from common import compute_bulk, finish_figure, load_tracers


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="tracers.txt", help="Path to tracers.txt")
    parser.add_argument("--output", default=None, help="Optional output PNG/PDF path")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    tr = load_tracers(args.file)
    if tr["ndust"] != 2:
        raise ValueError(f"plot_example_2_paper.py expects ndust=2, got ndust={tr['ndust']}")

    bulk = compute_bulk(tr)
    p = tr["p_bar"]
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()

    p_rc, = ax1.plot(bulk["r_c"], p, color=colours[0], label=r"$r_{\rm c}$")
    p_nc, = ax2.plot(bulk["N_c"], p, color=colours[1], label=r"$N_{\rm c}$", linestyle="--")

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.invert_yaxis()

    yticks = [100, 10, 1, 0.1, 0.01, 1.0e-3]
    yticks_lab = ["100", "10", "1", "0.1", "0.01", r"10$^{-3}$"]
    ax1.set_yticks(yticks, yticks_lab)
    ax1.set_ylim(300, 3.0e-3)

    ax1.tick_params(axis="both", which="major", labelsize=14)
    ax2.tick_params(axis="both", which="major", labelsize=14)
    ax1.set_xlabel(r"$r_{\rm c}$ [$\mu$m]", fontsize=16)
    ax2.set_xlabel(r"$N_{\rm c}$ [cm$^{-3}$]", fontsize=16)
    ax1.set_ylabel(r"$p_{\rm gas}$ [bar]", fontsize=16)

    ax2.legend([p_rc, p_nc], [p_rc.get_label(), p_nc.get_label()], fontsize=10, loc="upper left")
    finish_figure(fig, args.output, show=not args.no_show)


if __name__ == "__main__":
    main()
