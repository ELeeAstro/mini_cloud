#!/usr/bin/env python3
"""Generic diagnostic plotting script for scalar-q_2 mixed-grain results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from common import compute_bulk, finish_figure, load_tracers, setup_pressure_axis, species_label


def _save_path(save_dir: Path | None, name: str) -> Path | None:
    return None if save_dir is None else save_dir / name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="tracers.txt", help="Path to tracers.txt")
    parser.add_argument("--save-dir", default=None, help="Directory for PNG output")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args()

    tr = load_tracers(args.file)
    bulk = compute_bulk(tr)
    p = tr["p_bar"]
    species = tr["species"]
    save_dir = Path(args.save_dir) if args.save_dir else None
    show = not args.no_show
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots()
    for j, name in enumerate(species):
        c = colours[j % len(colours)]
        ax.plot(tr["q_v"][:, j], p, color=c, label=rf"{species_label(name)} $q_{{\rm v}}$")
        ax.plot(tr["q_1"][:, j], p, color=c, linestyle="--", label=rf"{species_label(name)} $q_1$")
    ax.plot(tr["q_0"], p, color=colours[len(species) % len(colours)], linestyle="-.", label=r"$q_0$")
    if tr["q_2"] is not None:
        ax.plot(tr["q_2"], p, color=colours[(len(species)+1) % len(colours)], linestyle=":", label=r"$q_{2,\rm tot}$")
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel("mixing ratio / moment variable")
    ax.legend(fontsize=8)
    finish_figure(fig, _save_path(save_dir, "generic_tracers.png"), show)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twiny()
    p1, = ax1.plot(bulk["r_c"], p, color=colours[0], label=r"$r_{\rm c}$")
    p2, = ax2.plot(bulk["N_c"], p, color=colours[1], linestyle="--", label=r"$N_{\rm c}$")
    setup_pressure_axis(ax1, p)
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.set_xlabel(r"$r_{\rm c}$ [$\mu$m]")
    ax2.set_xlabel(r"$N_{\rm c}$ [cm$^{-3}$]")
    ax1.legend([p1, p2], [p1.get_label(), p2.get_label()], loc="best")
    finish_figure(fig, _save_path(save_dir, "generic_radius_number.png"), show)

    fig, ax = plt.subplots()
    for j, name in enumerate(species):
        ax.plot(bulk["V_mix"][:, j], p, color=colours[j % len(colours)], label=species_label(name))
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel("volume fraction")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "generic_volume_fraction.png"), show)

    if "sig_g" in bulk:
        fig, ax = plt.subplots()
        ax.plot(bulk["sig_g"], p, color=colours[0], label=r"$\sigma_g$")
        setup_pressure_axis(ax, p)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\sigma_g$")
        ax.legend()
        finish_figure(fig, _save_path(save_dir, "generic_sigma_g.png"), show)

    if tr["dTdt"] is not None:
        fig, ax = plt.subplots()
        ax.plot(tr["dTdt"], p, color=colours[0], label=r"$dT/dt$")
        setup_pressure_axis(ax, p)
        ax.set_xlabel(r"$dT/dt$ [K s$^{-1}$]")
        ax.legend()
        finish_figure(fig, _save_path(save_dir, "generic_dTdt.png"), show)


if __name__ == "__main__":
    main()
