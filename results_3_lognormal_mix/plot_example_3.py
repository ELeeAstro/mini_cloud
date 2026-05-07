#!/usr/bin/env python3
"""Diagnostic plots for the four-species TiO2/Al2O3/Fe/Mg2SiO4 mixed-grain case.

Updated for the scalar total q_2 layout produced by the single-distribution
lognormal mixed-grain model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import (
    compute_bulk,
    equilibrium_mass_fraction,
    finish_figure,
    load_tracers,
    setup_pressure_axis,
    species_label,
    vapour_mol_weights,
)


def _save_path(save_dir: Path | None, name: str) -> Path | None:
    return None if save_dir is None else save_dir / name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", default="tracers.txt", help="Path to tracers.txt")
    parser.add_argument("--save-dir", default=None, help="Directory for PNG output")
    parser.add_argument("--no-show", action="store_true", help="Save/construct plots without opening windows")
    args = parser.parse_args()

    save_dir = Path(args.save_dir) if args.save_dir else None
    show = not args.no_show

    tr = load_tracers(args.file)
    if tr["ndust"] != 4:
        raise ValueError(f"plot_example_3.py expects ndust=4, got ndust={tr['ndust']}")

    bulk = compute_bulk(tr)
    p = tr["p_bar"]
    species = tr["species"]
    q_s = equilibrium_mass_fraction(
        species,
        tr["T"],
        tr["p_bar"],
        bulk["rho"],
        vapour_mol_weights(species, tr["mol_w_sp"]),
    )

    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots()
    ax.plot(tr["dTdt"], p, color=colours[0])
    setup_pressure_axis(ax, p)
    ax.set_xlabel(r"$dT/dt$ [K s$^{-1}$]")
    finish_figure(fig, _save_path(save_dir, "case3_dTdt.png"), show)

    fig, ax = plt.subplots()
    for j, name in enumerate(species):
        ax.plot(bulk["m_mix"][:, j], p, label=species_label(name), color=colours[j])
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel("mass fraction")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "case3_mass_fraction.png"), show)

    fig, ax = plt.subplots()
    for j, name in enumerate(species):
        ax.plot(bulk["V_mix"][:, j], p, label=species_label(name), color=colours[j])
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel("volume fraction")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "case3_volume_fraction.png"), show)

    fig, ax = plt.subplots()
    vf = tr["v_f"]
    ax.plot(vf[:, 0], p, color=colours[0], label=r"$v_{\rm f,0}$")
    for j, name in enumerate(species):
        ax.plot(vf[:, 1 + j], p, color=colours[1 + j], label=r"$v_{\rm f,1}$ " + species_label(name))
    ax.plot(vf[:, 1 + len(species)], p, color=colours[3], linestyle="--", label=r"$v_{\rm f,2}$")
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel(r"settling velocity [cm s$^{-1}$]")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "case3_settling_velocity.png"), show)

    fig, ax = plt.subplots()
    ax.plot(bulk["N_c"], p, color=colours[0], label=r"$N_{\rm c}$")
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel(r"$N_{\rm c}$ [cm$^{-3}$]")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "case3_number_density.png"), show)

    fig, ax = plt.subplots()
    ax.plot(bulk["r_c"], p, color=colours[1], label=r"$r_{\rm c}$")
    if "r_med" in bulk:
        ax.plot(bulk["r_med"], p, color=colours[2], linestyle="--", label=r"$r_{\rm med}$")
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel(r"radius [$\mu$m]")
    ax.legend()
    finish_figure(fig, _save_path(save_dir, "case3_radius.png"), show)

    fig, ax = plt.subplots()
    for j, name in enumerate(species):
        label = species_label(name)
        ax.plot(tr["q_v"][:, j], p, color=colours[j], label=rf"{label} $q_{{\rm v}}$")
        ax.plot(tr["q_1"][:, j], p, color=colours[j], linestyle="--", label=rf"{label} $q_1$")
        ax.plot(q_s[:, j], p, color=colours[j], linestyle=":", label=rf"{label} $q_{{\rm s}}$")
    ax.plot(tr["q_0"], p, color=colours[4], linestyle="-.", label=r"$q_0$")
    ax.plot(tr["q_2"], p, color=colours[5], linestyle=(0, (3, 1, 1, 1)), label=r"$q_{2,{\rm tot}}$")
    setup_pressure_axis(ax, p)
    ax.set_xscale("log")
    ax.set_xlabel("mixing ratio / moment variable")
    ax.legend(fontsize=8)
    finish_figure(fig, _save_path(save_dir, "case3_tracers.png"), show)

    if "sig_g" in bulk:
        fig, ax = plt.subplots()
        ax.plot(bulk["sig_g"], p, color=colours[4], label=r"$\sigma_g$")
        setup_pressure_axis(ax, p)
        ax.set_xscale("log")
        ax.set_xlabel(r"$\sigma_g$")
        ax.legend()
        finish_figure(fig, _save_path(save_dir, "case3_sigma_g.png"), show)


if __name__ == "__main__":
    main()
