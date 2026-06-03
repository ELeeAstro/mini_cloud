#!/usr/bin/env python3
"""Derive and plot latent heat from mini_cloud vapour-pressure fits.

Run with the Miniforge Python, for example:

    /Users/gl20y334/miniforge3/bin/python tools/latent_heat_from_pvap.py
    /Users/gl20y334/miniforge3/bin/python tools/latent_heat_from_pvap.py --species TiO2 H2O CO2 --tmin 50 --tmax 2000

The full expression is calculated from Clausius-Clapeyron:

    L = (R_gas / mol_w) * T**2 * d ln(p_vap) / dT

The linear expression keeps only the leading -A/T term in ln(p_vap), which
matches the approximation used by the current Fortran l_heat_sp routines.
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mini_cloud_matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


T = sp.symbols("T", positive=True)
R_GAS = 8.31446261815324e7  # erg mol-1 K-1
LOG10 = sp.log(10)


@dataclass(frozen=True)
class SpeciesFit:
    name: str
    mol_w: float
    ln_pvap: sp.Expr
    linear_a: sp.Expr
    note: str = ""

    @property
    def full_latent(self) -> sp.Expr:
        return sp.simplify((R_GAS / self.mol_w) * T**2 * sp.diff(self.ln_pvap, T))

    @property
    def linear_latent(self) -> sp.Expr:
        return sp.simplify((R_GAS / self.mol_w) * self.linear_a)


def exp_fit(name: str, mol_w: float, ln_pvap: sp.Expr, a: sp.Expr, note: str = "") -> SpeciesFit:
    return SpeciesFit(name=name, mol_w=mol_w, ln_pvap=ln_pvap, linear_a=a, note=note)


def log10_fit(name: str, mol_w: float, log10_pvap: sp.Expr, a: sp.Expr, note: str = "") -> SpeciesFit:
    return SpeciesFit(name=name, mol_w=mol_w, ln_pvap=LOG10 * log10_pvap, linear_a=LOG10 * a, note=note)


def build_species() -> dict[str, SpeciesFit]:
    mk_liquid_inner = (
        53.878
        - 1331.22 / T
        - 9.44523 * sp.log(T)
        + 0.014025 * T
    )
    mk_liquid = (
        54.842763
        - 6763.22 / T
        - 4.210 * sp.log(T)
        + 0.000367 * T
        + sp.tanh(0.0415 * (T - 218.8)) * mk_liquid_inner
    )
    h2o_ice = 9.550426 - 5723.265 / T + 3.53068 * sp.log(T) - 0.00728332 * T

    return {
        "C": exp_fit("C", 12.0107, 32.7860 - 8.65139e4 / (T + 4.80395e-1), 8.65139e4 * (T / (T + 4.80395e-1)) ** 2),
        "TiC": log10_fit("TiC", 59.8777, -33600.0 / T + 7.652, 33600.0),
        "SiC": exp_fit("SiC", 40.0962, -9.51431385e4 / T + 37.2019157 + 1.09809718e-3 * T - 5.63629542e-7 * T**2 + 6.97886017e-11 * T**3, 9.51431385e4),
        "CaTiO3": log10_fit("CaTiO3", 135.943, -72160.0 / T + 30.24, 72160.0, "pressure and metallicity terms are omitted because they differentiate to zero with respect to T"),
        "TiO2": exp_fit("TiO2", 79.866, -7.70443e4 / T + 40.3144 - 2.59140e-3 * T + 6.02422e-7 * T**2 - 6.86899e-11 * T**3, 7.70443e4),
        "VO": exp_fit("VO", 66.94090, -6.74603e4 / T + 38.2717 - 2.78551e-3 * T + 5.72078e-7 * T**2 - 7.41840e-11 * T**3, 6.74603e4),
        "Al2O3": log10_fit("Al2O3", 101.961, 17.7 - 45892.6 / T, 45892.6, "metallicity term is omitted because it differentiates to zero with respect to T"),
        "Fe": log10_fit("Fe", 55.845, 7.23 - 20995.0 / T, 20995.0),
        "FeS": exp_fit("FeS", 87.91, -5.69922e4 / T + 38.6753 - 4.68301e-3 * T + 1.03559e-6 * T**2 - 8.42872e-11 * T**3, 5.69922e4),
        "FeO": exp_fit("FeO", 71.8444, -6.30018e4 / T + 36.6364 - 2.42990e-3 * T + 3.18636e-7 * T**2, 6.30018e4),
        "Mg2SiO4": log10_fit("Mg2SiO4", 140.693, 14.88 - 32488.0 / T, 32488.0, "pressure and metallicity terms are omitted because they differentiate to zero with respect to T"),
        "MgSiO3": log10_fit("MgSiO3", 100.389, 13.43 - 28665.0 / T, 28665.0, "metallicity term is omitted because it differentiates to zero with respect to T"),
        "MgSiO3_amorph": log10_fit("MgSiO3_amorph", 100.389, 13.43 - 28665.0 / T, 28665.0, "same vapour-pressure fit as MgSiO3"),
        "MgO": exp_fit("MgO", 40.3044, -8.44591675e4 / T + 45.8952742 + 4.38797163e-3 * T - 9.84437451e-8 * T**2 - 7.30855823e-11 * T**3, 8.44591675e4),
        "SiO2": exp_fit("SiO2", 60.084, -7.28086e4 / T + 36.5312 - 2.56109e-4 * T - 5.24980e-7 * T**2 + 1.53343e-10 * T**3, 7.28086e4),
        "SiO2_amorph": exp_fit("SiO2_amorph", 60.084, -7.28086e4 / T + 36.5312 - 2.56109e-4 * T - 5.24980e-7 * T**2 + 1.53343e-10 * T**3, 7.28086e4, "same vapour-pressure fit as SiO2"),
        "SiO": exp_fit("SiO", 44.085, -49520.0 / T + 32.52, 49520.0),
        "Cr": exp_fit("Cr", 51.996, -4.78455e4 / T + 32.2423 - 5.28710e-4 * T - 6.17347e-8 * T**2 + 2.88469e-12 * T**3, 4.78455e4),
        "MnS": log10_fit("MnS", 87.003, 11.532 - 23810.0 / T, 23810.0, "metallicity term is omitted because it differentiates to zero with respect to T"),
        "Na2S": log10_fit("Na2S", 78.0445, 8.550 - 13889.0 / T, 13889.0, "metallicity term is omitted because it differentiates to zero with respect to T"),
        "ZnS": exp_fit("ZnS", 97.445, -4.75507888e4 / T + 36.6993865 - 2.49490016e-3 * T + 7.29116854e-7 * T**2 - 1.12734453e-10 * T**3, 4.75507888e4),
        "KCl": exp_fit("KCl", 74.551, -2.69250e4 / T + 33.9574 - 2.04903e-3 * T - 2.83957e-7 * T**2 + 1.82974e-10 * T**3, 2.69250e4),
        "NaCl": exp_fit("NaCl", 58.443, -2.79146e4 / T + 34.6023 - 3.11287e3 * T + 5.30965e-7 * T**2 - 2.59584e-12 * T**3, 2.79146e4),
        "S2": exp_fit("S2", 64.132, sp.Piecewise((27.0 - 18500.0 / T, T < 413.0), (16.1 - 14000.0 / T, True)), sp.Piecewise((18500.0, T < 413.0), (14000.0, True))),
        "S8": exp_fit("S8", 256.528, sp.Piecewise((20.0 - 11800.0 / T, T < 413.0), (9.6 - 7510.0 / T, True)), sp.Piecewise((11800.0, T < 413.0), (7510.0, True))),
        "NH4Cl": log10_fit("NH4Cl", 53.4915, 7.0220 - 4302.0 / T, 4302.0),
        "H2O": exp_fit("H2O", 18.015, sp.Piecewise((h2o_ice, T <= 273.16), (mk_liquid, T < 1048.0), (mk_liquid.subs(T, 1048.0), True)), sp.Piecewise((5723.265, T <= 273.16), (6763.22 + 1331.22 * sp.tanh(0.0415 * (T - 218.8)), T < 1048.0), (0.0, True))),
        "NH3": exp_fit("NH3", 17.031, 15.96 - 3537.0 / T - 3.310e4 / T**2 + 1.742e6 / T**3 - 2.995e7 / T**4, 3537.0),
        "CH4": exp_fit("CH4", 16.043, 10.51 - 1.110e3 / T - 4.341e3 / T**2 + 1.035e5 / T**3 - 7.910e5 / T**4, 1.110e3),
        "NH4SH": log10_fit("NH4SH", 51.1114, 7.8974 - 2409.4 / T, 2409.4),
        "H2S": exp_fit("H2S", 34.08, 12.98 - 2.707e3 / T, 2.707e3),
        "H2SO4": exp_fit("H2SO4", 98.0785, -1.01294e4 / T + 35.5465 - 8.34848e-3 * T, 1.01294e4),
        "CO": exp_fit("CO", 28.0101, sp.Piecewise((10.43 - 721.3 / T - 1.074e4 / T**2 + 2.341e5 / T**3 - 2.392e6 / T**4 + 9.478e6 / T**5, T < 61.55), (10.25 - 748.2 / T - 5.843e3 / T**2 + 3.939e4 / T**3, True)), sp.Piecewise((721.3, T < 61.55), (748.2, True))),
        "CO2": exp_fit("CO2", 44.0095, sp.Piecewise((14.76 - 2.571e3 / T - 7.781e4 / T**2 + 4.325e6 / T**3 - 1.207e8 / T**4 + 1.350e9 / T**5, T < 194.7), (18.61 - 4.154e3 / T + 1.041e5 / T**2, True)), sp.Piecewise((2.571e3, T < 194.7), (4.154e3, True))),
        "O2": exp_fit("O2", 31.9988, 15.29 - 1166.2 / T - 0.75587 * sp.log(T) + 0.14188 * T - 1.8665e-3 * T**2 + 7.582e-6 * T**3, 1166.2),
    }


def expression_text(expr: sp.Expr) -> str:
    return sp.sstr(sp.expand(sp.simplify(expr)), full_prec=False)


def print_expressions(species: list[SpeciesFit]) -> None:
    for fit in species:
        print(f"\n{fit.name}")
        print(f"  molecular_weight = {fit.mol_w:g} g mol-1")
        if fit.note:
            print(f"  note = {fit.note}")
        print("  ln_pvap(T) =")
        print(f"    {expression_text(fit.ln_pvap)}")
        print("  linear L(T) [erg g-1] =")
        print(f"    {expression_text(fit.linear_latent)}")
        print("  full L(T) [erg g-1] =")
        print(f"    {expression_text(fit.full_latent)}")


def finite_values(func, t_grid: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        vals = np.asarray(func(t_grid), dtype=float)
    if vals.ndim == 0:
        vals = np.full_like(t_grid, float(vals))
    vals[~np.isfinite(vals)] = np.nan
    return vals


def plot_species(species: list[SpeciesFit], tmin: float, tmax: float, ntemp: int, cols: int) -> None:
    t_grid = np.linspace(tmin, tmax, ntemp)
    rows = math.ceil(len(species) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(5.2 * cols, 3.5 * rows), squeeze=False)

    for ax, fit in zip(axes.flat, species):
        full_func = sp.lambdify(T, fit.full_latent, "numpy")
        linear_func = sp.lambdify(T, fit.linear_latent, "numpy")
        full = finite_values(full_func, t_grid)
        linear = finite_values(linear_func, t_grid)

        ax.plot(t_grid, linear, label="linear", lw=2)
        ax.plot(t_grid, full, label="full", lw=1.6, ls="--")
        ax.set_title(fit.name)
        ax.set_xlabel("T [K]")
        ax.set_ylabel("L [erg g$^{-1}$]")
        ax.grid(alpha=0.25)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    for ax in axes.flat[len(species):]:
        ax.set_visible(False)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right")
    fig.tight_layout(rect=(0, 0, 0.97, 1))
    plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", nargs="+", default=["TiO2", "Al2O3", "Fe", "Mg2SiO4", "MgSiO3", "SiO2", "H2O", "NH3", "CH4", "CO", "CO2", "O2"], help="Species to print and plot, or 'all'.")
    parser.add_argument("--tmin", type=float, default=50.0, help="Minimum plot temperature [K].")
    parser.add_argument("--tmax", type=float, default=2500.0, help="Maximum plot temperature [K].")
    parser.add_argument("--ntemp", type=int, default=500, help="Number of temperature samples.")
    parser.add_argument("--cols", type=int, default=3, help="Number of subplot columns.")
    parser.add_argument("--no-plot", action="store_true", help="Print expressions only.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fits = build_species()

    if args.species == ["all"]:
        selected_names = sorted(fits)
    else:
        selected_names = args.species

    missing = [name for name in selected_names if name not in fits]
    if missing:
        known = ", ".join(sorted(fits))
        raise SystemExit(f"Unknown species: {', '.join(missing)}\nKnown species: {known}")

    selected = [fits[name] for name in selected_names]
    print_expressions(selected)

    if not args.no_plot:
        plot_species(selected, args.tmin, args.tmax, args.ntemp, max(1, args.cols))


if __name__ == "__main__":
    main()
