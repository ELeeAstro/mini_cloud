"""Shared plotting/data helpers for mini_cloud 3-moment lognormal mixed-grain runs.

Updated for the single-scalar q_2 layout:

    columns after thermodynamic/background columns are
    q_v(1:ndust), q_0, q_1(1:ndust), q_2_total, v_f(0,1...,2), dTdt

The scripts also tolerate the old vector-q_2 layout, but prefer the new layout
when both would be possible.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

R_GAS = 8.31446261815324e7  # erg mol^-1 K^-1
KB = 1.380649e-16           # erg K^-1
AMU = 1.66053906660e-24     # g
R_SEED = 1.0e-7             # cm
V_SEED = 4.0 / 3.0 * np.pi * R_SEED**3
FLOOR = 1.0e-300

SPECIES_BY_MW = {
    74.551: "KCl",
    97.445: "ZnS",
    79.866: "TiO2",
    101.961: "Al2O3",
    55.845: "Fe",
    140.693: "Mg2SiO4",
}

# Vapour species used by the saturation fits in the original plotting scripts.
VAPOUR_MW = {
    "KCl": 74.5513,
    "ZnS": 97.4450,
    "TiO2": 79.8658,
    "Al2O3": 26.98153860,
    "Fe": 55.8450,
    "Mg2SiO4": 24.305,
}


@dataclass(frozen=True)
class Layout:
    nbg: int
    q2_mode: Literal["scalar", "vector", "none"]
    has_dtdt: bool
    vf_count: int


def species_label(name: str) -> str:
    return {
        "KCl": r"KCl",
        "ZnS": r"ZnS",
        "TiO2": r"TiO$_2$",
        "Al2O3": r"Al$_2$O$_3$",
        "Fe": r"Fe",
        "Mg2SiO4": r"Mg$_2$SiO$_4$",
    }.get(name, name)


def infer_species_names(mol_w_sp: np.ndarray) -> list[str]:
    species: list[str] = []
    for idx, mol in enumerate(np.atleast_1d(mol_w_sp)):
        match = None
        for ref, candidate in SPECIES_BY_MW.items():
            if np.isclose(mol, ref, rtol=0.0, atol=5.0e-4):
                match = candidate
                break
        species.append(match or f"sp{idx + 1}")
    return species


def _load_header(path: str | Path):
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        first = np.fromstring(handle.readline(), sep=" ")
    if first.size != 1 or abs(first[0] - round(first[0])) > 1.0e-8:
        raise ValueError(f"{path} does not look like a tracer file with a 3-line header.")

    ndust = int(round(first[0]))
    rho_d = np.atleast_1d(np.loadtxt(path, skiprows=1, max_rows=1)).astype(float)
    mol_w_sp = np.atleast_1d(np.loadtxt(path, skiprows=2, max_rows=1)).astype(float)
    data = np.loadtxt(path, skiprows=3)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return ndust, rho_d, mol_w_sp, data


def _infer_layout(ncols: int, ndust: int) -> Layout:
    # Prefer the new scalar-q_2 layout. The old layout is kept only for backwards
    # compatibility with existing results files.
    candidates: list[tuple[str, bool, int]] = [
        ("scalar", True, ndust + 2),
        ("scalar", False, ndust + 2),
        ("vector", True, 2 * ndust + 1),
        ("vector", False, 2 * ndust + 1),
        ("none", True, ndust + 1),
        ("none", False, ndust + 1),
    ]

    for q2_mode, has_dtdt, vf_count in candidates:
        q2_count = {"scalar": 1, "vector": ndust, "none": 0}[q2_mode]
        fixed_without_bg = 6 + ndust + 1 + ndust + q2_count + vf_count
        if has_dtdt:
            fixed_without_bg += 1
        nbg = ncols - fixed_without_bg
        if nbg >= 1:
            return Layout(nbg=nbg, q2_mode=q2_mode, has_dtdt=has_dtdt, vf_count=vf_count)

    raise ValueError(f"Could not infer tracer layout for {ncols} columns and ndust={ndust}.")


def load_tracers(path: str | Path = "tracers.txt") -> dict[str, np.ndarray | int | str | list[str] | Layout | None]:
    ndust, rho_d, mol_w_sp, data = _load_header(path)
    layout = _infer_layout(data.shape[1], ndust)

    i = 0
    t_step = data[:, i]
    i += 1
    time = data[:, i]
    i += 1
    T = data[:, i]
    i += 1
    P_raw = data[:, i]
    i += 1
    grav = data[:, i]
    i += 1
    mu = data[:, i]
    i += 1

    VMR = data[:, i : i + layout.nbg]
    i += layout.nbg
    q_v = data[:, i : i + ndust]
    i += ndust
    q_0 = data[:, i]
    i += 1
    q_1 = data[:, i : i + ndust]
    i += ndust

    q_2 = None
    if layout.q2_mode == "scalar":
        q_2 = data[:, i]
        i += 1
    elif layout.q2_mode == "vector":
        q_2 = data[:, i : i + ndust]
        i += ndust

    v_f = data[:, i : i + layout.vf_count]
    i += layout.vf_count
    dTdt = data[:, i] if layout.has_dtdt and i < data.shape[1] else None

    # The original example scripts divide the pressure column by 1e5 to obtain bar,
    # implying that the output pressure is in Pa. Use that convention here, while
    # also providing cgs pressure for density calculations.
    p_bar = P_raw / 1.0e5
    p_cgs = p_bar * 1.0e6

    return {
        "path": str(path),
        "layout": layout,
        "ndust": ndust,
        "nbg": layout.nbg,
        "species": infer_species_names(mol_w_sp),
        "rho_d": rho_d,
        "mol_w_sp": mol_w_sp,
        "time_step": t_step,
        "time": time,
        "T": T,
        "P_raw": P_raw,
        "p_bar": p_bar,
        "p_cgs": p_cgs,
        "grav": grav,
        "mu": mu,
        "VMR": VMR,
        "q_v": q_v,
        "q_0": q_0,
        "q_1": q_1,
        "q_2": q_2,
        "v_f": v_f,
        "dTdt": dTdt,
    }


def vapour_mol_weights(species: list[str], mol_w_sp: np.ndarray) -> np.ndarray:
    return np.array([VAPOUR_MW.get(name, mw) for name, mw in zip(species, mol_w_sp)], dtype=float)


def compute_bulk(tr: dict) -> dict[str, np.ndarray | float]:
    ndust = int(tr["ndust"])
    T = np.asarray(tr["T"], dtype=float)
    p_bar = np.asarray(tr["p_bar"], dtype=float)
    p_cgs = np.asarray(tr["p_cgs"], dtype=float)
    mu = np.asarray(tr["mu"], dtype=float)
    rho_d = np.asarray(tr["rho_d"], dtype=float)
    q_0 = np.maximum(np.asarray(tr["q_0"], dtype=float), 1.0e-30)
    q_1 = np.maximum(np.asarray(tr["q_1"], dtype=float), 1.0e-30)

    nd_atm = p_cgs / (KB * T)
    rho = (p_cgs * mu * AMU) / (KB * T)
    N_c = q_0 * nd_atm
    rho_c = q_1 * rho[:, None]
    rho_c_t = np.maximum(rho_c.sum(axis=1), 1.0e-300)

    m_seed = V_SEED * rho_d[0]
    m_c = np.maximum(rho_c_t / np.maximum(N_c, 1.0e-300), m_seed)

    # Additive-volume mixed-material density: rho_mix = M_tot / V_tot.
    vol = rho_c / rho_d[None, :]
    vol_tot = np.maximum(vol.sum(axis=1), 1.0e-300)
    rho_d_m = rho_c_t / vol_tot

    r_c = np.maximum(((3.0 * m_c) / (4.0 * np.pi * rho_d_m)) ** (1.0 / 3.0), R_SEED) * 1.0e4
    V_mix = vol / vol_tot[:, None]
    m_mix = rho_c / rho_c_t[:, None]

    out = {
        "nd_atm": nd_atm,
        "rho": rho,
        "N_c": N_c,
        "rho_c": rho_c,
        "rho_c_t": rho_c_t,
        "m_c": m_c,
        "rho_d_m": rho_d_m,
        "r_c": r_c,
        "V_mix": V_mix,
        "m_mix": m_mix,
        "m_seed": m_seed,
    }

    q_2 = tr.get("q_2")
    if q_2 is not None:
        q_2_arr = np.maximum(np.asarray(q_2, dtype=float), 1.0e-300)
        if q_2_arr.ndim == 1:
            Z_c_t = q_2_arr * rho**2
        elif q_2_arr.shape[1] == ndust:
            # Backwards-compatible old layout: species q_2 values were summed.
            Z_c_t = q_2_arr.sum(axis=1) * rho**2
        else:
            raise ValueError(f"Unsupported q_2 shape: {q_2_arr.shape}")

        width_argument = np.maximum(N_c * Z_c_t / rho_c_t**2, 1.0)
        lnsig_raw = np.sqrt(np.maximum(np.log(width_argument), 0.0))
        sig_g = np.clip(np.exp(lnsig_raw), 1.01, 3.0)
        lnsig2 = np.log(sig_g) ** 2
        m_med = np.maximum(m_c * np.exp(-0.5 * lnsig2), m_seed)
        r_med = np.maximum(((3.0 * m_med) / (4.0 * np.pi * rho_d_m)) ** (1.0 / 3.0), R_SEED) * 1.0e4
        out.update({
            "Z_c_t": Z_c_t,
            "sig_g": sig_g,
            "lnsig2": lnsig2,
            "m_med": m_med,
            "r_med": r_med,
        })
    return out


def equilibrium_mass_fraction(
    species: list[str],
    T: np.ndarray,
    p_bar: np.ndarray,
    rho: np.ndarray,
    mol_w_v: np.ndarray,
    met: float = 0.0,
) -> np.ndarray:
    q_s = np.full((T.size, len(species)), np.nan)
    for j, name in enumerate(species):
        if name == "KCl":
            p_v = np.exp(
                -2.69250e4 / T + 3.39574e1 - 2.04903e-3 * T
                - 2.83957e-7 * T**2 + 1.82974e-10 * T**3
            )
        elif name == "ZnS":
            p_v = np.exp(
                -4.75507888e4 / T + 3.66993865e1 - 2.49490016e-3 * T
                + 7.29116854e-7 * T**2 - 1.12734453e-10 * T**3
            )
        elif name == "TiO2":
            p_v = np.exp(
                -7.70443e4 / T + 4.03144e1 - 2.59140e-3 * T
                + 6.02422e-7 * T**2 - 6.86899e-11 * T**3
            )
        elif name == "Al2O3":
            p_v = 1.0e6 * 10.0 ** (17.7 - 45892.6 / T - 1.66 * met)
        elif name == "Fe":
            p_v = 1.0e6 * 10.0 ** (7.23 - 20995.0 / T)
        elif name == "Mg2SiO4":
            p_v = 1.0e6 * 10.0 ** (14.88 - 32488.0 / T - 1.4 * met - 0.2 * np.log10(p_bar))
        else:
            continue
        rd_v = R_GAS / mol_w_v[j]
        q_s[:, j] = (p_v / (rd_v * T)) / rho
    return q_s


def setup_pressure_axis(ax, p_bar: np.ndarray | None = None) -> None:
    ax.set_yscale("log")
    ax.invert_yaxis()
    ax.set_ylabel(r"$p_{\rm gas}$ [bar]")
    if p_bar is not None:
        finite = np.asarray(p_bar)[np.isfinite(p_bar) & (np.asarray(p_bar) > 0.0)]
        if finite.size:
            ax.set_ylim(finite.max(), finite.min())


def finish_figure(fig, save_path: str | Path | None = None, show: bool = True) -> None:
    fig.tight_layout()
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    if show:
        import matplotlib.pyplot as plt
        plt.show()


def load_opacity(prefix: str = "opac") -> dict[str, np.ndarray]:
    wl = np.loadtxt(f"{prefix}_k.txt", max_rows=1)
    data_k = np.loadtxt(f"{prefix}_k.txt", skiprows=1)
    data_a = np.loadtxt(f"{prefix}_a.txt", skiprows=1)
    data_g = np.loadtxt(f"{prefix}_g.txt", skiprows=1)
    return {
        "wl": np.atleast_1d(wl),
        "pl": data_k[:, 2] / 1.0e6,
        "k_ext": data_k[:, 3:],
        "ssa": data_a[:, 3:],
        "g": data_g[:, 3:],
    }
