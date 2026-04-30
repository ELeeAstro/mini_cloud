from pathlib import Path

import numpy as np


R_GAS = 8.31446261815324e7
KB = 1.380649e-16
AMU = 1.66053906660e-24
R_SEED = 1e-7
V_SEED = 4.0 / 3.0 * np.pi * R_SEED**3


SPECIES_BY_MW = {
    74.5513: "KCl",
    97.4450: "ZnS",
    79.8658: "TiO2",
    101.961: "Al2O3",
    55.8450: "Fe",
    140.693: "Mg2SiO4",
}

VAPOUR_MW = {
    "KCl": 74.5513,
    "ZnS": 97.4450,
    "TiO2": 79.8658,
    "Al2O3": 26.98153860,
    "Fe": 55.8450,
    "Mg2SiO4": 24.305,
}


def _load_header(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        first = np.fromstring(handle.readline(), sep=" ")
    if first.size != 1 or abs(first[0] - round(first[0])) > 1e-8:
        raise ValueError(f"{path} does not look like a multi-species tracer file with a 3-line header.")

    ndust = int(round(first[0]))
    rho_d = np.atleast_1d(np.loadtxt(path, skiprows=1, max_rows=1))
    mol_w_sp = np.atleast_1d(np.loadtxt(path, skiprows=2, max_rows=1))
    data = np.loadtxt(path, skiprows=3)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return ndust, rho_d, mol_w_sp, data


def infer_species_names(mol_w_sp):
    species = []
    for idx, mol in enumerate(np.atleast_1d(mol_w_sp)):
        name = None
        for ref, candidate in SPECIES_BY_MW.items():
            if np.isclose(mol, ref, rtol=0.0, atol=5e-4):
                name = candidate
                break
        species.append(name or f"sp{idx + 1}")
    return species


def _infer_layout(ncols, ndust):
    candidates = [
        (True, True),
        (False, True),
        (True, False),
        (False, False),
    ]
    for has_q2, has_dtdt in candidates:
        fixed = 6 + ndust + 1 + ndust + 1 + (ndust if has_q2 else 0) + (1 if has_dtdt else 0)
        nbg = ncols - fixed
        if nbg >= 1:
            return nbg, has_q2, has_dtdt
    raise ValueError(f"Could not infer tracer layout for {ncols} columns and ndust={ndust}.")


def load_tracers(path="tracers.txt"):
    ndust, rho_d, mol_w_sp, data = _load_header(path)
    nbg, has_q2, has_dtdt = _infer_layout(data.shape[1], ndust)

    i0 = 0
    time_index = slice(i0, i0 + 2)
    i0 += 2
    thermo_index = slice(i0, i0 + 4)
    i0 += 4
    bg_index = slice(i0, i0 + nbg)
    i0 += nbg
    qv_index = slice(i0, i0 + ndust)
    i0 += ndust
    q0_index = i0
    i0 += 1
    q1_index = slice(i0, i0 + ndust)
    i0 += ndust
    q2_index = slice(i0, i0 + ndust) if has_q2 else None
    if has_q2:
        i0 += ndust
    vf_index = i0
    i0 += 1
    dtdt_index = i0 if has_dtdt else None

    species = infer_species_names(mol_w_sp)

    return {
        "ndust": ndust,
        "nbg": nbg,
        "species": species,
        "rho_d": rho_d,
        "mol_w_sp": mol_w_sp,
        "time_step": data[:, time_index.start],
        "time": data[:, time_index.stop - 1],
        "T": data[:, thermo_index.start],
        "P": data[:, thermo_index.start + 1],
        "grav": data[:, thermo_index.start + 2],
        "mu": data[:, thermo_index.start + 3],
        "VMR": data[:, bg_index],
        "q_v": data[:, qv_index],
        "q_0": data[:, q0_index],
        "q_1": data[:, q1_index],
        "q_2": data[:, q2_index] if has_q2 else None,
        "v_f": data[:, vf_index],
        "dTdt": data[:, dtdt_index] if has_dtdt else None,
    }


def compute_bulk(tr):
    T = tr["T"]
    P = tr["P"]
    mu = tr["mu"]
    rho_d = np.asarray(tr["rho_d"])
    q_0 = np.maximum(tr["q_0"], 1e-30)
    q_1 = np.maximum(tr["q_1"], 1e-30)

    nd_atm = P / (KB * T)
    rho = (P * mu * AMU) / (KB * T)
    N_c = q_0 * nd_atm
    rho_c = q_1 * rho[:, None]
    rho_c_t = np.maximum(rho_c.sum(axis=1), 1e-30)

    m_seed = V_SEED * rho_d[0]
    m_c = np.maximum(rho_c_t / np.maximum(N_c, 1e-30), m_seed)
    rho_d_m = (rho_c / rho_c_t[:, None]) @ rho_d
    r_c = np.maximum(((3.0 * m_c) / (4.0 * np.pi * rho_d_m)) ** (1.0 / 3.0), R_SEED) * 1e4

    vol = rho_c / rho_d[None, :]
    vol_tot = np.maximum(vol.sum(axis=1), 1e-30)
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

    if tr["q_2"] is not None:
        q_2 = np.maximum(tr["q_2"], 1e-30)
        Z_c_t = q_2.sum(axis=1) * rho**2
        sig2 = np.maximum(Z_c_t / np.maximum(N_c, 1e-30) - m_c**2, m_seed**2)
        nu = np.clip(m_c**2 / sig2, 1.0e-2, 1.0e2)
        lam = m_c / nu
        out["Z_c_t"] = Z_c_t
        out["sig2"] = sig2
        out["nu"] = nu
        out["lam"] = lam

    return out


def vapour_mol_weights(species, mol_w_sp):
    return np.array([VAPOUR_MW.get(name, mw) for name, mw in zip(species, mol_w_sp)], dtype=float)


def equilibrium_mass_fraction(species, T, p_bar, rho, mol_w_v, met=0.0):
    species = list(species)
    q_s = np.full((T.size, len(species)), np.nan)

    for j, name in enumerate(species):
        if name == "KCl":
            p_v = np.exp(
                -2.69250e4 / T + 3.39574e1 - 2.04903e-3 * T - 2.83957e-7 * T**2 + 1.82974e-10 * T**3
            )
        elif name == "ZnS":
            p_v = np.exp(
                -4.75507888e4 / T + 3.66993865e1 - 2.49490016e-3 * T + 7.29116854e-7 * T**2
                - 1.12734453e-10 * T**3
            )
        elif name == "TiO2":
            p_v = np.exp(
                -7.70443e4 / T + 4.03144e1 - 2.59140e-3 * T + 6.02422e-7 * T**2 - 6.86899e-11 * T**3
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


def load_opacity(prefix="opac"):
    wl = np.loadtxt(f"{prefix}_k.txt", max_rows=1)
    data_k = np.loadtxt(f"{prefix}_k.txt", skiprows=1)
    data_a = np.loadtxt(f"{prefix}_a.txt", skiprows=1)
    data_g = np.loadtxt(f"{prefix}_g.txt", skiprows=1)
    return {
        "wl": wl,
        "pl": data_k[:, 2] / 1e6,
        "k_ext": data_k[:, 3:],
        "ssa": data_a[:, 3:],
        "g": data_g[:, 3:],
    }
