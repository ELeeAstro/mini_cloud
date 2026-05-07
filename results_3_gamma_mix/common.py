from pathlib import Path
import numpy as np

R_GAS = 8.31446261815324e7
KB = 1.380649e-16
AMU = 1.66053906660e-24
PA = 10.0
R_SEED = 1e-7
V_SEED = 4.0 / 3.0 * np.pi * R_SEED**3

SPECIES_BY_MW = {
    74.5513: "KCl",
    97.4450: "ZnS",
    79.8658: "TiO2",
    101.961: "Al2O3",
    55.8450: "Fe",
    140.693: "Mg2SiO4",
    18.01528: "H2O",
    17.03052: "NH3",
}

VAPOUR_MW = {
    "KCl": 74.5513,
    "ZnS": 97.4450,
    "TiO2": 79.8658,
    "Al2O3": 26.98153860,
    "Fe": 55.8450,
    "Mg2SiO4": 24.305,
    "H2O": 18.01528,
    "NH3": 17.03052,
}

SPECIES_TEX = {
    "KCl": r"${\rm KCl}$",
    "ZnS": r"${\rm ZnS}$",
    "TiO2": r"${\rm TiO_2}$",
    "Al2O3": r"${\rm Al_2O_3}$",
    "Fe": r"${\rm Fe}$",
    "Mg2SiO4": r"${\rm Mg_2SiO_4}$",
    "H2O": r"${\rm H_2O}$",
    "NH3": r"${\rm NH_3}$",
}


def _load_header(path):
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        first = np.fromstring(handle.readline(), sep=" ")
        if first.size != 1 or abs(first[0] - round(first[0])) > 1e-8:
            raise ValueError(
                f"{path} does not look like a multi-species tracer file with a 3-line header."
            )
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
    # Updated scalar-q2 layout:
    # t, time, T, P, grav, mu, VMR_bg(:), q_v(:), q_0, q_1(:), q_2, v_f(:), dTdt
    # with v_f length ndust+2.
    fixed_new_dtdt = 6 + ndust + 1 + ndust + 1 + (ndust + 2) + 1
    nbg = ncols - fixed_new_dtdt
    if nbg >= 1:
        return {"nbg": nbg, "layout": "scalar_q2", "has_dtdt": True}

    fixed_new_no_dtdt = 6 + ndust + 1 + ndust + 1 + (ndust + 2)
    nbg = ncols - fixed_new_no_dtdt
    if nbg >= 1:
        return {"nbg": nbg, "layout": "scalar_q2", "has_dtdt": False}

    # Legacy fallback: q_2(:) and v_f length 2*ndust+1.
    fixed_old_dtdt = 6 + ndust + 1 + ndust + ndust + (2 * ndust + 1) + 1
    nbg = ncols - fixed_old_dtdt
    if nbg >= 1:
        return {"nbg": nbg, "layout": "legacy_vector_q2", "has_dtdt": True}

    fixed_old_no_dtdt = 6 + ndust + 1 + ndust + ndust + (2 * ndust + 1)
    nbg = ncols - fixed_old_no_dtdt
    if nbg >= 1:
        return {"nbg": nbg, "layout": "legacy_vector_q2", "has_dtdt": False}

    raise ValueError(f"Could not infer tracer layout for {ncols} columns and ndust={ndust}.")


def load_tracers(path="tracers.txt"):
    ndust, rho_d, mol_w_sp, data = _load_header(path)
    layout = _infer_layout(data.shape[1], ndust)
    nbg = layout["nbg"]

    i0 = 0
    time_step = data[:, i0]
    i0 += 1
    time = data[:, i0]
    i0 += 1

    T = data[:, i0]
    P_pa = data[:, i0 + 1]
    grav = data[:, i0 + 2]
    mu = data[:, i0 + 3]
    i0 += 4

    VMR = data[:, i0 : i0 + nbg]
    i0 += nbg

    q_v = data[:, i0 : i0 + ndust]
    i0 += ndust

    q_0 = data[:, i0]
    i0 += 1

    q_1 = data[:, i0 : i0 + ndust]
    i0 += ndust

    if layout["layout"] == "scalar_q2":
        q_2 = data[:, i0]
        i0 += 1
        v_f = data[:, i0 : i0 + ndust + 2]
        i0 += ndust + 2
    else:
        q_2 = np.sum(data[:, i0 : i0 + ndust], axis=1)
        i0 += ndust
        v_f = data[:, i0 : i0 + 2 * ndust + 1]
        i0 += 2 * ndust + 1

    dTdt = data[:, i0] if layout["has_dtdt"] else None

    species = infer_species_names(mol_w_sp)
    return {
        "ndust": ndust,
        "nbg": nbg,
        "layout": layout["layout"],
        "species": species,
        "species_tex": [SPECIES_TEX.get(s, s) for s in species],
        "rho_d": rho_d,
        "mol_w_sp": mol_w_sp,
        "time_step": time_step,
        "time": time,
        "T": T,
        "P_pa": P_pa,
        "P_bar": P_pa / 1.0e5,
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


def compute_bulk(tr, nu_min=1.0e-2, nu_max=1.0e2):
    T = tr["T"]
    P_pa = tr["P_pa"]
    mu = tr["mu"]
    rho_d = np.asarray(tr["rho_d"])

    q_0 = np.maximum(tr["q_0"], 1.0e-300)
    q_1 = np.maximum(tr["q_1"], 1.0e-300)
    q_2 = np.maximum(tr["q_2"], 1.0e-300)

    p_cgs = P_pa * PA
    nd_atm = p_cgs / (KB * T)
    rho = p_cgs * mu * AMU / (KB * T)

    N_c = np.maximum(q_0 * nd_atm, 1.0e-300)
    rho_c = q_1 * rho[:, None]
    rho_c_t = np.sum(rho_c, axis=1)

    m_seed = V_SEED * rho_d[0]
    m_c = np.maximum(rho_c_t / N_c, m_seed)

    vol = rho_c / rho_d[None, :]
    vol_tot = np.sum(vol, axis=1)

    rho_d_m = np.full_like(N_c, rho_d[0])
    good = (rho_c_t > 0.0) & (vol_tot > 0.0)
    rho_d_m[good] = rho_c_t[good] / vol_tot[good]

    r_c = np.maximum(
        ((3.0 * m_c) / (4.0 * np.pi * rho_d_m)) ** (1.0 / 3.0),
        R_SEED,
    ) * 1.0e4

    V_mix = np.zeros_like(rho_c)
    V_mix[good, :] = vol[good, :] / vol_tot[good, None]

    m_mix = np.zeros_like(rho_c)
    good_mass = rho_c_t > 0.0
    m_mix[good_mass, :] = rho_c[good_mass, :] / rho_c_t[good_mass, None]

    Z_c_t = q_2 * rho**2
    sig2 = np.maximum(Z_c_t / N_c - m_c**2, m_seed**2)
    nu = np.clip(m_c**2 / sig2, nu_min, nu_max)
    lam = m_c / nu

    return {
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
        "Z_c_t": Z_c_t,
        "sig2": sig2,
        "nu": nu,
        "lam": lam,
    }


def vapour_mol_weights(species, mol_w_sp):
    return np.array([VAPOUR_MW.get(name, mw) for name, mw in zip(species, mol_w_sp)], dtype=float)


def equilibrium_mass_fraction(species, T, p_bar, rho, mol_w_v, met=0.0):
    species = list(species)
    q_s = np.full((T.size, len(species)), np.nan)
    for j, name in enumerate(species):
        if name == "KCl":
            p_v = np.exp(
                -2.69250e4 / T
                + 3.39574e1
                - 2.04903e-3 * T
                - 2.83957e-7 * T**2
                + 1.82974e-10 * T**3
            )
        elif name == "ZnS":
            p_v = np.exp(
                -4.75507888e4 / T
                + 3.66993865e1
                - 2.49490016e-3 * T
                + 7.29116854e-7 * T**2
                - 1.12734453e-10 * T**3
            )
        elif name == "TiO2":
            p_v = np.exp(
                -7.70443e4 / T
                + 4.03144e1
                - 2.59140e-3 * T
                + 6.02422e-7 * T**2
                - 6.86899e-11 * T**3
            )
        elif name == "Al2O3":
            p_v = 1.0e6 * 10.0 ** (17.7 - 45892.6 / T - 1.66 * met)
        elif name == "Fe":
            p_v = 1.0e6 * 10.0 ** (7.23 - 20995.0 / T)
        elif name == "Mg2SiO4":
            p_v = 1.0e6 * 10.0 ** (14.88 - 32488.0 / T - 1.4 * met - 0.2 * np.log10(p_bar))
        elif name == "H2O":
            p_v = np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T) * PA
        elif name == "NH3":
            p_v = np.exp(1.596e1 - 3.537e3 / T - 3.310e4 / T**2 + 1.742e6 / T**3 - 2.995e7 / T**4) * 1.0e6
        else:
            continue
        rd_v = R_GAS / mol_w_v[j]
        q_s[:, j] = p_v / (rd_v * T * rho)
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
