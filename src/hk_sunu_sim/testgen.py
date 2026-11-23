"""Test primaries for the wcprop propagation (no neutrino generator needed).

Each sample is a DataFrame with columns event_id, pdg, px, py, pz (MeV/c),
the input format of ``wcprop``. Used to calibrate the visible-energy proxy and
to check the low-energy physics (Michel spectrum, muon capture, neutron
capture) before running generator events.
"""

import numpy as np
import pandas as pd

MASS = {11: 0.51099895, 13: 105.6583755, 211: 139.57039, 2112: 939.56542, 22: 0.0}


def isotropic(n, rng):
    cost = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    sint = np.sqrt(1 - cost**2)
    return np.column_stack([sint * np.cos(phi), sint * np.sin(phi), cost])


def momentum(pdg, ekin):
    m = MASS[abs(pdg)]
    return np.sqrt(ekin * (ekin + 2 * m))


def single(pdg, ekin, n, rng, first_id=0):
    """n events with one particle of kinetic energy ekin (scalar or array)."""
    ekin = np.broadcast_to(np.asarray(ekin, float), (n,))
    u = isotropic(n, rng) * momentum(pdg, ekin)[:, None]
    return pd.DataFrame({"event_id": np.arange(first_id, first_id + n), "pdg": pdg,
                         "px": u[:, 0], "py": u[:, 1], "pz": u[:, 2]})


def ncqe_like(n, rng, e_gamma=6.18, e_neutron=10.0, first_id=0):
    """A de-excitation gamma plus a knocked-out neutron, both isotropic."""
    g = single(22, e_gamma, n, rng, first_id)
    nn = single(2112, e_neutron, n, rng, first_id)
    return pd.concat([g, nn]).sort_values("event_id", kind="stable").reset_index(drop=True)


def write(df, path):
    df.to_csv(path, index=False, float_format="%.6g")
