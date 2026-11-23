"""Read wcprop outputs and build per-event observables.

Light is stored in 10 ns bins with photon-weighted sums. Bins are grouped into
time clusters (sub-events): the prompt cluster, delayed decay electrons and
neutron-capture gammas. Each cluster gets a photon count, a mean time, a
centroid and a direction-isotropy proxy |sum n u| / sum n (1 for a single
straight track, small for isotropic light).
"""

import numpy as np
import pandas as pd


def read(prefix):
    light = pd.read_csv(f"{prefix}_light.csv")
    ncap = pd.read_csv(f"{prefix}_ncap.csv")
    ends = pd.read_csv(f"{prefix}_ends.csv")
    return light, ncap, ends


def clusters(light, gap_bins=5):
    """Group light bins of each event into clusters separated by > gap_bins empty bins."""
    df = light.sort_values(["event_id", "tbin"]).reset_index(drop=True)
    new = (df["event_id"].diff() != 0) | (df["tbin"].diff() > gap_bins)
    df["cluster"] = new.cumsum() - 1
    g = df.groupby("cluster")
    out = g[["n", "sx", "sy", "sz", "sux", "suy", "suz", "st"]].sum()
    out["event_id"] = g["event_id"].first()
    out["t"] = out["st"] / out["n"]
    for c in "xyz":
        out[c] = out[f"s{c}"] / out["n"]
    out["iso"] = np.sqrt(out["sux"]**2 + out["suy"]**2 + out["suz"]**2) / out["n"]
    out = out[["event_id", "t", "n", "x", "y", "z", "iso"]].reset_index(drop=True)
    out["rank"] = out.groupby("event_id").cumcount()
    return out


# Physical reference values for stopped negative muons in water (all bound to O
# after mu p -> mu O transfer): capture fraction and mean lifetime in 16O.
MU_MINUS_CAPTURE_O = 0.184
MU_MINUS_TAU_O_NS = 1795.0


def hydrogen_binding_weights(ends):
    """Per-event weights that undo the binding of stopped mu-/pi- to hydrogen.

    Geant4 11.0 binds ~30 % of stopped mu- and pi- in water to H and has no
    transfer of the mu p / pi p atom to oxygen. In real water these atoms
    transfer to O, so H binding (with ~0 % mu capture, or pi- p charge exchange)
    is unphysical. Events in which any primary-chain mu- or pi- stops on H get
    weight 0; the rest are scaled by 1 / (fraction of stops on O), separately
    for mu- and pi-, so totals are preserved.

    Returns a Series indexed by event_id (events without such stops have weight 1).
    """
    stops = ends[ends["pdg"].isin([13, -211]) & ends["target_z"].isin([1, 8])]
    w = pd.Series(1.0, index=pd.Index(ends["event_id"].unique(), name="event_id"))
    for pdg in (13, -211):
        s = stops[stops["pdg"] == pdg]
        if s.empty:
            continue
        f_o = (s["target_z"] == 8).mean()
        on_h = s.loc[s["target_z"] == 1, "event_id"].unique()
        on_o = s.loc[s["target_z"] == 8, "event_id"].unique()
        w.loc[on_o] *= 1.0 / f_o
        w.loc[on_h] = 0.0
    return w
