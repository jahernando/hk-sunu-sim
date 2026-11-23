"""Build DSNB-like candidates from wcprop output and label them by truth.

Event building (simplified SK low-energy logic):

- a *trigger* is a light cluster with E_rec >= ``e_trigger``;
- the *candidate* is the first trigger of the event, which is the decay
  electron when the primary muon/pion is invisible;
- *pre-activity*: light before the candidate within ``pre_window_ns`` with
  E_rec >= ``e_pre`` (an invisible muon can still leave some light);
- *decay electron after*: a later cluster within ``post_window_ns`` with
  E_rec >= ``e_post``.

Truth category of the candidate (SK-IV names):

- ``decay_e``: the candidate is the end product of a stopped mu/pi (matched
  in time to the ends table);
- ``nue_cc``, ``ncqe``, ``mupi``, ``nc_other``, ``other``: prompt candidates,
  from the generator mode and the visible primaries.
"""

import numpy as np
import pandas as pd

from . import response

CHERENKOV_P_THRESHOLD = {13: 120.4, 211: 159.2}   # MeV/c in water (n = 1.333)


def candidates(clusters, ends, det, rng, e_trigger=8.0, e_pre=6.0, e_post=6.0,
               pre_window_ns=50e3, post_window_ns=50e3, match_ns=50.0):
    cl = clusters.copy()
    cl["e_rec"] = response.energy(cl["n"], det, rng)
    trig = cl[cl["e_rec"] >= e_trigger]
    cand = trig.sort_values(["event_id", "t"]).groupby("event_id").head(1).set_index("event_id")

    merged = cl.merge(cand[["t"]].rename(columns={"t": "t0"}), left_on="event_id", right_index=True)
    dt = merged["t"] - merged["t0"]
    pre = merged[(dt < -match_ns) & (dt > -pre_window_ns) & (merged["e_rec"] >= e_pre)]
    post = merged[(dt > match_ns) & (dt < post_window_ns) & (merged["e_rec"] >= e_post)]
    cand["pre_activity"] = cand.index.isin(pre["event_id"])
    cand["decay_e_after"] = cand.index.isin(post["event_id"])

    # is the candidate the end product (decay electron) of a stopped mu/pi?
    stops = ends[ends["fate"] == "decay"][["event_id", "t"]].rename(columns={"t": "t_end"})
    m = stops.merge(cand[["t"]], left_on="event_id", right_index=True)
    is_decay = m[np.abs(m["t"] - m["t_end"]) < match_ns]["event_id"].unique()
    cand["from_decay"] = cand.index.isin(is_decay)
    return cand


def visible_mupi(primaries):
    """Events with a primary mu or charged pi above the Cherenkov threshold."""
    p = np.sqrt(primaries["px"]**2 + primaries["py"]**2 + primaries["pz"]**2)
    thr = primaries["pdg"].abs().map(CHERENKOV_P_THRESHOLD)
    return primaries.loc[p > thr, "event_id"].unique()


def category(cand, gen, primaries):
    """Truth category per candidate.

    gen: generator table indexed by event_id with columns nu_pdg, cc (bool),
    mode (string; 'qel' for quasi-elastic).
    """
    g = gen.reindex(cand.index)
    vis = cand.index.isin(visible_mupi(primaries))
    cat = np.where(cand["from_decay"], "decay_e",
          np.where(vis, "mupi",
          np.where(g["cc"] & g["nu_pdg"].abs().eq(12), "nue_cc",
          np.where(~g["cc"] & g["mode"].eq("qel"), "ncqe",
          np.where(~g["cc"], "nc_other", "other")))))
    return pd.Series(cat, index=cand.index, name="category")
