"""Convert GENIE ``gst`` summary trees into the generator-independent tables.

GENIE writes ``ghep`` event records; ``gntpc -f gst`` flattens them into the
``gst`` tree (one entry per interaction, GeV units). This module reads that
tree with ``uproot`` and writes:

- ``atm_gen``: one row per interaction (neutrino, target, process, common
  channel, lepton kinematics, pion/neutron counts);
- ``atm_particles``: the primary lepton and the hadrons leaving the nucleus
  (after FSI);
- the ``wcprop`` input (``event_id,pdg,px,py,pz`` in MeV/c), without neutrinos.

Units on output: MeV, MeV/c, MeV^2, cm^2. The common channel table
(:data:`CHANNELS`) is shared with the NEUT converter so that the downstream
chain does not depend on the generator.

Notes on GENIE content:

- the target is taken from the ``Z``/``A`` branches: gst writes ``tgt = 2212``
  (not 1000010010) for the free proton;
- the nuclear remnant is not in the final-state list, but the 16O
  de-excitation gammas are (seen in G18_10a_02_11b, ~45 % of NC elastic events
  in the smoke test); coherent events carry the recoiling 16O ion;
- ``xsec`` is not written by gntpc 3.04 (``w_xsec`` = 0); the normalisation
  comes from the splines;
- GENIE pseudo-particles (pdg >= 2e9) are dropped from the ``wcprop`` input.
"""

import json

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .schema import ATM_GEN, ATM_PARTICLES

GEV = 1000.0

# Common interaction channels (generator independent). NEUT modes map to the
# same names in neut.py.
CHANNELS = {
    "ccqe": "CC quasi-elastic (incl. CC on free proton for nubar)",
    "cc2p2h": "CC 2p2h / MEC",
    "cc1pi": "CC resonant, one pion before FSI",
    "cccoh": "CC coherent pion",
    "ccmpi_dis": "CC multi-pion and DIS (incl. resonant with != 1 pion)",
    "ncel": "NC elastic (NCQE on 16O)",
    "nc1pi": "NC resonant, one pion before FSI",
    "nccoh": "NC coherent pion",
    "ncother": "NC multi-pion, DIS and NC 2p2h",
    "nue_elastic": "neutrino-electron elastic / inverse muon decay",
    "other": "anything else (diffractive, charm-only, ...)",
}

PROCESS_FLAGS = ("qel", "mec", "res", "dis", "coh", "dfr", "nuel", "imd", "imdanh", "em")

NEUTRINOS = {12, 14, 16}
PIONS = {111, 211}

# Branches read from gst (missing optional ones are tolerated)
SCALARS = ("iev", "neu", "fspl", "tgt", "Z", "A", "hitnuc", "cc", "nc", "Ev",
           "pxv", "pyv", "pzv", "El", "pxl", "pyl", "pzl", "Q2", "xsec") + PROCESS_FLAGS
JAGGED = ("pdgf", "Ef", "pxf", "pyf", "pzf", "pdgi")


# --- reading ----------------------------------------------------------------

def read_gst(path, tree="gst", entry_stop=None):
    """Read the branches of a gst file into a dict of numpy arrays / lists of arrays."""
    import uproot
    with uproot.open(path) as f:
        t = f[tree]
        names = [b for b in SCALARS + JAGGED if b in t.keys()]
        return {b: t[b].array(library="np", entry_stop=entry_stop) for b in names}


def _jagged(x):
    """(counts, flat values) from a list/object-array of 1D arrays."""
    counts = np.array([len(v) for v in x], dtype=np.int64)
    flat = np.concatenate([np.asarray(v) for v in x]) if len(x) else np.array([])
    return counts, flat


# --- classification ---------------------------------------------------------

def process(gst):
    """Process name per event from the gst boolean flags."""
    n = len(gst["Ev"])
    out = np.full(n, "other", dtype=object)
    for name in ("dfr", "em", "imdanh", "imd", "nuel", "coh", "dis", "res", "mec", "qel"):
        if name in gst:
            flag = np.asarray(gst[name]).astype(bool)
            out[flag] = "imd" if name == "imdanh" else name
    return out


def count_pdg(x, codes):
    """Number of entries with |pdg| in codes, per event."""
    counts, flat = _jagged(x)
    hit = np.isin(np.abs(flat.astype(np.int64)), list(codes)).astype(np.int64)
    idx = np.repeat(np.arange(len(counts)), counts)
    return np.bincount(idx, weights=hit, minlength=len(counts)).astype(np.int16)


def channel(cc, mode, n_pi_primary):
    """Common channel from CC/NC, process and pion multiplicity before FSI."""
    cc = np.asarray(cc, bool)
    mode = np.asarray(mode, dtype=object)
    one_pi = np.asarray(n_pi_primary) == 1
    out = np.full(len(cc), "other", dtype=object)
    for c, pre in ((cc, "cc"), (~cc, "nc")):
        out[c & (mode == "dis")] = "ccmpi_dis" if pre == "cc" else "ncother"
        out[c & (mode == "res") & ~one_pi] = "ccmpi_dis" if pre == "cc" else "ncother"
        out[c & (mode == "res") & one_pi] = f"{pre}1pi"
        out[c & (mode == "coh")] = f"{pre}coh"
    out[cc & (mode == "qel")] = "ccqe"
    out[~cc & (mode == "qel")] = "ncel"
    out[cc & (mode == "mec")] = "cc2p2h"
    out[~cc & (mode == "mec")] = "ncother"
    out[np.isin(mode, ["nuel", "imd"])] = "nue_elastic"
    return out


# --- conversion --------------------------------------------------------------

def convert(gst, batch=0):
    """gst dict -> (atm_gen, atm_particles, wcprop primaries) DataFrames."""
    n = len(gst["Ev"])
    event_id = np.asarray(gst["iev"], np.int64) if "iev" in gst else np.arange(n)
    cc = np.asarray(gst["cc"]).astype(bool)
    mode = process(gst)
    n_pi_pre = count_pdg(gst["pdgi"], PIONS) if "pdgi" in gst else np.full(n, -1, np.int16)
    tgt = np.asarray(gst["tgt"], np.int64)

    pv = np.column_stack([gst["pxv"], gst["pyv"], gst["pzv"]]).astype(float)
    u = pv / np.linalg.norm(pv, axis=1)[:, None]
    pl = np.column_stack([gst["pxl"], gst["pyl"], gst["pzl"]]).astype(float)
    pl_norm = np.linalg.norm(pl, axis=1)
    cos_l = np.divide((pl * u).sum(1), pl_norm, out=np.ones(n), where=pl_norm > 0)

    gen = pd.DataFrame({
        "batch": np.int32(batch),
        "event_id": event_id,
        "nu_pdg": np.asarray(gst["neu"], np.int32),
        "e_nu": np.asarray(gst["Ev"], float) * GEV,
        "ux": u[:, 0], "uy": u[:, 1], "uz": u[:, 2],
        "target_z": np.asarray(gst["Z"], np.int16) if "Z" in gst else ((tgt // 10000) % 1000).astype(np.int16),
        "target_a": np.asarray(gst["A"], np.int16) if "A" in gst else ((tgt // 10) % 1000).astype(np.int16),
        "hit_nucleon": np.asarray(gst.get("hitnuc", np.zeros(n)), np.int32),
        "cc": cc,
        "mode": mode,
        "channel": channel(cc, mode, n_pi_pre),
        "gen_code": [f"{'CC' if c else 'NC'}:{m}" for c, m in zip(cc, mode)],
        "n_pi_primary": n_pi_pre,
        "n_pi": count_pdg(gst["pdgf"], PIONS),
        "n_neutron": count_pdg(gst["pdgf"], {2112}),
        "lepton_pdg": np.asarray(gst["fspl"], np.int32),
        "e_lepton": np.asarray(gst["El"], float) * GEV,
        "cos_lepton": cos_l,
        "q2": np.asarray(gst.get("Q2", np.zeros(n)), float) * GEV**2,
        "w_xsec": np.asarray(gst.get("xsec", np.zeros(n)), float) * 1e-38,
    })

    counts, pdg = _jagged(gst["pdgf"])
    had = pd.DataFrame({
        "batch": np.int32(batch),
        "event_id": np.repeat(event_id, counts),
        "pdg": pdg.astype(np.int32),
        "e": _jagged(gst["Ef"])[1] * GEV,
        "px": _jagged(gst["pxf"])[1] * GEV,
        "py": _jagged(gst["pyf"])[1] * GEV,
        "pz": _jagged(gst["pzf"])[1] * GEV,
        "is_lepton": False,
    })
    lep = pd.DataFrame({
        "batch": np.int32(batch), "event_id": event_id, "pdg": gen["lepton_pdg"].to_numpy(),
        "e": gen["e_lepton"].to_numpy(), "px": pl[:, 0] * GEV, "py": pl[:, 1] * GEV,
        "pz": pl[:, 2] * GEV, "is_lepton": True,
    })
    parts = (pd.concat([lep, had], ignore_index=True)
             .sort_values(["event_id", "is_lepton"], ascending=[True, False], kind="stable")
             .reset_index(drop=True))

    prim = parts[~parts["pdg"].abs().isin(NEUTRINOS) & (parts["pdg"].abs() < 2_000_000_000)]
    prim = prim[["event_id", "pdg", "px", "py", "pz"]].reset_index(drop=True)
    return gen, parts, prim


# --- writing -------------------------------------------------------------------

def metadata(version, tune, splines, flux, n_events, seed, targets, source, **extra):
    """Generator metadata stored in the Parquet schema (key ``hk_sunu_sim``)."""
    return {"generator": "genie", "version": version, "tune": tune, "splines": splines,
            "flux": flux, "n_events": int(n_events), "seed": seed, "targets": targets,
            "source": source, **extra}


def _write(df, schema, path, meta):
    table = pa.Table.from_pandas(df[schema.names], schema=schema, preserve_index=False)
    table = table.replace_schema_metadata({**(table.schema.metadata or {}),
                                           b"hk_sunu_sim": json.dumps(meta).encode()})
    pq.write_table(table, path)


def write(gen, parts, prim, prefix, meta):
    """Write <prefix>_atm_gen.parquet, <prefix>_atm_particles.parquet, <prefix>_primaries.csv."""
    _write(gen, ATM_GEN, f"{prefix}_atm_gen.parquet", meta)
    _write(parts, ATM_PARTICLES, f"{prefix}_atm_particles.parquet", meta)
    prim.to_csv(f"{prefix}_primaries.csv", index=False, float_format="%.6g")


def read_metadata(path):
    return json.loads(pq.read_schema(path).metadata[b"hk_sunu_sim"])
