import numpy as np
import pytest

from hk_sunu_sim import genie

O16 = 1000080160


def synthetic_gst():
    """Five events: nu_mu CCQE on O, nubar_e CCQE on H, NC elastic on O,
    nu_mu CC res with one pion, nu_e NC DIS. GeV units, as in gst."""
    a = np.array
    pdgf = [a([2212]), a([2112]), a([2112, 22]), a([2212, 211]), a([211, -211, 111, 2112])]
    pdgi = [a([2212]), a([2112]), a([2112]), a([2212, 211]), a([211, -211, 111])]
    ef = [np.full(len(p), 1.0) for p in pdgf]
    pf = [np.full(len(p), 0.1) for p in pdgf]
    z = np.zeros(5)
    return {
        "iev": a([0, 1, 2, 3, 4]),
        "neu": a([14, -12, 14, 14, 12]),
        "fspl": a([13, -11, 14, 13, 12]),
        "tgt": a([O16, 2212, O16, O16, O16]),        # gst: free proton -> 2212
        "Z": a([8, 1, 8, 8, 8]), "A": a([16, 1, 16, 16, 16]),
        "hitnuc": a([2112, 2212, 2112, 2212, 2212]),
        "cc": a([1, 1, 0, 1, 0]), "nc": a([0, 0, 1, 0, 1]),
        "qel": a([1, 1, 1, 0, 0]), "mec": z, "res": a([0, 0, 0, 1, 0]),
        "dis": a([0, 0, 0, 0, 1]), "coh": z, "dfr": z, "nuel": z, "imd": z, "imdanh": z, "em": z,
        "Ev": a([0.5, 0.03, 0.4, 1.5, 3.0]),
        "pxv": z, "pyv": z, "pzv": a([0.5, 0.03, 0.4, 1.5, 3.0]),
        "El": a([0.4, 0.02, 0.35, 1.0, 2.0]),
        "pxl": a([0.1, 0.0, 0.0, 0.0, 0.0]), "pyl": z, "pzl": a([0.37, 0.02, 0.35, 0.99, 2.0]),
        "Q2": a([0.2, 0.0, 0.05, 0.5, 1.0]),
        "xsec": a([1.0, 0.1, 0.5, 2.0, 3.0]),
        "pdgf": pdgf, "Ef": ef, "pxf": pf, "pyf": pf, "pzf": pf, "pdgi": pdgi,
    }


def test_convert_channels_and_targets():
    gen, parts, prim = genie.convert(synthetic_gst(), batch=7)
    assert list(gen["channel"]) == ["ccqe", "ccqe", "ncel", "cc1pi", "ncother"]
    assert list(gen["target_z"]) == [8, 1, 8, 8, 8]
    assert list(gen["n_pi"]) == [0, 0, 0, 1, 3]
    assert list(gen["n_neutron"]) == [0, 1, 1, 0, 1]
    assert np.allclose(gen["e_nu"], [500, 30, 400, 1500, 3000])
    assert (gen["batch"] == 7).all()
    assert set(gen["channel"]) <= set(genie.CHANNELS)


def test_particles_and_primaries():
    gen, parts, prim = genie.convert(synthetic_gst())
    assert len(parts) == 5 + 10                       # 5 leptons + 10 hadrons
    assert parts.groupby("event_id")["is_lepton"].first().all()   # lepton first
    assert not prim["pdg"].abs().isin([12, 14, 16]).any()          # NC neutrinos dropped
    assert len(prim) == len(parts) - 2
    assert prim["px"].max() == pytest.approx(100.0)                # GeV -> MeV


def test_events_category_contract():
    """events.category needs nu_pdg, cc and mode == 'qel' for NCQE."""
    gen, _, _ = genie.convert(synthetic_gst())
    g = gen.set_index("event_id")
    assert g.loc[2, "mode"] == "qel" and not g.loc[2, "cc"]


def test_write_roundtrip(tmp_path):
    gen, parts, prim = genie.convert(synthetic_gst())
    meta = genie.metadata("3.04.00", "G18_10a_02_11b", "test", "hkkm2014", 5, 1,
                          {"O16": 0.888, "H1": 0.112}, "synthetic")
    genie.write(gen, parts, prim, tmp_path / "t", meta)
    assert genie.read_metadata(tmp_path / "t_atm_gen.parquet")["tune"] == "G18_10a_02_11b"


def test_read_gst_with_uproot(tmp_path):
    uproot = pytest.importorskip("uproot")
    ak = pytest.importorskip("awkward")
    gst = synthetic_gst()
    branches = {k: (ak.Array([list(x) for x in v]) if isinstance(v, list) else v)
                for k, v in gst.items()}
    with uproot.recreate(tmp_path / "t.gst.root") as f:
        f.mktree("gst", {k: v.type if isinstance(v, ak.Array) else v.dtype
                         for k, v in branches.items()})
        f["gst"].extend(branches)
    back = genie.read_gst(tmp_path / "t.gst.root")
    gen, _, _ = genie.convert(back)
    assert list(gen["channel"]) == ["ccqe", "ccqe", "ncel", "cc1pi", "ncother"]
