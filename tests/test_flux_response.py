import numpy as np
import pandas as pd
import pytest

from hk_sunu_sim import events, flux, response


@pytest.mark.skipif(not (flux.DEFAULT_DIR / "kam-ally-20-01-solmin.d.gz").exists(),
                    reason="HKKM tables not downloaded")
def test_hkkm_integrals_and_ratio():
    t = flux.load("min")
    e = np.geomspace(100, 1e4, 500)
    numu = np.trapezoid(t.flux_solid_angle_integrated(14, e), e)
    assert 1.2 < numu < 1.8           # cm^-2 s^-1 above 0.1 GeV
    assert t.flux(14, 50.0, 0.5)[0] == 0.0
    assert t.flux(14, 50.0, 0.5, extrapolate=True)[0] > t.flux(14, 100.0, 0.5)[0]


def test_energy_scale_is_unbiased():
    rng = np.random.default_rng(0)
    det = response.SK4
    e = response.energy(np.full(20000, 20 * det.photons_per_mev), det, rng)
    assert abs(e.mean() - 20) < 0.05
    assert abs(e.std() - np.sqrt(20 / det.pe_per_mev)) < 0.05


def test_neutron_tags_window_and_efficiency():
    rng = np.random.default_rng(1)
    det = response.Detector("t", 6, 50, 10, (18.0, 535.0), eps_tag=1.0, mistag=0.0)
    t0 = pd.Series([0.0, 0.0], index=pd.Index([0, 1], name="event_id"))
    ncap = pd.DataFrame({"event_id": [0, 0, 1], "t": [100e3, 10e3, 600e3]})
    n = response.neutron_tags(t0, ncap, det, rng)
    assert n.tolist() == [1, 0]


def test_category_logic():
    cand = pd.DataFrame({"from_decay": [True, False, False, False]}, index=pd.Index([0, 1, 2, 3], name="event_id"))
    gen = pd.DataFrame({"nu_pdg": [14, 12, 14, 14], "cc": [True, True, False, True],
                        "mode": ["qel", "qel", "qel", "res"]}, index=cand.index)
    prim = pd.DataFrame({"event_id": [3], "pdg": [13], "px": [0.0], "py": [0.0], "pz": [300.0]})
    assert events.category(cand, gen, prim).tolist() == ["decay_e", "nue_cc", "ncqe", "mupi"]
