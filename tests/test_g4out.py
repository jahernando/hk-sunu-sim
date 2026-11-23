import pandas as pd

from hk_sunu_sim.g4out import hydrogen_binding_weights


def test_hydrogen_binding_weights_preserve_total():
    ends = pd.DataFrame({
        "event_id": [0, 1, 2, 3, 4],
        "pdg": [13, 13, 13, 13, -13],
        "target_z": [1, 8, 8, 8, 0],
    })
    w = hydrogen_binding_weights(ends)
    assert w[0] == 0.0
    assert abs(w[[0, 1, 2, 3]].sum() - 4.0) < 1e-12
    assert w[4] == 1.0
