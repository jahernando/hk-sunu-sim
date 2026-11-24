import pytest

pytest.importorskip("matplotlib")
from matplotlib.colors import to_rgb  # noqa: E402

from hk_sunu_sim import plotting  # noqa: E402


def test_palette_tables_align():
    n = len(plotting.SERIES)
    assert len(plotting.LINESTYLES) == len(plotting.MARKERS) == len(plotting.HATCHES) == n
    assert len(set(plotting.SERIES)) == n and len(set(map(str, plotting.LINESTYLES))) == n
    assert len(plotting.BACKGROUND_FILLS) <= n
    for colour in plotting.SERIES + plotting.BACKGROUND_FILLS:
        to_rgb(colour)


def test_style_always_sets_linestyle():
    assert plotting.style(2) == {"color": plotting.SERIES[2], "ls": plotting.LINESTYLES[2]}
    assert plotting.style(1, ls="-", marker=True) == {"color": plotting.SERIES[1], "ls": "-", "marker": plotting.MARKERS[1]}
