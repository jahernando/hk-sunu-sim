"""Shared, colour-blind-safe plotting style for the notebooks.

Copy of ``hk_sunu.plotting`` (kept identical so both repositories draw the
same way, without making hk-sunu a dependency).

Colours follow Okabe & Ito (2008), distinguishable under protan, deutan and
tritan vision. Colour is never the only cue: each series index also has its own
line style, marker and hatch, so plots stay readable in greyscale.

Usage in a notebook::

    from hk_sunu_sim.plotting import INK, MUTED, GRID, SERIES, style, stair, use_style
    use_style()
    ax.plot(x, y, **style(0), label="...")
"""

import matplotlib.pyplot as plt
from cycler import cycler

INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e6e5df"

# Okabe-Ito: blue, orange, bluish green, reddish purple, sky blue, vermillion, yellow.
# Yellow goes last: low contrast on a white background.
SERIES = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442")
LINESTYLES = ("-", "--", "-.", ":", (0, (3, 1, 1, 1, 1, 1)), (0, (8, 3)), (0, (1, 3)))
MARKERS = ("o", "s", "^", "D", "v", "P", "X")
HATCHES = ("", "//", "..", "xx", "\\\\", "--", "++")

# Neutral fills for stacked background components (paired with HATCHES).
BACKGROUND_FILLS = ("#c9c7bd", "#a6a499", "#7b7a72", "#e6e5df", "#d9b77a")


def style(i, ls=None, marker=False):
    """Line keywords for series ``i``: colour, line style and optionally marker.

    ``ls=None`` gives the series' own line style. Pass an explicit style (e.g.
    ``"-"``) when line style encodes something else and the marker identifies
    the series; always set, so the default cycle never picks one silently.
    """
    kwargs = {"color": SERIES[i], "ls": LINESTYLES[i] if ls is None else ls}
    if marker:
        kwargs["marker"] = MARKERS[i]
    return kwargs


def fill(i, background=False):
    """Bar keywords for stacked component ``i``: face colour and hatch."""
    colors = BACKGROUND_FILLS if background else SERIES
    return {"color": colors[i], "hatch": HATCHES[i]}


def stair(ax, values, bin_edges, **kwargs):
    """Histogram outline without the vertical edges down to zero."""
    if "linewidth" not in kwargs:
        kwargs.setdefault("lw", 2.0)  # patches default to 1 pt, too thin for dash patterns
    return ax.stairs(values, bin_edges, baseline=None, **kwargs)


def use_style():
    """Set the notebook rcParams; the default cycle pairs colour and line style."""
    plt.rcParams.update({
        "figure.figsize": (7.2, 4.4), "figure.dpi": 110, "font.size": 11,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
        "lines.linewidth": 2.0, "lines.markersize": 5.5, "legend.frameon": False,
        "legend.handlelength": 3.0,  # long enough to show dash patterns
        "hatch.linewidth": 0.8,
        "axes.prop_cycle": cycler(color=SERIES) + cycler(linestyle=LINESTYLES),
    })
