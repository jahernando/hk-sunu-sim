"""Atmospheric neutrino flux tables (HKKM 2014, Honda et al., PRD 92, 023004).

Tables are downloaded from the authors' public page and are not redistributed:

    http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/index.html

Files used: ``kam-ally-20-01-sol{min,max}.d.gz`` (Kamioka, 20 bins in cos(zenith),
azimuth-averaged). Each block has E_nu [GeV] and the flux of nu_mu, anti-nu_mu,
nu_e, anti-nu_e in (m^2 s sr GeV)^-1, from 0.1 GeV to 10 TeV.

Below 0.1 GeV the tables are not defined. ``flux`` returns 0 there unless
``extrapolate=True``, which continues the log-log slope of the first two points:
a placeholder, not a prediction (use a dedicated low-energy flux, e.g. the FLUKA
calculation of Battistoni et al. 2005, for systematic studies).
"""

import gzip
import re
from pathlib import Path

import numpy as np

COLUMNS = {14: 0, -14: 1, 12: 2, -12: 3}
DEFAULT_DIR = Path(__file__).resolve().parents[2] / "data" / "flux" / "hkkm2014"
URL = "http://www-rccn.icrr.u-tokyo.ac.jp/mhonda/public/nflx2014/"


class HKKMTable:
    """Flux table in cos(zenith) bins; values in cm^-2 s^-1 sr^-1 MeV^-1."""

    def __init__(self, path):
        self.path = Path(path)
        opener = gzip.open if self.path.suffix == ".gz" else open
        blocks, bounds, current = [], [], None
        with opener(self.path, "rt") as f:
            for line in f:
                m = re.search(r"cosZ\s*=\s*([-\d.]+)\s*--\s*([-\d.]+)", line)
                if m:
                    bounds.append((float(m.group(1)), float(m.group(2))))
                    current = []
                    blocks.append(current)
                    continue
                parts = line.split()
                if current is not None and parts and parts[0][0].isdigit():
                    current.append([float(x) for x in parts[:5]])
        data = np.array(blocks)                       # (ncos, nE, 5)
        order = np.argsort([b[0] for b in bounds])    # increasing cos(zenith)
        data = data[order]
        self.cos_edges = np.array([bounds[i][0] for i in order] + [bounds[order[-1]][1]])
        self.e = data[0, :, 0] * 1e3                  # MeV
        # (m^2 s sr GeV)^-1 -> (cm^2 s sr MeV)^-1
        self.table = data[:, :, 1:] * 1e-4 * 1e-3      # (ncos, nE, 4)

    def flux(self, pdg, e_mev, cosz, extrapolate=False):
        """d Phi / dE dOmega at energies e_mev [MeV] and cos(zenith) cosz."""
        e = np.atleast_1d(np.asarray(e_mev, float))
        cz = np.broadcast_to(np.asarray(cosz, float), e.shape)
        icos = np.clip(np.searchsorted(self.cos_edges, cz, side="right") - 1, 0, len(self.cos_edges) - 2)
        logf = np.log(self.table[:, :, COLUMNS[pdg]])
        loge, loge_tab = np.log(e), np.log(self.e)
        out = np.empty_like(e)
        for i in np.unique(icos):
            sel = icos == i
            out[sel] = np.exp(np.interp(loge[sel], loge_tab, logf[i]))
            low = sel & (e < self.e[0])
            if extrapolate:
                slope = (logf[i, 1] - logf[i, 0]) / (loge_tab[1] - loge_tab[0])
                out[low] = np.exp(logf[i, 0] + slope * (loge[low] - loge_tab[0]))
            else:
                out[low] = 0.0
        out[e > self.e[-1]] = 0.0
        return out

    def flux_solid_angle_integrated(self, pdg, e_mev):
        """d Phi / dE integrated over all directions [cm^-2 s^-1 MeV^-1]."""
        dcos = np.diff(self.cos_edges)
        mids = 0.5 * (self.cos_edges[1:] + self.cos_edges[:-1])
        return sum(2 * np.pi * dc * self.flux(pdg, e_mev, c) for dc, c in zip(dcos, mids))


def load(solar="min", directory=DEFAULT_DIR):
    return HKKMTable(Path(directory) / f"kam-ally-20-01-sol{solar}.d.gz")
