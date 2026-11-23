"""Parametrised water-Cherenkov detector response.

Input is the light produced by wcprop (mean Cherenkov photons per time cluster).
The chain is:

1. photon detection: n_pe ~ Poisson(eps * n_photons), with eps such that an
   electron gives ``pe_per_mev`` photoelectrons per MeV (SK-IV: ~6 effective
   hits per MeV);
2. energy: E_rec = n_pe / (eps * k), with k the photons per MeV of electrons
   (calibrated with wcprop, ``validation/g4_lowe_checks.py``); the electron
   resolution then comes from photostatistics plus the light-yield fluctuations
   already in Geant4;
3. vertex: Gaussian smearing with sigma(E) = a / sqrt(E / 10 MeV) + b;
4. neutron tagging: each true capture in the time window is tagged with
   probability ``eps_tag``; fake tags are Poisson with mean ``mistag``.

All parameters are detector configurations (``SK4``, ``HK``); the values are
first guesses to be calibrated against published SK-IV numbers.
"""

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

PHOTONS_PER_MEV = 165.2   # wcprop electrons, n = 1.333, 300-600 nm (g4_lowe_checks)


@dataclass(frozen=True)
class Detector:
    name: str
    pe_per_mev: float          # photoelectrons (effective hits) per MeV of electron
    vertex_a_cm: float         # vertex resolution at 10 MeV (stochastic part)
    vertex_b_cm: float         # constant part
    ntag_window_us: tuple      # neutron search window after the prompt event
    eps_tag: float             # tagging efficiency of a capture on H in the window
    mistag: float              # mean fake tags per event
    photons_per_mev: float = PHOTONS_PER_MEV


# SK-IV DSNB analysis (arXiv:2109.11174): window 18-535 us; the tagging point is
# the one quoted for the spectral fit (~25 % efficiency, ~1 % mistag per event):
# values to be checked against the paper before validation.
SK4 = Detector("SK-IV", pe_per_mev=6.0, vertex_a_cm=50.0, vertex_b_cm=10.0,
               ntag_window_us=(18.0, 535.0), eps_tag=0.25, mistag=0.01)

# HK: same effective light yield assumed (20 % coverage with higher-QE PMTs);
# tagging scenarios A/B/C are set by replacing eps_tag and mistag.
HK = replace(SK4, name="HK")


def energy(n_photons, det, rng):
    """Reconstructed energy [MeV] from mean Cherenkov photons."""
    eps = det.pe_per_mev / det.photons_per_mev
    npe = rng.poisson(eps * np.asarray(n_photons, float))
    return npe / det.pe_per_mev


def smear_vertex(xyz, e_rec, det, rng):
    """Gaussian vertex smearing; xyz is (N, 3) in cm."""
    e = np.maximum(np.asarray(e_rec, float), 1.0)
    sigma = det.vertex_a_cm / np.sqrt(e / 10.0) + det.vertex_b_cm
    return np.asarray(xyz, float) + rng.normal(size=(len(e), 3)) * sigma[:, None]


def neutron_tags(t_prompt_ns, ncap, det, rng):
    """Number of tagged neutrons per candidate.

    t_prompt_ns: Series indexed by event_id with the candidate time.
    ncap: wcprop capture table (event_id, t, target_z, ...).
    Returns a Series indexed like t_prompt_ns.
    """
    lo, hi = (w * 1e3 for w in det.ntag_window_us)
    c = ncap.merge(t_prompt_ns.rename("t0"), left_on="event_id", right_index=True)
    dt = c["t"] - c["t0"]
    c = c[(dt >= lo) & (dt <= hi)]
    tagged = c[rng.random(len(c)) < det.eps_tag].groupby("event_id").size()
    n = tagged.reindex(t_prompt_ns.index, fill_value=0)
    return n + rng.poisson(det.mistag, len(n))
