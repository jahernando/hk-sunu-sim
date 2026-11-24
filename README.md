# hk-sunu-sim

Simulation of the low-energy backgrounds to the Diffuse Supernova Neutrino
Background (DSNB) search in a pure-water Hyper-Kamiokande tank:

1. **Cosmic-muon spallation**: isotope yields (⁹Li, ⁸He, ¹²B, …) per muon in
   infinite water, folded with the muon flux through the tank.
2. **Atmospheric neutrinos**: NCQE, invisible-muon decay electrons, ν_e CC and
   visible μ/π, generated per unit mass and propagated in infinite water.

The output is a set of event tables that the analysis package
[hk-sunu](https://github.com/jahernando/hk-sunu) reads to build
Hyper-Kamiokande background templates.

"Fast" refers to the detector: the response is parametrised (energy, vertex,
Cherenkov angle, neutron tagging) instead of a full optical simulation. The
particle transport (FLUKA, Geant4) is a full physics simulation.

This repository uses only public information. It is not an official product of
the Hyper-Kamiokande collaboration.

## Design

- **Infinite water, local coordinates.** Events are simulated in an infinite
  volume of pure water and stored in local coordinates: relative to the
  interaction vertex (atmospheric) or to the muon track (spallation). They are
  placed in the tank geometry afterwards, so the fiducial volume, geometry and
  cuts can change without re-running the transport.
- **Tables, not seeds.** Every production writes Parquet tables plus a metadata
  file (tool, version, physics list, seed, number of primaries, normalisation).
  Seeds are kept for provenance only.
- **Weights, not regeneration.** Atmospheric events carry flux × cross-section
  and oscillation weights, so fluxes (HKKM, FLUKA) and oscillation parameters
  can be changed by re-weighting.
- **Muon-energy grid.** Spallation is simulated on a grid of muon energies and
  weighted afterwards by the muon spectrum at SK (validation) or HK (projection).
- **Validate on SK-IV first.** Every component is checked against published
  Super-Kamiokande IV results before it is projected to Hyper-Kamiokande.

## Public and collaboration tools

This repository is public. It contains **only our own code and inputs**:

| Tool | Access | Role here |
|---|---|---|
| Geant4 | public | propagation in water, alternative spallation yields |
| GENIE | public | default atmospheric-neutrino generator |
| FLUKA | free, with registration | reference spallation yields |
| NEUT | Hyper-K / Super-K collaboration | optional generator backend |
| WCSim | public code | calibration of the parametrised response |

- The source code of NEUT, collaboration geometries, internal fluxes and
  internal data are **never** committed.
- NEUT-specific inputs that are not public go in
  `production/atmospheric/neut/private/`, which is ignored by git.
- Every step has a public backend, so the chain can run without collaboration
  access.
- Results obtained with collaboration tools or inputs fall under the
  Hyper-Kamiokande publication policy.

## Layout

| Path | Content |
|---|---|
| `src/hk_sunu_sim/` | Python package: table schemas and I/O, parametrised response, tank placement, spallation folding |
| `NB/` | Validation notebooks (see below) |
| `production/spallation/` | FLUKA and Geant4 inputs for muons in infinite water; muon flux at depth |
| `production/atmospheric/` | Flux tables, generator configurations (GENIE, NEUT), Geant4 propagation |
| `config/` | One YAML file per production |
| `env/` | Container definitions (Apptainer) with pinned versions |
| `slurm/` | Job-array templates for the cluster |
| `validation/` | Comparisons with SK-IV published results |
| `tests/` | Unit and consistency checks |
| `data/` | Local output, not tracked |

## Setup

```bash
pip install -e ".[dev,nb]"
pytest
```

## Notebooks

| Notebook | Content | Needs |
|---|---|---|
| `NB01_wcprop_lowe_validation` | Geant4 low-energy physics in water: e⁻ light yield and resolution, μ± at rest, π± at rest, neutron capture | `data/g4checks` from `validation/g4_lowe_checks.py` |
| `NB02_hkkm_flux` | HKKM 2014 flux at Kamioka: spectra, ratios, solar and zenith dependence, the gap below 100 MeV | HKKM tables in `data/flux/hkkm2014` |
| NB03 (planned) | Generator: ν / ν̄ spectra, interaction modes vs E_ν, lepton energy, neutron multiplicity, ¹⁶O de-excitation γ | generator events |
| NB04 (planned) | Final-state particles above / below the Cherenkov threshold, per event and per particle | generator + wcprop |
| NB05 (planned) | Response, candidates and comparison with SK-IV Table XIV | all of the above |

## Status

First milestone (atmospheric neutrinos in the SK-IV configuration, compared
with SK-IV Table XIV) in progress:

- [x] `wcprop`: Geant4 propagation in water, validated at low energy
  (`validation/g4_lowe_checks.py`)
- [x] HKKM 2014 flux reader (`flux.py`; tables downloaded from the authors' page)
- [x] parametrised response: energy, vertex, neutron tagging (`response.py`)
- [x] candidate building and truth categories (`events.py`)
- [ ] neutrino generator (GENIE / NEUT) and converter to wcprop input
- [ ] Cherenkov-angle proxy calibration
- [ ] comparison with SK-IV
