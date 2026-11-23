"""Column schemas of the production tables.

Units: MeV, cm, ns, g/cm2. Positions are local:

- spallation tables are in the muon-track frame: ``s`` along the track from the
  start of the recording segment, ``lt`` transverse distance, ``phi`` azimuth;
- atmospheric tables are relative to the interaction vertex.

Each table is a Parquet file; ``batch`` identifies the production batch whose
metadata (tool, version, physics list, seed, number of primaries,
normalisation) is stored alongside.
"""

import pyarrow as pa

# --- Spallation -------------------------------------------------------------

MUONS = pa.schema([
    ("batch", pa.int32()),
    ("muon_id", pa.int64()),
    ("e_mu", pa.float64()),          # muon energy at the segment start [MeV]
    ("segment_length", pa.float64()),  # recorded track length [cm]
    ("edep_total", pa.float64()),    # energy deposited along the segment [MeV]
    ("dedx_profile", pa.list_(pa.float32())),  # deposit per bin along s [MeV]
    ("dedx_bin", pa.float32()),      # bin width of dedx_profile [cm]
    ("n_ncapture", pa.int32()),      # neutron captures from this muon
])

ISOTOPES = pa.schema([
    ("batch", pa.int32()),
    ("muon_id", pa.int64()),
    ("z", pa.int16()),
    ("a", pa.int16()),
    ("isomer", pa.int16()),          # 0 = ground state
    ("s", pa.float64()),             # [cm]
    ("lt", pa.float64()),            # [cm]
    ("phi", pa.float32()),           # [rad]
    ("t", pa.float64()),             # production time after muon [ns]
    ("parent_pdg", pa.int32()),      # particle that produced the isotope
    ("process", pa.string()),        # producing process, tool naming
])

NCAPTURES = pa.schema([
    ("batch", pa.int32()),
    ("muon_id", pa.int64()),
    ("s", pa.float64()),
    ("lt", pa.float64()),
    ("phi", pa.float32()),
    ("t", pa.float64()),             # capture time after muon [ns]
    ("target_z", pa.int16()),        # 1 = H, 8 = O
    ("e_gamma", pa.float32()),       # total gamma energy [MeV]
])

# --- Atmospheric neutrinos --------------------------------------------------

ATM_EVENTS = pa.schema([
    ("batch", pa.int32()),
    ("event_id", pa.int64()),
    ("nu_pdg", pa.int32()),
    ("e_nu", pa.float64()),          # [MeV]
    ("cos_zenith", pa.float32()),
    ("target_z", pa.int16()),        # 1 = H, 8 = O
    ("mode", pa.string()),           # generator interaction mode
    ("category", pa.string()),       # nue_cc, decay_e, mupi, ncqe, nc_other, other
    ("w_gen", pa.float64()),         # 1 / generation density, per unit mass
    ("w_flux", pa.float64()),        # flux x cross section, nominal flux
    ("w_osc", pa.float64()),         # oscillation weight, nominal parameters
    ("e_vis", pa.float64()),         # Cherenkov-weighted visible energy proxy [MeV]
    ("n_ncapture", pa.int32()),
    ("n_michel", pa.int32()),
    ("extent", pa.float32()),        # max distance of any deposit from vertex [cm]
])

ATM_PARTICLES = pa.schema([
    ("batch", pa.int32()),
    ("event_id", pa.int64()),
    ("pdg", pa.int32()),
    ("ekin", pa.float64()),          # [MeV]
    ("x", pa.float32()), ("y", pa.float32()), ("z", pa.float32()),  # start [cm]
    ("ux", pa.float32()), ("uy", pa.float32()), ("uz", pa.float32()),
    ("t", pa.float64()),             # [ns]
    ("above_cherenkov", pa.bool_()),
    ("origin", pa.string()),         # primary, michel, deexcitation, secondary
])

ATM_NCAPTURES = pa.schema([
    ("batch", pa.int32()),
    ("event_id", pa.int64()),
    ("x", pa.float32()), ("y", pa.float32()), ("z", pa.float32()),
    ("t", pa.float64()),
    ("target_z", pa.int16()),
    ("e_gamma", pa.float32()),
])

TABLES = {
    "muons": MUONS,
    "isotopes": ISOTOPES,
    "ncaptures": NCAPTURES,
    "atm_events": ATM_EVENTS,
    "atm_particles": ATM_PARTICLES,
    "atm_ncaptures": ATM_NCAPTURES,
}
