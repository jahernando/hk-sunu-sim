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

# Light emission from wcprop: photon-weighted sums per 10 ns bin (see wcprop.cc)
LIGHT = pa.schema([
    ("batch", pa.int32()),
    ("event_id", pa.int64()),
    ("tbin", pa.int64()),            # floor(t / 10 ns)
    ("n", pa.float64()),             # mean Cherenkov photons, 300-600 nm
    ("sx", pa.float64()), ("sy", pa.float64()), ("sz", pa.float64()),     # sum n*x [cm]
    ("sux", pa.float64()), ("suy", pa.float64()), ("suz", pa.float64()),  # sum n*u
    ("st", pa.float64()),            # sum n*t [ns]
])

# End of stopped/killed mu+- and pi+- (truth)
ENDS = pa.schema([
    ("batch", pa.int32()),
    ("event_id", pa.int64()),
    ("track_id", pa.int32()),
    ("parent_id", pa.int32()),
    ("pdg", pa.int32()),
    ("process", pa.string()),
    ("fate", pa.string()),           # decay, capture, absorption, inelastic, other
    ("t", pa.float64()),             # time of the end products [ns]
    ("x", pa.float32()), ("y", pa.float32()), ("z", pa.float32()),
    ("e_charged", pa.float32()),     # kinetic energy of charged end products [MeV]
    ("n_neutron", pa.int16()),
    ("n_gamma", pa.int16()),         # includes muonic/pionic X-rays
    ("e_gamma", pa.float32()),
    ("n_nu", pa.int16()),
    ("target_z", pa.int16()),        # nucleus that took a stopped mu-/pi- (0 if none)
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
    "light": LIGHT,
    "ends": ENDS,
    "atm_ncaptures": ATM_NCAPTURES,
}
