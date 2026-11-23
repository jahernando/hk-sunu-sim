# wcprop

Geant4 propagation of primary particles in an infinite volume of water.

Optical photons are not tracked. For every charged step the mean Cherenkov
yield (constant n = 1.333, 300–600 nm) is accumulated in 10 ns bins, with
photon-weighted position, direction and time. The parametrised detector response
turns this light into visible energy, vertex and a Cherenkov-angle proxy.
Neutron captures and the end of stopped μ± and π± are stored as truth.

Physics: reference list (default `QGSP_BIC_HP`, which already includes radioactive
decay) plus `G4ThermalNeutrons` with the S(α,β) data of hydrogen in water.

## Build and run

```bash
mkdir build && cd build
cmake -DGeant4_DIR=$(geant4-config --prefix)/lib/Geant4-<version> ..
make
./wcprop primaries.csv out QGSP_BIC_HP <seed> [tmax_ns]
```

Input: `event_id,pdg,px,py,pz` (MeV/c), vertex at the origin, t = 0.
Outputs: `out_light.csv`, `out_ncap.csv`, `out_ends.csv`; columns are described
at the top of `wcprop.cc` and in `hk_sunu_sim.schema`.

## Known issue: μ⁻ and π⁻ bound to hydrogen

Geant4 11.0 binds ~30 % of stopped μ⁻ (and π⁻) in water to hydrogen and does not
transfer the μp (πp) atom to oxygen. In real water the transfer happens, so all
μ⁻ end on O:

| | Geant4 11.0 | Physical |
|---|---|---|
| μ⁻ stops on H | 0.30 | ≈ 0 |
| μ⁻ capture fraction | 0.12 overall, 0.185 on O | 0.184 |
| μ⁻ mean time on O | 1785 ns | 1795 ns |

The nucleus is stored (`target_z` in the ends table), and
`hk_sunu_sim.g4out.hydrogen_binding_weights` gives event weights that drop the
H-bound stops and rescale the O-bound ones. On O, Geant4 reproduces the capture
fraction and lifetime.

## Checks

`validation/g4_lowe_checks.py` runs the low-energy checks (electron light
yield, μ± and π± at rest, neutron capture on H).
