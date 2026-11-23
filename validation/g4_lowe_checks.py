"""Low-energy physics checks of the wcprop propagation in water.

Generates test primaries, runs wcprop and compares with reference values:
electron light yield (calibration of the visible-energy proxy), mu+ lifetime and
Michel spectrum, mu- capture fraction and lifetime on O, pi+- at rest, and
neutron capture on H (time constant, gamma energy).

Usage:
    python validation/g4_lowe_checks.py <path/to/wcprop> [outdir]
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from hk_sunu_sim import g4out, testgen as tg

ELECTRON_E = [5, 10, 15, 20, 30, 40, 50, 60, 80]


def make_inputs(out, rng):
    tg.write(pd.concat([tg.single(11, e, 200, rng, first_id=i * 1000)
                        for i, e in enumerate(ELECTRON_E)]), out / "ecal.csv")
    tg.write(tg.single(13, rng.uniform(1, 50, 2000), 2000, rng), out / "mum.csv")
    tg.write(tg.single(-13, rng.uniform(1, 50, 2000), 2000, rng), out / "mup.csv")
    tg.write(tg.single(211, 30.0, 1000, rng), out / "pip.csv")
    tg.write(tg.single(-211, 30.0, 1000, rng), out / "pim.csv")
    tg.write(tg.ncqe_like(1000, rng), out / "ncqe.csv")


def run(wcprop, out, names, seed0=100):
    procs = [subprocess.Popen([str(wcprop), f"{n}.csv", n, "QGSP_BIC_HP", str(seed0 + i)], cwd=out,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
             for i, n in enumerate(names)]
    for p in procs:
        p.wait()


def main():
    wcprop = Path(sys.argv[1]).resolve()
    out = Path(sys.argv[2] if len(sys.argv) > 2 else "data/g4checks")
    out.mkdir(parents=True, exist_ok=True)
    names = ["ecal", "mum", "mup", "pip", "pim", "ncqe"]
    make_inputs(out, np.random.default_rng(2))
    run(wcprop, out, names)
    rows = []

    light, _, _ = g4out.read(out / "ecal")
    p = g4out.clusters(light).query("rank == 0").copy()
    p["E"] = np.array(ELECTRON_E)[p.event_id // 1000]
    m = p.groupby("E")["n"].mean()
    k, b = np.polyfit(m.index, m.values, 1)
    rows.append(("e- light yield [photons/MeV]", f"{k:.1f} (offset {b:.0f})", "calibration"))

    _, _, ends = g4out.read(out / "mup")
    d = ends.query("parent_id == 0")
    rows.append(("mu+ mean decay time [ns]", f"{d.t.mean():.0f} +- {d.t.std() / np.sqrt(len(d)):.0f}", "2197"))
    rows.append(("mu+ Michel <E>, max [MeV]", f"{d.e_charged.mean():.1f}, {d.e_charged.max():.1f}", "~37, 52.8"))

    _, _, ends = g4out.read(out / "mum")
    d = ends.query("parent_id == 0")
    o = d[d.target_z == 8]
    rows.append(("mu- stops bound to H (Geant4)", f"{(d.target_z == 1).mean():.3f}", "~0 (transfer to O)"))
    rows.append(("mu- capture fraction on O", f"{(o.fate == 'capture').mean():.3f}", f"{g4out.MU_MINUS_CAPTURE_O}"))
    rows.append(("mu- mean time on O [ns]", f"{o.t.mean():.0f} +- {o.t.std() / np.sqrt(len(o)):.0f}",
                 f"{g4out.MU_MINUS_TAU_O_NS:.0f}"))
    rows.append(("mu- capture on O: <neutrons>", f"{o[o.fate == 'capture'].n_neutron.mean():.2f}", "to check"))

    for n, ref in [("pip", "decay ~100 %"), ("pim", "absorption ~100 %")]:
        _, ncap, ends = g4out.read(out / n)
        d = ends.query("parent_id == 0")
        rows.append((f"{n} at rest: fates", str(d.fate.value_counts(normalize=True).round(3).to_dict()), ref))

    light, ncap, _ = g4out.read(out / "ncqe")
    h = ncap[ncap.target_z == 1]
    rows.append(("n capture on H: fraction", f"{len(h) / len(ncap):.3f}", "~1 (O: ~3e-4)"))
    rows.append(("n capture on H: tau from median [us]", f"{h.t.median() / 1e3 / np.log(2):.0f}", "~204"))
    rows.append(("n capture on H: E_gamma [MeV]", f"{h.e_gamma.mean():.3f}", "2.224"))

    print(pd.DataFrame(rows, columns=["check", "wcprop", "reference"]).to_string(index=False))


if __name__ == "__main__":
    main()
