#!/usr/bin/env python3
"""
    Uses the raw data found in buildmaster to compute a kinematics file for E906 for Vrap
"""
from glob import glob
from pathlib import Path
import numpy as np

masses, _ = np.loadtxt("./input_kinematics/E906R.dat", unpack=True)
new_kin = Path("./input_kinematics/new_E906R.dat")

buildmaster_path = "."
_, _, _, _, _, all_pt, _, _, _, _ = np.loadtxt(f"{buildmaster_path}/data_paper.dat", unpack=True)
bins = sorted(glob(f"{buildmaster_path}/DYE906R_BIN*"))


def compute_rapidity(xb, xt, M=None, pT=None, simplified=False):
    if simplified:
        return 0.5 * np.log(xb / xt)
    if M is None or pT is None:
        raise Exception("oh no!")
    beam_energy = 120
    proton_mass = 0.938
    s = 2.0 * proton_mass ** 2 + 2.0 * beam_energy * proton_mass
    tau = M ** 2 / s
    xF = xb - xt
    pt = pT / M
    return 0.5 * np.log(
        (np.sqrt(xF ** 2 + 4 * tau * (1 + pt ** 2)) + xF)
        / (np.sqrt(xF ** 2 + 4 * tau * (1 + pt ** 2)) - xF)
    )


index = 0
ret = []
for bin_file in bins:
    xbxt = np.loadtxt(bin_file)
    for pT, (xb, xt) in zip(all_pt, xbxt):
        m = masses[index]
        y = compute_rapidity(xb, xt, simplified=True, M=m, pT=pT)
        if m == 7.53:
            y = 0.0
        ret.append(f"{m:.8f} {y:.8f}")
        index += 1

new_kin.write_text("\n".join(ret))
