#!/usr/bin/env python3
"""
    The datapoints for E886 have a factor of M**3 / s**3 with respect to E605.
    Since vrap uses the same code for both, one of the two needs to be rescaled.

    Note that only E886P needs to be rescaled since E886R is a ratio and the factor cancels out.
    The same is true for the k-factors and therefore those are correct even if they were
    computed in an inconsistent manner between E605 and E886.
"""
import numpy as np
from pineappl.bin import BinRemapper
from pineappl.grid import Grid
import lhapdf

lhapdf.setVerbosity(0)
lpdf = lhapdf.mkPDF("NNPDF40_nnlo_as_01180")

input_fktable = "E886Pnlo.pineappl.lz4"
output_fktable = "rescaled_E886Pnlo.pineappl.lz4"
comE = 38.8

grid = Grid.read(input_fktable)
res = grid.convolute_with_one(2212, lpdf.xfxQ2, lpdf.alphasQ2)
bl1 = grid.raw.bin_left(0)
bl2 = grid.raw.bin_left(1)
br1 = grid.raw.bin_left(0)
br2 = grid.raw.bin_left(1)
all_bins = []
for ba, bb in zip(list(zip(bl1, br1)), list(zip(bl2, br2))):
    all_bins.append(ba)
    all_bins.append(bb)

ratio_mass = bl1 ** 3 / comE ** 3

rebin = BinRemapper(1 / ratio_mass, all_bins)
grid.set_remapper(rebin)
new_res = grid.convolute_with_one(2212, lpdf.xfxQ2, lpdf.alphasQ2)

# Check that the values are correct before writing
np.testing.assert_allclose(res * ratio_mass, new_res)
grid.write_lz4(output_fktable)
