from glob import glob
import numpy as np
import pandas as pd
import lhapdf
from pineappl.grid import Grid
from validphys.api import API
from validphys.convolution import central_predictions

import NNPDF, lhapdf

NNPDF.SetVerbosity(0)
lhapdf.setVerbosity(0)


def run_comparison(
    dataset_vp,
    dataset_pin,
    pdf_name="NNPDF40_nnlo_as_01180",
    E906_factor=False,
    verbose=True,
    remove_points=False,
):
    dat = API.dataset(dataset_input={"dataset": dataset_vp}, theoryid=200, use_cuts="internal")
    pdf = API.pdf(pdf=pdf_name)
    lpdf = lhapdf.mkPDF(pdf_name)
    res_vp = central_predictions(dat, pdf)

    if E906_factor:
        # Fix the result computed by vp for E906 with (Energy_605 / Energy_906)
        res_vp *= (38.8/15.06)**2

    if isinstance(dataset_pin, str):
        if "*" in dataset_pin:
            dataset_pin = [
                i.replace(".pineappl.lz4", "") for i in sorted(glob(dataset_pin + ".pineappl.lz4"))
            ]
        else:
            dataset_pin = [dataset_pin]

    res_pine_partial = []
    for d in dataset_pin:
        grid = Grid.read(f"{d}.pineappl.lz4")

        remove_factor = 1.0
        if remove_points:
            remove_factor *= np.less(grid.raw.bin_left(0), 7.52) | np.greater(
                grid.raw.bin_left(0), 7.54
            )
        res_pine_partial.append(
            grid.convolute_with_one(2212, lpdf.xfxQ2, lpdf.alphasQ2) * remove_factor
        )

    res_pine_raw = pd.DataFrame(res_pine_partial[0])
    if len(res_pine_partial) == 2:
        res_pine_raw /= pd.DataFrame(res_pine_partial[1])
    elif len(res_pine_partial) == 10:
        # Two wrong bins:
        #         res_pine_partial[3][1] = 0.0
        #         res_pine_partial[5][0] = 0.0
        ##############
        res_pine_raw = pd.DataFrame(np.sum(res_pine_partial, axis=0))
    elif len(res_pine_partial) == 20:
        res_pine_raw = pd.DataFrame(np.sum(res_pine_partial[:10], axis=0)) / pd.DataFrame(
            np.sum(res_pine_partial[10:], axis=0)
        )

    # And now put them both side by side
    res_pine = res_pine_raw.loc[res_vp.index]
    all_res = pd.concat([res_vp, res_pine, res_vp / res_pine], axis=1)
    all_res.to_csv(f"out_csv/{dataset_vp}.csv", sep="\t", header=["vp", "pine", "ratio"])
    if verbose:
        print(all_res)

    if np.allclose(res_vp, res_pine, rtol=1.6e-2):
        print(f"The pineappl and NNPDF version for {dataset_vp} are compatible!")
    else:
        print(f"XXX The pineappl and NNPDF version for {dataset_vp} are NOT compatible!")
    return res_vp, res_pine


if __name__ == "__main__":
    a, b = run_comparison("DYE605", "E605nlo", verbose=False)
    a, b = run_comparison("DYE886P", "E866nlo", verbose=False)
    a, b = run_comparison("DYE886R", ["E866deutRnlo", "E866Rnlo"], verbose=False)
#     a, b = run_comparison("DYE906_D", "E906deutnlo_bin_*", E906_factor=True, verbose=False)
    a, b = run_comparison("DYE906R", "E906*nlo_bin_*", verbose=False)
