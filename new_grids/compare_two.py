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
    ratio=1.0,
    verbose=True,
    mass_divide=False,
    remove_points=False,
):
    dat = API.dataset(dataset_input={"dataset": dataset_vp}, theoryid=200, use_cuts="internal")
    pdf = API.pdf(pdf=pdf_name)
    lpdf = lhapdf.mkPDF(pdf_name)
    res_vp = central_predictions(dat, pdf)
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
        if mass_divide:
            # For some reason not all grids need that
            ratio_mass = grid.raw.bin_left(0) ** 3 / 38.8 ** 3
        else:
            ratio_mass = 1.0

        if remove_points:
            ratio_mass *= np.less(grid.raw.bin_left(0), 7.52) | np.greater(
                grid.raw.bin_left(0), 7.54
            )
        res_pine_partial.append(
            grid.convolute_with_one(2212, lpdf.xfxQ2, lpdf.alphasQ2) * ratio_mass
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
    if ratio != 1.0:
        all_res_corrected = pd.concat([res_vp, res_pine * ratio, res_vp / res_pine / ratio], axis=1)
        all_res_corrected.to_csv(
            f"out_csv/{dataset_vp}_corr.csv", sep="\t", header=["vp", f"pine*{ratio}", "ratio"]
        )
    if verbose and ratio != 1.0:
        print(all_res_corrected)
    elif verbose:
        print(all_res)

    if np.allclose(res_vp, res_pine * ratio, rtol=1e-2):
        print(f"The pienappl and NNPDF version are compatible!")
    else:
        print(f"XXX The pineappl and NNPDF version are NOT compatible!")
    return res_vp, res_pine


if __name__ == "__main__":
    #     a, b = run_comparison("DYE605", "E605nlo", ratio=54.35, verbose=False)
    #     a, b = run_comparison("DYE886P", "rescaled_E886Pnlo", ratio=54.35, mass_divide=False, verbose=True)
    a, b = run_comparison("DYE886P", "E886Pnlo", ratio=54.35, mass_divide=True, verbose=True)
#     a, b = run_comparison("DYE886R", ["E886deutRnlo", "E886Rnlo"], verbose=False)
#     a, b = run_comparison("DYE906_D", "E906deutnlo_bin_*", ratio=54.35, mass_divide=True, remove_points=False)
#     a, b = run_comparison("DYE906R", "E906*nlo_bin_*", mass_divide=False, remove_points=False)
#     i = 9
#     a, b = run_comparison("DYE906_D", f"E906deutnlo_bin_0{i}", ratio=54.35, mass_divide=True, remove_points=False)
# problems at: i=3,5
#         a, b = run_comparison("DYE906_D", f"E906nlo_bin_0{i}", ratio=54.35, mass_divide=True, remove_points=False)
