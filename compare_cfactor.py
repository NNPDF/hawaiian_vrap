#!/usr/bin/env python3
"""
    Script to compare two cfactors
"""
import sys
from pathlib import Path
from argparse import ArgumentParser, ArgumentTypeError

import numpy as np
import pandas as pd

pd.set_option("display.max_rows", None)


def exiting_path(value):
    val = Path(value)
    if not val.exists():
        raise ArgumentTypeError(f"The path {value} doesn't exist")
    return val


def read_cfac(cfac_path):
    """Try to read a vrap cfactor, if it fails it means that this is an NNPDF cfactor"""
    try:
        cfac_val = np.loadtxt(cfac_path)
    except ValueError:
        print(f"This is an NNPDF cfactor, right? ({cfac_path}) skipping the 9 first rows")
        cfac_val, _ = np.loadtxt(cfac_path, skiprows=9, unpack=True)
    return cfac_val


def show_ratios(values, ratio_col=-1):
    """Generate a df of the ratios of values with respect to the last val
    (by default) and print it as a df"""
    new_df = pd.DataFrame([i / values[ratio_col] for i in values])
    print(new_df.T)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("cfactor_paths", type=exiting_path, nargs="+")
    parser.add_argument(
        "-r", "--ratios", help="Show ratio with respect to the last column", action="store_true"
    )
    args = parser.parse_args()

    cfactor_vals = [read_cfac(cpath) for cpath in args.cfactor_paths]

    df = pd.DataFrame(cfactor_vals).T

    if args.ratios:
        show_ratios(cfactor_vals)
    else:
        print(df)
        if sys.flags.interactive:
            show_me = lambda: show_ratios(cfactor_vals)
            print("Run `show_me()` if you want to see the ratios as well")
