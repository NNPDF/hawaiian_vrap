# Hawaian Vrap

This program is an extension of [Vrap-v0.9](https://www.slac.stanford.edu/~lance/Vrap/). Vrap is a program that computes rapidity distributions for production of lepton-pairs via virtual photons, W or Z bosons at hadron colliders at NNLO in QCD. The program was originally developed by Frank Petriello and Lance Dixon. The original version of Vrap can be obtained in [Lance Dixon's website](https://www.slac.stanford.edu/~lance/Vrap/).

The code in this repository modifies version 0.9 of Vrap and it is redistributed with permission of the authors.

## Changes with respect to Vrap

The code in this repository adds and option to use isoscalar targets as a hadron type (`piso`) in addition to the original `pp` and `ppbar` options. It also interfaces the code with [pineappl](https://n3pdf.github.io/pineappl/), a library that produces fast-interpolation grids for fitting parton distribution functions.

:warning: The goal of these modifications is to generate `pineappl` grids for fixed-target Drell-Yan. These grids are used by the NNPDF Collaboration to determine parton distribution functions. As such, this use-case is the only one that has been tested and benchmarked (see [regression tests](https://github.com/scarlehoff/Hawaiian_vrap/tree/main/regression_test)).

## How to run

In order to run the code, clone the repository and run the following commands:

```bash
git clone https://github.com/NNPDF/Hawaiian_vrap.git
cd Hawaiian_vrap
cd src && autoreconf -fiv && cd ..
mkdir -p build && cd build
../src/configure --prefix=$PWD
make -j4
./Vrap ../regression_test/inputE605nlo.dat 7 0.2 # needs to have NNPDF40_nnlo_as_01180 installed
```

## Disclaimer

As warned in Lance Dixon's webpage, for fixed-target kinematics the results might be slightly off due to an integration problem at high $\tau=\frac{M^2}{s}$.

> For large values of M^2/s, as typically encountered in fixed-target production, the version of the program given here is a bit off, due to an integration issue.

A cut on $\tau$ is typically enforced on the data to avoid this issue.

## Citations

The numerical program, "Vrap", is based on the paper
Phys.Rev.D69:094008,2004 [hep-ph/0312266]
by C. Anastasiou, L. Dixon, K. Melnikov and F. Petriello.

Please cite this paper whenever the "Vrap" program is used.
