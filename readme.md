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

---

Note: the documentation below is an exact copy of that in [Lance Dixon's website](https://www.slac.stanford.edu/~lance/Vrap/doc.txt)

============= Introduction ===========================================

The numerical program, "Vrap", is based on the paper 
Phys.Rev.D69:094008,2004 [hep-ph/0312266]
by C. Anastasiou, L. Dixon, K. Melnikov and F. Petriello.

It can produce rapidity distributions for lepton pairs produced
at hadron colliders via gamma^*, Z and W intermediate vector bosons, 
at LO, NLO, and NNLO in perturbative QCD.
 
Electroweak and QED corrections are not included.
Some of the short-distance EW corrections can be approximated
by running alpha_QED up to around the Z.

[Should one want to change the quark/lepton couplings to that of a
particular beyond-the-Standard-Model Z' or W', that can be done
within "Vlumifns.C".]

The lepton-pair mass M in the paper is generally referred to as Q in the code.
The lepton-pair rapidity Y is generally referred to as y in the code.

The typical quantity produced has the form (as in hep-ph/0312266)
d^2sigma/dM/dY [pb/GeV].
If for example you run at the Z peak (Q = m_Z),
you will get d^2sigma/dM/dY |_{M = m_Z} [pb/GeV].
If you want to work in the narrow-resonance approximation,
you can just multiply the result by pi/2 * Gamma_Z = 3.919 GeV.

But it is also possible (and may be more accurate) to integrate over 
a range in M (in GeV), (see the "int_M" options).  
Then the output has the form
d sigma/dY [pb].

Comments on the code or the documentation are welcome.

=========== Unpacking & Compiling =======================================

The program is written in C++, but it does not use much in the way of classes.

The uu-encoded .tar.gz file is at:

http://www.slac.stanford.edu/~lance/Vrap.uu

Run csh Vrap.uu to unpack it.

To compile it, use

$ make Vrap

I have compiled it under linux (g++).  There may be problems with at
least one linux compiler, however, so let me know if yours complains.

=========== Usage =======================================

There are a lot of files, but for the most part, the only one you
should have to touch is "Vrap.C", in particular the "control" section,
beginning after

"//===========  MAIN PROGRAM  =========================================="

And much of the control can be implemented by commenting/uncommenting
various lines there.  (A new, improved user interface may be coming
soon, courtesy of the CMS group at ETH Zurich, Guenther Dissertori,
Andre Holzner, and company.  I have been lazy in not incorporating their
improvements yet....)

The first few lines after "MAIN PROGRAM" choose the collider.
Q sets the lepton pair invariant mass.
But this is over-ridden by the "scan mass" option (see below).

We use "Nf = 5" light quark flavors for all the non-fixed-target runs.

The line 
muF = 0.5*Q;   muR = muF;  // ren. & fact. scales
is an example of how to set the renormalization & factorization scales,
independently if you like.

In case you want to isolate some of the contributions of individual
partonic channels to the cross section, 
replace the command "compute_all()" by "compute_qg()" 
for the quark-gluon channel only, etc.

The command
"setV(Zgamma,Q,alphat,Nf,0);"
determines that the computation is for both Z + gamma^* together, 
loads Q in, as well as "alphat" = alpha_QED(M_Z) 
("alphat" should just affect overall normalization 
--- alphat = 1./128. is a "stand-in" for the fact 
that no true EW corrections are included here).

Besides "Zgamma", you can also use "gamma_only",
"Z_only", "Zgamma_interf" (just the interference terms
between Z and gamma exchange),
and in the charged-current case, "Wplus" and "Wminus".

Then the next few lines are where the all-important
"order" of the calculation is set.
order_flag = 0,1,2    corresponds respectively to    LO,NLO,NNLO
in QCD perturbation theory.

Typically, in the same line the strong coupling constant, 
alpha_s(M_Z) = alpha_s_Z,
and the parton distribution function (pdf) set are initialized.
That's because for consistency a given pdf set requires 
a particular value of alpha_s(M_Z), used in its Q^2 evolution.

For example, to set the program to run the leading-order (LO) computation
with LO MRST2002 pdfs (which runs very quickly), use the line
"order_flag = 0;  alpha_s_Z = 0.130;  pdf_init(mrst,1);"

To go to NNLO MRST2002 (which might take several hours), use instead
"order_flag = 2;  alpha_s_Z = 0.1155;  pdf_init(mrst,6);"

The more recent MRST2004 pdfs are obtained as follows:

 // MRST2004 NLO:
 order_flag = 1;  alpha_s_Z = 0.120;  pdf_init(mrst,15);  

 // MRST2004 NNLO:
 order_flag = 2;  alpha_s_Z = 0.1167;  pdf_init(mrst,16);  

This is followed by a way to get the Alekhin02 distributions,
and some "non-standard" MRST options.

The next block of code is typing out the parameter settings
to try to remind you what you did when you read the output.

Next there are a bunch of options for scanning in rapidity y.

"sym" commands, as in "sym_scan_all", can be used for distributions
known to be symmetric in y, thus saving you half the time.
Basically most distributions are symmetric in y <-> -y,
except for W bosons at a p-pbar collider like the Tevatron; 
use "asym_scan_all" in that case.

These scans in rapidity produce d^2sigma/dM/dY curves like in the
plot, for different values of the renormalization and factorization
(muR,muF) scales.  "sym_scan_all" does all orders at once
(thus overwriting the order_flag, etc., chosen earlier)
It does it at two different values of muR = muF, namely mu1*Q and mu2*Q.
In this way you can produce "error" bands for the rapidity distributions.

You can also scan the renormalization and factorization (muR,muF) scales
at any fixed rapidity with "scan_mu"  (Various scans are possible with
slight modifications inside "scan_mu".)

"sym_scan_rap_y_M" does an additional integral over the lepton-pair 
mass, Q=M, in the symmetric case, and
"asym_scan_rap_y_M()" can be used for the asymmetric case.

Included in the package is a companion program, "Vtot", which 
implements the total cross section formulae of Hamberg, van Neerven
and Matsuura (Nucl. Phys. B359:343-405,1991, Erratum-ibid. B644:403-404,2002)
and uses very similar parameter controls as "Vrap".
This program can be quite useful for cross-checking that
d sigma/dY integrates back up to sigma under various conditions
(bug checking...)

=========== Warnings =======================================

The NLO and NNLO computations (order_flag = 0,1)
will sometimes spit out the phrase, 
"Born-type integral = ..."
part way through the computation.
This is just a piece of the full result, and should not be
used for anything.

Note that the file "QCDbasics.h" contains m_W, and "Vlumnifns.C" 
has mW, Gamma_W, and Br_l_W.
These parameters control the location, width and strength of the
W resonance (similarly for the Z).  But "Q" in the main program
is what controls the lepton-pair mass requested, that is,
which point on the resonance curve is being computed.
