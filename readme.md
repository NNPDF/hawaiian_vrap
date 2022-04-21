## How to run

```bash
mkdir -p build && cd build
../src/configure --prefix=$PWD
make -j4
./Vrap ../src/input.dat # needs to have NNPDF31_nnlo_as_0118 installed
```

The older standard we can use with more recent versions of `lhapdf` is `c++11` so some changes have been necessary in order to get `vrap` to compile: https://github.com/NNPDF/external/pull/58

## How to reproduce NNPDF40 c-factors

The scripts to compute the 4.0 c-factors are stored in the `cfactors_nnpdf40` folder.
These will use system calls to call `Vrap` with various arguments. If you compiled `vrap` as above the following should work (in the `cfactors_nnpdf40` folder):

```bash
export PATH=$PATH:$PWD/../build:
./runE605 1 1 # This will run NNLO, to run NLO just d `./runE605
```

Once we've run NNLO and NLO we can create the cfactor with the `runCfactors` script, commenting out the necessary datasets.

## Compare cfactors from NNPDF and vrap

There's a small python script that compares side by side the NNPDF and vrap cfactors:


```bash
python -i compare_cfactor.py cfactors_nnpdf40/output/E605_Cfactors.dat cfactors_nnpdf40/nnpdf40_theory200/CF_QCD_DYE605.dat
```

TODO: a script that compiles and runs and checks that the results are compatible with the NNPDF40 cfactors (i.e., regression tests)



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
