#pragma once
/*
 * Settings for Vrap, taken from `Vrap.C`
 */

#include "Vlumifns.h"
#include "pdf.h"
#include "pineappl_interface.h"
#include <string>

double E_CM;           // hadron-hadron center-of-mass energy
double Q, muR, muF;    // DY-mass; renormalization & factorization scales
double muRrel, muFrel; // global muR/Q, muF/Q (needed to do M (Q) integral)
double alphat;         // local value of alpha_QED(Q)
exchange
    exchM; // value of "exchange" for doing M (Q) integral, from `Vlumifns.h`
double Nf; // Nf = number of light quarks flavors; should ~ depend on Q.

// rapidity & integration ranges:
double y, y_lower, y_upper, xi_l, xi_u, ymax, Ml, Mu;

// pienappl interface
pinerap::CheffPanopoulos piner = pinerap::CheffPanopoulos();

// collider type (pp, ppbar, piso), from `pdf.h`
collider coll;

// technical settings
int ranseed, n_points;
bool NNLO_only;
int f_NNLO_only; // fudge parameter to isolate NNLO terms
                 // normally should be 1
                 // Set to 0 for only the NNLO hard cross section terms

int f_quiet = 0; // Set to 1 to suppress intermediate printing by rap_y()
                 // and by setV().

std::string pdfMode;
std::string pdfFile;

bool useMyAlphaRunning = false;
bool useOtherPDF;
bool jacobianTau2M = false;
int pdfSet;
int parton_flag; // 1: all, 2:qqbar 3: qg 4: gg 5:qq 6: qqbar_plus_qq

int o_f;
int direction;
int nbrYPnts;
