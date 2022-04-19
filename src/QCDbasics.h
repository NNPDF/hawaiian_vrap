/* Basic constants and functions for QCD calculations.
  Inherited from "kinematics/kinematics.C(h)"
  but without maple/kinematics stuff.
  Also, uses different way of flagging order in alpha_s (order_flag), 
  and alpha_s incorporates 3-loop case.  Includes L_i and Ls_i functions */
#ifndef QCDBASICS_H
#define QCDBASICS_H

#include <iomanip>
#include "complex.h"
#include "thbasics.h"
#include "random.h"

const double GeV2fb = 389379292000.;   // GeV^2*fb conversion constant
const double GeV2pb = 389379292.;   // GeV^2*pb conversion constant
// QED coupling, and Fermi constant
const double alpha0 = 1/137.036;    const double G_F = 1.16639e-05;
// sin^2 theta_W(M_Z)_{MS-bar}, "effective" (Z pole asymmetry) value
const double s2W = 0.23143;   const double sW = sqrt(s2W);  
const double c2W = 1.-s2W;   const double cW = sqrt(c2W);
// QCD color constants
const double N_c = 3.;
const double C_F = (N_c*N_c-1.)/2./N_c;
const double C_A = N_c;
const double T_R = 0.5;
const double N_f = 5.;        // for u+c + d+s+b --- watch out, this could change!
const double sumQsq = 2. * 4./9. + 3. * 1./9.;  // (charge)^2 sum for:  u+c + d+s+b 
const double sumQ4 = 2. * 16./81. + 3. * 1./81.;  
// (charge)^4 sum for:  u+c + d+s+b   (don't forget leptons too, for lbyl)
const double sumQ6 = 2. * 64./729. + 3. * 1./729.;
// (charge)^4 sum for:  u+c + d+s+b   (don't forget leptons too, for lbyl)
// individual quark charges and weak isospins:
const double e_u = 2./3.;  const double e_c = 2./3.;  const double e_t = 2./3.; 
const double e_d = -1./3.;  const double e_s = -1./3.;  const double e_b = -1./3.; 
const double I_3u = 0.5;  const double I_3c = 0.5;  const double I_3t = 0.5; 
const double I_3d = -0.5;  const double I_3s = -0.5;  const double I_3b = -0.5; 
// lepton masses
const double m_e = 0.0005109989;  const double m_mu = 0.105658;  const double m_tau = 1.77703;
// "best (?) light quark masses" 
const double m_u = 0.005;  const double m_d = 0.005;  const double m_s = 0.100;   
// current set of heavy quark masses, m_q(m_q):
const double m_c = 1.25; const double m_b = 4.24;  
// value of m_t(m_t) corresponding to  m_t(pole) = 174.3,
// using Melnikov-van Ritbergen 3-loop relation, for alpha_s(M_Z) = 0.119
const double m_t = 164.6;

// W and Z masses:
const double m_W = 80.398;   
const double m_Z = 91.1876;
// running couplings:
const double alpha_QED_Z = 1/128.9;   // alpha_QED(m_Z).
double alpha_QED_low(double mu); // 1-loop running of alpha_QED, valid
                             // between 0 and (about) 2*m_mu.
double alpha_QED_OLD(double mu); // 1-loop running of alpha_QED, valid
                             // between 2*m_b and 2*m_W.
double alpha_QED(double mu); // 1-loop running of alpha_QED, should be
                             // ~ valid between 0 and 2*m_W.

extern int order_flag;  // 0 for LO (as1); 1 for NLO (as2), 2 for NNLO (as3). 
                 // (Make sure to use consistent pdf's!)
                 // Also controls which alpha_s and M_q evolution is used. 

extern double alpha_s_Z;   // alpha_s(m_Z).  Should depend on "order_flag".
double alpha_s(double mu); // Up to 3-loop PDG formula for alpha_s(mu)
                           // alpha_s(m_Z) = alpha_s_Z is adjustable.
// Calls these functions, alpha_s^{1,2,3-loop}(mu)/Pi.
double as1(double mu);
double as2(double mu);
double as3(double mu);
double iterlnx3(double lnx, double alpha_s, double N_f);  // used by as3
double as3_fixed(double mu); // Three-loop PDG formula for alpha_s(mu)/PI
                             // with alpha_s(m_Z) = 0.118 fixed.
// Formulae to run c and b quark MS-bar masses
double cq(double x, int Nf);
double M_q(double m_q, double mu1, double mu2, int Nf);
double M_t(double mu);
double M_b(double mu);
double M_c(double mu);
//
//inline double ran(int n=0);
int signum(double a);
complex Clog(double s1, double s2); // complex log of a ratio
complex Clog1(double s);   // complex log of a single s
complex L0(double s1, double s2);
complex L0(double s1, double s2);
complex L1(double s1, double s2);
complex L2(double s1, double s2);
// Li_2(1-r), for r = s1/s2:
complex CLi2r(double s1, double s2);
complex Lsm1(double s1, double s2, double s3, double s4);
complex Ls0(double s1, double s2, double s3, double s4);
complex Ls1(double s1, double s2, double s3, double s4);
complex Ls2(double s1, double s2, double s3, double s4);
complex Ls3(double s1, double s2, double s3, double s4);

inline double ran(int n=0) {return ran2(n);}

# endif    /*  QCDBASICS_H  */
