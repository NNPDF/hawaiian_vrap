/* =============================================================
   pdf luminosity functions required for Z-gamma^* total cross section
   and rapidity distributions.
   ============================================================= */
#ifndef ZGLUMIFNS_H
#define ZGLUMIFNS_H
#include "pdf.h"

enum exchange { gamma_only = 1, Zgamma_interf = 2, Z_only = 3, Zgamma = 4,
                Wplus = 5, Wminus = 6 } ;
// procedures to decide which piece(s) are being computed:
void compute_all();
void compute_qqbar();
void compute_qg();
void compute_gg(); 
void compute_qq();
void compute_qqbar_plus_qq();
// procedure to decide which of Z/gamma^* exchanges are being computed:
void setV(exchange E, double Q, double alpha, double Nf, int f_quiet);
// the gamma-Z interference, pure-Z, and pure-W normalization factors:
double N_Zgamma(double Q, double alphat);
double N_Z(double Q, double alphat);
double N_W(double Q, double alphat);
// Computation of coupling parameters:
double VsqAsq_u();
double VsqAsq_d();
double VuVu();
double VdVd();
double VuVd();
double AAf_u(int Nf);
double AAf_d(int Nf);
double AuAu();
double AdAd();
double AuAd();
//======== Luminosity functions ============================
// Nf = number of light flavors; should depend on DY mass Q.
double NFf(int Nf);  // NFf = \sum_{i=1}^Nf e_f^2, generalized for (Z,gamma)
double Sum_a_f(int Nf);  // Sum_a_f = \sum_{f=1}^{Nf} a_f_Z.
// Sum_W_f =  F_W * ( v_W^2 + a_W^2) * \sum_{f,f'=1}^{Nf} |V_{ff'}|^2 :
double Sum_W_f(int Nf);
// QWf = F_W * ( v_W^2 + a_W^2) * \sum_{f=1}^{Nf} |V_{Qf}|^2 :
double UWf(int Nf);
double CWf(int Nf);
double DWf(int Nf);
double SWf(int Nf);
double BWf(int Nf);
double qqbar_lumi(double x1, double x2, double muF, process p, collider c);
double qqbar_BC_lumi(double x1, double x2, double muF, collider c);
double qqbar_ax_lumi(double x1, double x2, double muF, collider c);
double qg_lumi(double x1, double x2, double muF, collider c);
double gq_lumi(double x1, double x2, double muF, collider c);
double gg_lumi(double x1, double x2, double muF, collider c);
double qq_11_lumi(double x1, double x2, double muF, collider c);
double qq_22_lumi(double x1, double x2, double muF, collider c);
double qq_12_lumi(double x1, double x2, double muF, collider c);
double qq_12_ax_lumi(double x1, double x2, double muF, collider c);
double qq_CE1_lumi(double x1, double x2, double muF, collider c);
double qq_CE2_lumi(double x1, double x2, double muF, collider c);
double qq_CF_lumi(double x1, double x2, double muF, collider c);
# endif    /*  ZGLUMIFNS_H  */
