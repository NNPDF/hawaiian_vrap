/* =============================================================
   The functions required for the qg initial-state 
   NNLO contributions to the DY rapidity distributions.
   ============================================================= */

#ifndef QGFNS_H
#define QGFNS_H

// main "hard" function, and limits:
double qg(double y, double z);
double qg_u1(double z);              // u -> 1
double qg_uz(double y, double z);    // u -> z
double qg_u1mz(double y, double z);    // u -> 1 - z
double qg_z1(double y, double z);    // z -> 1
double qg_full(double y, double z);  // fully patched result
// The mu-dependent qg hard term:
double qg_mu(double y, double z, double muFQ); 
double qg_muR_hard(double y, double z, double Nf, double muFQ, double muRQ);
// "Boost" soft limits ("tot" now includes mu-dep. terms):
double qg_boost_mu(double z, double muFQ);  // the mu-dependent boost term
double qg_boost_muR(double z, double Nf, double muFQ, double muRQ); // muR
double qg_boost_tot(double z, double Nf, double muFQ, double muRQ);
// "Real" soft limits:
double qg_mu_aux(double y, double z, double muFQ);
double qg_muR_soft(double y, double z, double Nf, double muFQ, double muRQ);
// now includes mu- and muR-dep
double qg_soft_aux(double y, double z, double muFQ); 
double qg_real_soft(double y, double z, double Nf, double muFQ, double muRQ);
// remainder, treated as hard :
double qg_soft_to_hard(double y,double z, double muFQ);
//====================================================================
// The full q-\bar{q} real "hard" term, including conversion factor,
// and adding in all other terms multiplying "qg_lumi_ys":
double qg_real_hard(double y, double z, double Nf, double muFQ, double muRQ);
// EXPANDED RESULTS:
// expanded version of qg_full, from Frank:
double qg_full_EXP(double y, double z);
double qg_real_hard_EXP(double y, double z, double Nf, double muFQ, double muRQ);
# endif    /*  QGFNS_H  */
