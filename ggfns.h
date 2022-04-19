/* =============================================================
   The functions required for the gg initial-state 
   NNLO contributions to the DY rapidity distributions.
   ============================================================= */

#ifndef GGFNS_H
#define GGFNS_H

// main "hard" function, and limits, for patching:
double gg_full_1(double y, double z);
double gg_u1(double z);             // u -> 1
double gg_uz(double y, double z);   // u -> z
double gg_z1(double y, double z);   // z -> 1
double gg_full(double y, double z); // fully patched result
// The mu-dependent gg terms:
double gg_mu_1(double y, double z, double muFQ); 
double gg_mu(double y, double z, double muFQ);
// includes conversion factor and mu-dependence:
double gg_real_hard(double y, double z, double muFQ);
// EXPANDED RESULTS:
// expanded version of gg_full, from Frank:
double gg_full_EXP(double y, double z);
double gg_real_hard_EXP(double y, double z, double muFQ);
# endif    /*  GGFNS_H  */
