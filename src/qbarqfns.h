/* =============================================================
   The functions required for the q-\bar{q} initial-state 
   NNLO contributions to the DY rapidity distributions.
   ============================================================= */

#ifndef QBARQFNS_H
#define QBARQFNS_H

// no-NF function, and limits:
double qbarq_no_NF_1(double y, double z);
double qbarq_no_NF_u1(double z);  // u -> 1
double qbarq_no_NF_uz(double y, double z);  // u -> z
double qbarq_no_NF_z1(double y, double z);  // z -> 1
double qbarq_no_NF(double y, double z);   // fully patched result
// NF function, and limits:
double qbarq_NF_1(double y, double z);
double qbarq_NF_z1(double y, double z);  // z -> 1
double qbarq_NF(double y, double z);   // fully patched result
// "NFf" function, and limits:
double qbarq_NFf_1(double y, double z);
double qbarq_NFf_z1(double y, double z);  // z -> 1
double qbarq_NFf(double y, double z);   // fully patched result
// terms proportional to NFf = \sum_{i=1}^Nf e_f^2 (or generalization)
double qbarq_NFf_real_hard(double y, double z);
// "Boost" soft limits:
double qbarq_boost_mu_soft(double z, double NF, double muFQ); 
double qbarq_boost_muR_soft(double z, double Nf, double muFQ, double muRQ);
double qbarq_boost_soft(double z, double NF, double muFQ, double muRQ); 
double qbarq_boost_mu_hard(double z, double NF, double muFQ); 
double qbarq_boost_muR_hard(double z, double Nf, double muFQ, double muRQ);
double qbarq_boost_hard(double z, double NF, double muFQ, double muRQ); 
// "Real" soft limits (muF-dep. terms now included here):
double qbarq_real_soft_A(double y, double z, double NF, double muFQ, double muRQ);
double qbarq_plus2_A(double y, double z, double Nf, double muFQ);
double qbarq_real_soft_extra_A(double y, double z, double NF, double muFQ);
double qbarq_soft_to_hard_A(double y, double z, double NF, double muFQ);
double qbarq_real_soft_z1_A(double y, double z, double NF, double muFQ);
double qbarq_real_soft(double y, double z, double NF, double muFQ, double muRQ);
double qbarq_plus2(double y, double z, double Nf, double muFQ);
double qbarq_real_soft_extra(double y, double z, double NF, double muFQ);
double qbarq_soft_to_hard(double y, double z, double NF, double muFQ);
double qbarq_real_soft_0(double y, double z, double NF, double muFQ, double muRQ); 
double qbarq_real_soft_1(double y, double z, double NF, double muFQ, double muRQ);
double qbarq_real_soft_z1(double y, double z, double NF, double muFQ);
// Just the NF coefficients of qbarq_real_soft terms (muF=Q only):
double qbarq_NF_real_soft_A(double y, double z);
double qbarq_NF_real_soft_extra_A(double y, double z);
double qbarq_NF_soft_to_hard_A(double y, double z);
double qbarq_NF_real_soft_z1_A(double y, double z);
double qbarq_NF_real_soft(double y, double z);
double qbarq_NF_real_soft_extra(double y, double z);
double qbarq_NF_soft_to_hard(double y, double z);
double qbarq_NF_real_soft_0(double y, double z);
double qbarq_NF_real_soft_1(double y, double z);
double qbarq_NF_real_soft_z1(double y, double z);
//====================================================================
// For W production, we have to separate out the "BC" interference terms.
// They are only hard terms, have no factorization scale dependence,
// and have a y -> 1-y symmetry.
// BC function, and limits:
double qbarq_BC_1(double y, double z);  // before adding y -> 1-y
double qbarq_BC_u1(double z);  // u -> 1
double qbarq_BC_z1(double y, double z);  // z -> 1
double qbarq_BC(double y, double z);   // fully patched result
double qbarq_BC_real_hard(double y, double z);  // includes phase-space
//====================================================================
// The full q-\bar{q} real "hard" term, including conversion factor,
// and adding in all other terms multiplying "lumi_ys",
// as well as muF-dependent terms:
// just the mu-dependent hard terms
double qbarq_mu_1(double y, double z, double NF, double muFQ);
double qbarq_mu(double y, double z, double NF, double muFQ);
double qbarq_muR_hard(double y, double z, double Nf, double muFQ, double muRQ);
// grand total hard term:
double qbarq_real_hard(double y, double z, double NF, double muFQ, double muRQ);
// expanded version of qbarq_real_hard, from Frank:
double qbarq_no_NF_NS_EXP(double y, double z);
double qbarq_NF_EXP(double y, double z);
double qbarq_NFf_EXP(double y, double z);
double qbarq_NFf_real_hard_EXP(double y, double z);
double qbarq_BC_EXP(double y, double z);
double qbarq_BC_real_hard_EXP(double y, double z);
double qbarq_real_hard_EXP(double y, double z, double NF, double muFQ, double muRQ);
# endif    /*  QBARQFNS_H  */
