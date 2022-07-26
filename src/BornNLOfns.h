/* =============================================================
   The functions required for the Born kinematics and
   NLO contributions to the DY rapidity distributions.
   ============================================================= */

#ifndef BORNNLOFNS_H
#define BORNNLOFNS_H

// The NLO & NNLO Born kinematics terms:
double Born_NLO(double muFQ);
double Born_NNLO(double Nf, double muFQ, double muRF, double *logterms);
// The NLO boost terms, and "soft" (z=1) subtraction (for q-\bar{q]):
double NLO_qbarq_boost(double z, double muFQ);
double NLO_qbarq_boost_soft(double z, double muFQ);
double NLO_qg_boost(double z, double muFQ);
// The NLO real terms, and y=0, y=1, subtractions:
double NLO_qbarq_real(double y, double z);
double NLO_qbarq_real_soft(double y, double z);
double NLO_qg_real(double y, double z);
double NLO_qg_real_soft(double y, double z);
# endif    /*  BORNNLOFNS_H  */
