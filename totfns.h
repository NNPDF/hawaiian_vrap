/*==============================================================
 Functions for the DY total cross section; 
 formulae from Hamberg, van Neerven, & Matsuura
==============================================================*/
#ifndef TOTFNS_H
#define TOTFNS_H
// "plus" distributions (without the subtraction):
double DD1(double n, double y);
// Functions appearing in extra term required by the change of 
// variables to x1h and z1h.  To be evaluated at z = x2:
double DD1e(double n, double y);
// the NLO delta(1-z) terms:
double SV_Born_NLO(double muFQ);
// the NLO plus distribution terms:
double SV_boost_NLO(double z, double muFQ);
// Extra term required by the change of variables to x1h and z1h.
// To be evaluated at z = x2:
double SV_boost_NLOe(double z, double muFQ);
// the rest of the NLO q-\bar{q} terms:
double qbarq_hard_NLO(double z, double muFQ);
// NLO qg terms:
double qg_NLO(double z, double muFQ);
// NNLO functions ================================
// the NNLO delta(1-z) terms:
double SV_Born_NNLO(double NF, double muFQ, double muRQ);
// the NNLO plus distribution terms:
double SV_boost_NNLO(double z, double NF, double muFQ, double muRQ);
// Extra term required by the change of variables to x1h and z1h.
// To be evaluated at z = x2:
double SV_boost_NNLOe(double z, double NF, double muFQ, double muRQ);
// SOFT TEST VERSIONS of the NNLO plus distribution terms:
double SV_boost_NNLO_TEST(double z, double NF, double muFQ, double muRQ);
double SV_boost_NNLOe_TEST(double z, double NF, double muFQ, double muRQ);
// the q-\bar{q} "NFf" terms:
double qbarq_NFf(double Z);
// the q-\bar{q} "NF" terms:
double qbarq_NF(double Z, double muFQ, double muRQ);
// TEST function:  Integral of just the "hard" part of the NF terms.
double diffNF(double z);
// the q-\bar{q} C_A terms:
double qbarq_C_A(double Z, double muFQ, double muRQ);
// prior to Harlander-Kilgore correction:
double qbarq_C_A_preHK(double Z, double muFQ, double muRQ);
// the q-\bar{q} C_F terms:
double qbarq_C_F(double Z, double muFQ);
// prior to Harlander-Kilgore correction:
double qbarq_C_F_preHK(double Z, double muFQ);
// q-\bar{q} AC interference terms:
double qbarq_AC(double Z, double muFQ);
// q-\bar{q} BC interference terms:
double qbarq_BC(double Z);
// q-\bar{q} (or qq) C^2 terms:
// This is "CC_tot" in HvNM's NPB paper
// (or CC_tot_F, in FORTRAN program) -- see HvNMresults.
double qbarq_CC(double Z, double muFQ);
// q-\bar{q} (or qq) CD (vector) interference terms.
// This is "CD_V_tot" in HvNM's NPB paper
// (or -2*CD_V_tot_F, in FORTRAN program) -- see HvNMresults.
double qbarq_CD_V(double Z);
// The total q-\bar{q} "hard" cross section:
double qbarq_hard_NNLO(double z, double NF, double muFQ, double muRQ);
// In HvNM's notation, just the nonsinglet (NS) hard terms:
double qbarq_NS_hard_NNLO(double z, double NF, double muFQ, double muRQ);
// prior to Harlander-Kilgore correction:
double qbarq_NS_hard_NNLO_preHK(double z, double NF, double muFQ, double muRQ);
// TEST function:  Integral of just the "hard" part of the NO-NF terms.
double diffnoNF(double z);
// TEST function:  Integral of the "lonely" 1/(1-z)_+ parts of the NF terms.
double lonely(double z, double NF);
// to be evaluated at z = x2:
double lonely_e(double z, double NF);
//=========== quark-gluon terms ============================
double qg_C_A(double Z, double muFQ, double muRQ);
// prior to Harlander-Kilgore correction:
double qg_C_A_preHK(double Z, double muFQ, double muRQ);
double qg_C_F(double Z, double muFQ, double muRQ);
// prior to Harlander-Kilgore correction:
double qg_C_F_preHK(double Z, double muFQ, double muRQ);
double qg_NF(double Z, double muFQ, double muRQ); // just for muR <> muF
//=========== gluon-gluon terms ============================
double gg_C_A(double Z); 
double gg_C_F(double Z, double muFQ);
//=========== quark-quark terms ============================
// Interference exchange terms, for identical final state quarks.
// This is qq_CE_tot = qq_CE_tot_F(z): 
double qq_CE_tot(double Z, double muFQ);
// This is qq_CF_tot = 2 * qq_CF_tot_F(z): 
double qq_CF(double Z);
// Axial vector contributions:
// For AB, we insert a T_F with respect to the HvNM NPB formula !!!
// This is AB_ax_tot = AB_ax_tot_F(z):
double AB_ax(double Z);
// This is CD_ax_tot = 2 * CD_ax_tot_F(z):
double CD_ax(double Z);
# endif    /*  TOTFNS_H  */
