/* =============================================================
   The functions required for the qi \neq qj quark-(anti)quark initial-state 
   NNLO contributions to the DY rapidity distributions.
   eq1 and eq2 are the electric charges of the two incoming quarks.
   The "11" term is proportional to eq1^2 (C^2 of HvNM)
   The "12" term is proportional to eq1*eq2 (C*D of HvNM).
   The "22" term is proportional to eq2^2 (D^2 of HvNM),
   but is obtained from "11" using y -> 1-y.

   The file also includes the two identical-quark exchange contributions,
   CE and CF,
   and the two "double fermion loop" contributions for which the 
   axial vector coupling must be treated differently, AB and CD.
   ============================================================= */

#ifndef QIQJFNS_H
#define QIQJFNS_H

// "11" terms:
// main "hard" functions, and limits, for patching:
double qq11(double y, double z);
double qq11_u1(double z);             // u -> 1
double qq11_uz(double y, double z);   // u -> z
double qq11_u1mz(double y, double z); // u -> 1 - z
double qq11_z1(double y, double z);   // z -> 1
double qq11_full(double y, double z); // fully patched result
// The mu-dependent qq11 hard term:
double qq11_mu(double y, double z, double muFQ); 
// "11" "boost" terms --- delta(y) ("tot" now includes mu-dep. terms):
double qq11_boost_mu(double z, double muFQ);  // the mu-dependent boost term
double qq11_boost_tot(double z, double muFQ);   
// "11" "real soft" terms:
double qq11_mu_aux(double y, double z, double muFQ);
double qq11_soft_aux(double y, double z, double muFQ); // now includes mu-dep
double qq11_real_soft(double y, double z, double muFQ); 
// remainder, treated as hard:
double qq11_soft_to_hard(double y, double z, double muFQ);
//====================================================================
// The full qq11 real "hard" term, including conversion factor,
// and adding in all other terms multiplying "qq_11_lumi_ys":
// includes conversion factor:
double qq11_real_hard(double y, double z, double muFQ);
//=============================================
// "12" terms(there are no "soft" terms here):
// Main "hard" functions, and limits, for patching:
double qq12_full_1(double y, double z);  // one term, using symmetry
double qq12_u1(double z);             // u -> 1
double qq12_uz(double y, double z);   // u -> z
double qq12_z1(double y, double z);   // z -> 1
double qq12_full(double y, double z); // fully patched result
// includes conversion factor:
double qq12_real_hard(double y, double z);
//====================================================================
//====================================================================
// quark-quark exchange "CE" terms (mu-dependence only in boost term): 
double qq_CE(double y, double z);
double qq_CE_u1(double z);                    // u -> 1
double qq_CE_uz(double y, double z);          // u -> z
double qq_CE_z1(double y, double z);          // z -> 1
double qq_CE_full(double y, double z);
// "boost" soft limits
double qq_CE_boost_mu(double z, double muFQ);  // mu dep. term
double qq_CE_boost_tot(double z, double muFQ); // total now includes mu-dep.
double qq_CE_soft_aux(double y, double z);   // "real" soft limits
double qq_CE_real_soft(double y, double z);
double qq_CE_soft_to_hard(double y, double z);  // treated as hard
//====================================================================
// The full qq_CE real "hard" term, including conversion factor,
// and adding in all other "CE" terms multiplying "qq_ident_lumi_ys":
double qq_CE_real_hard(double y, double z);
//====================================================================
// quark-quark exchange "CF" terms. 
double qq_CF_full_1(double y, double z);
double qq_CF_u1(double z);            // u -> 1
double qq_CF_z1(double y, double z);  // z -> 1
double qq_CF_full(double y, double z);
double qq_CF_real_hard(double y, double z);
//====================================================================
// quark-(anti)-quark axial vector "AB" terms. 
double qq_AB_ax_full_1(double y, double z);
double qq_AB_ax_u1(double z);  // u -> 1 
double qq_AB_ax_z1(double y, double z);  // z -> 1 
double qq_AB_ax_z0(double y, double z);  // z -> 0 
double qq_AB_ax_full(double y, double z);
double qq_AB_ax_real_hard(double y, double z);
//====================================================================
// quark-(anti)-quark axial vector "CD" terms. 
double qq_CD_ax_full_1(double y, double z);
double qq_CD_ax_z1(double y, double z);  // z -> 1 
double qq_CD_ax_u1(double z);  // u -> 1 
double qq_CD_ax_full(double y, double z);
double qq_CD_ax_real_hard(double y, double z);
// expanded version of qq11_full and qq12_full, from Frank:
double qq11_full_EXP(double y, double z);
double qq11_real_hard_EXP(double y, double z, double muFQ);
double qq12_full_EXP(double y, double z);
double qq12_real_hard_EXP(double y, double z);
double qq_AB_ax_full_EXP(double y, double z);
double qq_AB_ax_real_hard_EXP(double y, double z);
double qq_CD_ax_full_EXP(double y, double z);
double qq_CD_ax_real_hard_EXP(double y, double z);
# endif    /*  QIQJFNS_H  */
