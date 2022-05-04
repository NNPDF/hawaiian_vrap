/* =============================================================
   pdf luminosity functions required for Z-gamma^* total cross section
   and rapidity distributions.
   ============================================================= */
#include <iomanip>
#include "NAClasses.h"
#include "pdf.h"  // pdf distributions
#include "Vlumifns.h"
#include "LHApdf.h"
#include "EWparams.h"

using namespace std;

// defined in Vlimufns.C
extern int f_qqbar, f_qg, f_gg, f_qq;

// Parameters to weight gamma^* vs. Z vs. W exchange
// (all terms, LO, NLO, NNLO!):
extern double F_gamma, F_Zgamma, F_Z, F_W;

// flag W+, W-:
extern int fWp ;   extern int fWm ;

// storage of coupling parameters:
extern double vsqasq_u, vsqasq_d, NF_f, vuvu, vdvd, vuvd; 
extern double auau, adad, auad, aaf_u, aaf_d;
// these contain CKM angles (for W production):
extern double sum_W_f, uWf, cWf, dWf, sWf, bWf;



//======== Luminosity functions ============================

/* The (v^2+a^2)-weighted parton luminosity function 
   that always appears in main q-qbar PLUS qbar-q 
   contribution to Z-gamma^* production.
   For archaic reasons we use a "process" flag to weight the "charges"
   "Drell-Yan" (for process = 1) -> (v^2+a^2)
   or "gluon" (for "NFf" term) (for process = 0) -> no-charge weighting here.
   Collider can be pp (for collider = 1) or p-pbar (for collider = 2). */

double qqbar_lumi(const pdfArray& X1, const pdfArray& X2, process p, collider c){
  if ((f_qqbar == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  double pref, pref_u, pref_d;
// The "DY" type W^\pm cases are cumbersome because of CKM elements:
  if ( ((fWp==1) || (fWm==1)) && (p==DY) ){
    pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWp==1) && (c==pp)){ return pref * ( 
       X1.u * ( V_ud_sq * X2.dbar 
	           + V_us_sq * X2.s + V_ub_sq * X2.b )
     + X1.c * ( V_cd_sq * X2.dbar 
	           + V_cs_sq * X2.s + V_cb_sq * X2.b )
     + X2.u * ( V_ud_sq * X1.dbar 
	           + V_us_sq * X1.s + V_ub_sq * X1.b )
     + X2.c * ( V_cd_sq * X1.dbar 
		     + V_cs_sq * X1.s + V_cb_sq * X1.b ) );
    }
//
    if ((fWp==1) && (c==ppbar)){ return pref * ( 
       X1.u * ( V_ud_sq * X2.d 
	           + V_us_sq * X2.s + V_ub_sq * X2.b )
     + X1.c * ( V_cd_sq * X2.d 
	           + V_cs_sq * X2.s + V_cb_sq * X2.b )
     + X2.ubar * ( V_ud_sq * X1.dbar 
	           + V_us_sq * X1.s + V_ub_sq * X1.b )
     + X2.c * ( V_cd_sq * X1.dbar 
		     + V_cs_sq * X1.s + V_cb_sq * X1.b ) );
    }
//
    if ((fWm==1) && (c==pp)){ return pref * ( 
       X1.ubar * ( V_ud_sq * X2.d 
	           + V_us_sq * X2.s + V_ub_sq * X2.b )
     + X1.c * ( V_cd_sq * X2.d 
	           + V_cs_sq * X2.s + V_cb_sq * X2.b )
     + X2.ubar * ( V_ud_sq * X1.d 
	           + V_us_sq * X1.s + V_ub_sq * X1.b )
     + X2.c * ( V_cd_sq * X1.d 
		     + V_cs_sq * X1.s + V_cb_sq * X1.b ) );
    }
//
    if ((fWm==1) && (c==ppbar)){ return pref * ( 
       X1.ubar * ( V_ud_sq * X2.dbar 
	           + V_us_sq * X2.s + V_ub_sq * X2.b )
     + X1.c * ( V_cd_sq * X2.dbar 
	           + V_cs_sq * X2.s + V_cb_sq * X2.b )
     + X2.u * ( V_ud_sq * X1.d 
	           + V_us_sq * X1.s + V_ub_sq * X1.b )
     + X2.c * ( V_cd_sq * X1.d 
		     + V_cs_sq * X1.s + V_cb_sq * X1.b ) );
    }
  }
// The "gluon" W^\pm and all the (Z,gamma) cases are the same up to 
// overall prefactor:
  if (p==gluon){ 
    if ((fWp==1) || (fWm==1)){ 
      pref_u = sum_W_f;   
      pref_d = pref_u;
    }
    else { pref_u = NF_f;    pref_d = pref_u; }
  }
  else if (p==DY){ pref_u = vsqasq_u;  pref_d = vsqasq_d; }
  else return 0.;
//  
  if (c==pp) { return 
     pref_u * ( X1.u * X2.ubar + X1.ubar * X2.u 
		                + 2. * X1.c * X2.c )
   + pref_d * ( X1.d * X2.dbar + X1.dbar * X2.d 
	                        + X1.s * X2.sbar + X1.sbar * X2.s 
		                + 2. * X1.b * X2.b )
 ; }
  else if (c==ppbar) { return 
     pref_u * ( X1.u * X2.u + X1.ubar * X2.ubar 
                                + 2. * X1.c * X2.c )
   + pref_d * ( X1.d * X2.d + X1.dbar * X2.dbar 
	                        + X1.s * X2.s + X1.sbar * X2.sbar 
                                + 2. * X1.b * X2.b ) ; }
  else if (c==piso) { return 
      pref_u * ( X1.u * (X2.ubar+X2.dbar)/2. + X1.ubar * (X2.u+X2.d)/2.
                                + 2. * X1.c * X2.c )
      + pref_d * ( X1.d * (X2.dbar+X2.ubar)/2. + X1.dbar * (X2.d+X2.u)/2.
	                        + X1.s * X2.sbar + X1.sbar * X2.s 
                                + 2. * X1.b * X2.b ) ; }
  else return 0.;
}
double qqbar_lumi_dy(const pdfArray& X1, const pdfArray& X2, collider c){
    return qqbar_lumi(X1, X2, DY, c);
}

/* For W^\pm production, the "BC" terms get a different type of luminosity
   function.  This function is equal to 2 * qbarq_lumi(DY) for the (Z,gamma)
   case. */

double qqbar_BC_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qqbar == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
// The W^+ and W^- cases are the same here:
  if ((fWp==1) || (fWm==1)) { 
    if (c==pp) { return
        uWf * ( X1.u * X2.ubar + X1.ubar * X2.u ) 
      + 2. * cWf * X1.c * X2.c
      + dWf * ( X1.d * X2.dbar + X1.dbar * X2.d )
      + 2. * sWf * X1.s * X2.s 
      + 2. * bWf * X1.b * X2.b ;
    }
    if (c==ppbar) { return
        uWf * ( X1.u * X2.u + X1.ubar * X2.ubar ) 
      + 2. * cWf * X1.c * X2.c
      + dWf * ( X1.d * X2.d + X1.dbar * X2.dbar )
      + 2. * sWf * X1.s * X2.s 
      + 2. * bWf * X1.b * X2.b ;
    }
  }
  else if ((c==pp) || (c==ppbar) || c==piso) { return 2. * qqbar_lumi(X1,X2,DY,c); }
  return 0.;
}
    
// q-\bar{q} luminosity function required for the "AB axial" terms.
// It should return 0 for W^\pm production.

double qqbar_ax_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qqbar == 0) || (fWp==1) || (fWm==1) || (X1.x >= 1.) || (X2.x >= 1.)) { 
    return 0.; 
  }
  if (c==pp) { return 
     aaf_u * ( X1.u * X2.ubar + X1.ubar * X2.u 
                                + 2. * X1.c * X2.c )
   + aaf_d * ( X1.d * X2.dbar + X1.dbar * X2.d 
                                + X1.s * X2.sbar + X1.sbar*X2.s 
                                + 2. * X1.b * X2.b ) ; }
  else if (c==ppbar) { return 
     aaf_u * ( X1.u * X2.u + X1.ubar * X2.ubar 
                                + 2. * X1.c * X2.c )
   + aaf_d * ( X1.d * X2.d + X1.dbar * X2.dbar 
	                        + X1.s * X2.s + X1.sbar * X2.sbar
                                + 2. * X1.b * X2.b ) ; }
  else if (c==piso) { return 
      aaf_u * ( X1.u * (X2.ubar+X2.dbar)/2. + X1.ubar * (X2.u+X2.d)/2. 
                                + 2. * X1.c * X2.c )
      + aaf_d * ( X1.d * (X2.ubar+X2.dbar)/2. + X1.dbar * (X2.u+X2.d)/2.
                                + X1.s * X2.sbar + X1.sbar*X2.s 
                                + 2. * X1.b * X2.b ) ; }
  else return 0.;
}

/* The quark-charge weighted parton luminosity function 
   that always appears in the qg PLUS qbar-g (assuming C invariance!) 
   initiated contribution to W+, W-, or gamma^* / Z production
   for pp (for collider = 1) or p-pbar (for collider = 2). */ 

double qg_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qg == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if (fWp==1){ return 
      ( uWf * X1.u + cWf * X1.c
      + dWf * X1.dbar + sWf * X1.s + bWf * X1.b )
      * X2.gluon;
  }
  else if (fWm==1){ return 
      ( uWf * X1.ubar + cWf * X1.c
      + dWf * X1.d + sWf * X1.s + bWf * X1.b )
      * X2.gluon;
  }
  else if ((c==pp) || (c==ppbar) || c==piso) { return 
     ( vsqasq_u * (X1.u + X1.ubar + 2. * X1.c)
     + vsqasq_d * (X1.d + X1.dbar 
                  + X1.s + X1.sbar + 2. * X1.b) )
      * X2.gluon;
  }
  else return 0.;
}

/* And similarly for gq PLUS g-qbar (not quite so simple in W case;
   can't just use x1 <-> x2 because of p vs. pbar distinction). */

double gq_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qg == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if (((fWp==1) && (c==pp)) || ((fWm==1) && (c==ppbar))) { return 
      ( uWf * X2.u + cWf * X2.c
      + dWf * X2.dbar + sWf * X2.s + bWf * X2.b )
      * X1.gluon;
  }
  else if (((fWm==1) && (c==pp)) || ((fWp==1) && (c==ppbar))) { return 
      ( uWf * X2.ubar + cWf * X2.c
      + dWf * X2.d + sWf * X2.s + bWf * X2.b )
      * X1.gluon;
  }
  else if ((c==pp) || (c==ppbar)) { return 
     ( vsqasq_u * (X2.u + X2.ubar + 2. * X2.c)
     + vsqasq_d * (X2.d + X2.dbar 
                  + X2.s + X2.sbar + 2. * X2.b) )
      * X1.gluon;
  }
  else if (c==piso) { return 
      ( vsqasq_u * ((X2.u + X2.ubar + X2.d+X2.dbar)/2. + 2. * X2.c)
	+ vsqasq_d * ((X2.u + X2.ubar + X2.d+X2.dbar)/2. 
		      + X2.s + X2.sbar + 2. * X2.b) )
      * X1.gluon;
  }
  else return 0.;
}

/* The gluon-gluon parton luminosity function (same for pp or p-pbar).
 For (Z,gamma) production it always appears multiplied by "NFf";
 for W^\pm production it is multiplied by sum_W_f.  */

double gg_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_gg == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  double pref;
  if ((fWp==1) || (fWm==1)) { pref = sum_W_f; }
  else { pref = NF_f; }
  return pref * X1.gluon * X2.gluon ;
}

/* The (generalized) quark-charge weighted parton luminosity function 
   that appears in the "C^2" (v_1^2 + a_1^2) q_i q_j (qbar_j) PLUS charge-conj
   initiated contributions to Drell-Yan production.
   Here (v_1^2 +a_1^2) means both quark charges are associated with hadron 1.
   The function is the same for both pp (collider = 1) 
   and p-pbar (collider = 2). */

double qq_11_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if (fWp==1){ return 
      ( uWf * X1.u + cWf * X1.c
      + dWf * X1.dbar + sWf * X1.s + bWf * X1.b )
   * ( X2.u + X2.ubar + 2.*X2.c
     + X2.d + X2.dbar + 2.*X2.s + 2.*X2.b ); 
  }
  else if (fWm==1){ return 
      ( uWf * X1.ubar + cWf * X1.c
      + dWf * X1.d + sWf * X1.s + bWf * X1.b )
   * ( X2.u + X2.ubar + 2.*X2.c
     + X2.d + X2.dbar + 2.*X2.s + 2.*X2.b ); 
  }
  else if ((c==pp) || (c==ppbar) || c==piso){ return 
    ( vsqasq_u * ( X1.u + X1.ubar + 2.*X1.c )
    + vsqasq_d * ( X1.d + X1.dbar + X1.s + X1.sbar + 2.*X1.b ) )  
   * ( X2.u + X2.ubar + 2.*X2.c
     + X2.d + X2.dbar + X2.s + X2.sbar + 2.*X2.b );
  }
  else return 0.;
}

/* Same function, except for e_2^2, i.e. both quark charges are associated
   with hadron 2.  Get by letting (x1 <-> x2) in qq_11_lumi, except for
   the case of W^\pm production... */

double qq_22_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if (((fWp==1) && (c==pp)) || ((fWm==1) && (c==ppbar))) { return 
      ( uWf * X2.u + cWf * X2.c
      + dWf * X2.dbar + sWf * X2.s + bWf * X2.b )
   * ( X1.u + X1.ubar + 2.*X1.c
     + X1.d + X1.dbar + 2.*X1.s + 2.*X1.b ); 
  }
  else if (((fWm==1) && (c==pp)) || ((fWp==1) && (c==ppbar))) { return 
      ( uWf * X2.ubar + cWf * X2.c
      + dWf * X2.d + sWf * X2.s + bWf * X2.b )
   * ( X1.u + X1.ubar + 2.*X1.c
     + X1.d + X1.dbar + 2.*X1.s + 2.*X1.b ); 
  }
  else if ((c==pp) || (c==ppbar)){ return 
    ( vsqasq_u * ( X2.u + X2.ubar + 2.*X2.c )
    + vsqasq_d * ( X2.d + X2.dbar + X2.s + X2.sbar + 2.*X2.b ) )  
   * ( X1.u + X1.ubar + 2.*X1.c
     + X1.d + X1.dbar + X1.s + X1.sbar + 2.*X1.b );
  }
  else if (c==piso){ return 
      ( vsqasq_u * ( (X2.u + X2.ubar + X2.d + X2.dbar)/2. + 2.*X2.c )
    + vsqasq_d * ( (X2.u + X2.ubar + X2.d + X2.dbar)/2. + X2.s + X2.sbar + 2.*X2.b ) )  
   * ( X1.u + X1.ubar + 2.*X1.c
     + X1.d + X1.dbar + X1.s + X1.sbar + 2.*X1.b );
  }
  else return 0.;
}

/* Same function, except for v_1*v_2, 
   i.e. one quark vector charge is associated with each hadron. 
   Note minus signs under q <-> qbar (vector case).
   Hence the nonvalence quarks do not contribute.
   Overall (-) is so q-qbar enters with plus sign.
   Like qqbar_ax_lumi, this function should return 0 for W^\pm production. */

double qq_12_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (fWp==1) || (fWm==1) || (X1.x >= 1.) || (X2.x >= 1.)) { 
    return 0.;
  }
  double temp = - ( 
      vuvu * (X1.u - X1.ubar) * (X2.u - X2.ubar)
    + vuvd * ( (X1.u - X1.ubar) * (X2.d - X2.dbar)
	     + (X1.d - X1.dbar) * (X2.u - X2.ubar) )
    + vdvd * (X1.d - X1.dbar) * (X2.d - X2.dbar) ) ;
  if (c==pp){ return temp; }
  else if (c==ppbar) { return - temp; }
  else if(c==piso) {
    return - ( 
	      vuvu * (X1.u - X1.ubar) 
	      + vuvd * ( (X1.u - X1.ubar + X1.d - X1.dbar)
			 + vdvd * (X1.d - X1.dbar))
	      * (X2.u - X2.ubar + X2.d - X2.dbar)/2. );
  }
  else return 0.;
}


/* Same function as qq_12_lumi, except contains axial vector coupling
  product, a_1*a_2. It also returns 0 for W^\pm production. */

double qq_12_ax_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (fWp==1) || (fWm==1) || (X1.x >= 1.) || (X2.x >= 1.)) {
     return 0.;
  }
  if ((c==pp) || (c==ppbar)){ return 
      auau * (X1.u + X1.ubar + 2. * X1.c) 
           * (X2.u + X2.ubar + 2. * X2.c)
    + auad * ( (X1.u + X1.ubar + 2. * X1.c)
             * (X2.d + X2.dbar + X2.s + X2.sbar + 2. * X2.b)
             + (X1.d + X1.dbar + X1.s + X1.sbar + 2. * X1.b)
	     * (X2.u + X2.ubar + 2. * X2.c) )
    + adad * (X1.d + X1.dbar + X1.s +X1.sbar + 2. * X1.b)
           * (X2.d + X2.dbar + X2.s +X2.sbar + 2. * X2.b) ;
  }
  else if (c==piso){ return 
      auau * (X1.u + X1.ubar + 2. * X1.c) 
      * ((X2.u + X2.ubar +X2.d + X2.dbar)/2. + 2. * X2.c)
      + auad * ( (X1.u + X1.ubar + 2. * X1.c)
		 * ((X2.u + X2.ubar +X2.d + X2.dbar)/2. + X2.s + X2.sbar + 2. * X2.b)
		 + (X1.d + X1.dbar + X1.s + X1.sbar + 2. * X1.b)
		 * ((X2.u + X2.ubar +X2.d + X2.dbar)/2. + 2. * X2.c) )
      + adad * (X1.d + X1.dbar + X1.s +X1.sbar + 2. * X1.b)
      * ((X2.u + X2.ubar +X2.d + X2.dbar)/2. + X2.s +X2.sbar + 2. * X2.b) ;
  }
  else return 0.;
}

/* "CE1" and "CE2" luminosities for identical-quark (or identical antiquark)
   exchange contribution. 
   "CE2" is same as "CE1" in (Z,gamma) case; but not W^\pm case. */

double qq_CE1_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
  double pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWp==1) && (c==pp)) { return pref * (
      V_ud_sq * ( X1.u * X2.d + X1.dbar * X2.ubar )
    + V_us_sq * ( X1.u * X2.s + X1.s * X2.ubar )
    + V_ub_sq * ( X1.u * X2.b + X1.b * X2.ubar )
    + V_cd_sq * ( X1.c * X2.d + X1.dbar * X2.c )
    + V_cs_sq * ( X1.c * X2.s + X1.s * X2.c )
    + V_cb_sq * ( X1.c * X2.b + X1.b * X2.c ) );
    }
    if ((fWp==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( X1.u * X2.dbar + X1.dbar * X2.u )
    + V_us_sq * ( X1.u * X2.s + X1.s * X2.u )
    + V_ub_sq * ( X1.u * X2.b + X1.b * X2.u )
    + V_cd_sq * ( X1.c * X2.dbar + X1.dbar * X2.c )
    + V_cs_sq * ( X1.c * X2.s + X1.s * X2.c )
    + V_cb_sq * ( X1.c * X2.b + X1.b * X2.c ) );
    }
    if ((fWm==1) && (c==pp)) { return pref * (
      V_ud_sq * ( X1.d * X2.u + X1.ubar * X2.dbar )
    + V_us_sq * ( X1.s * X2.u + X1.ubar * X2.s )
    + V_ub_sq * ( X1.b * X2.u + X1.ubar * X2.b )
    + V_cd_sq * ( X1.d * X2.c + X1.c * X2.dbar )
    + V_cs_sq * ( X1.s * X2.c + X1.c * X2.s )
    + V_cb_sq * ( X1.b * X2.c + X1.c * X2.b ) );
    }
    if ((fWm==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( X1.d * X2.ubar + X1.ubar * X2.d )
    + V_us_sq * ( X1.s * X2.ubar + X1.ubar * X2.s )
    + V_ub_sq * ( X1.b * X2.ubar + X1.ubar * X2.b )
    + V_cd_sq * ( X1.d * X2.c + X1.c * X2.d )
    + V_cs_sq * ( X1.s * X2.c + X1.c * X2.s )
    + V_cb_sq * ( X1.b * X2.c + X1.c * X2.b ) );
    }
    return 0.;
  }
  else if (c==pp){ return
     vsqasq_u * ( X1.u * X2.u + X1.ubar * X2.ubar 
                + 2. * X1.c * X2.c )
   + vsqasq_d * ( X1.d * X2.d + X1.dbar * X2.dbar 
        + X1.s * X2.s + X1.sbar * X2.sbar + 2. * X1.b * X2.b ); }
  else if (c==ppbar) { return
     vsqasq_u * ( X1.u * X2.ubar + X1.ubar * X2.u 
               + 2. * X1.c * X2.c )
   + vsqasq_d * ( X1.d * X2.dbar + X1.dbar * X2.d 
        + 2. * X1.s * X2.s + 2. * X1.b * X2.b ); }
  else if (c==piso){ return
      vsqasq_u * ( X1.u * (X2.u+X2.d)/2. + X1.ubar * (X2.ubar+X2.dbar)/2.
		   + 2. * X1.c * X2.c )
      + vsqasq_d * ( X1.d * (X2.d+X2.u)/2. + X1.dbar * (X2.dbar+X2.ubar)/2.
		     + X1.s * X2.s + X1.sbar * X2.sbar + 2. * X1.b * X2.b ); }
  else return 0.;
}

double qq_CE2_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
  double pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWm==1) && (c==pp)) { return pref * (
      V_ud_sq * ( X1.u * X2.d + X1.dbar * X2.ubar )
    + V_us_sq * ( X1.u * X2.s + X1.s * X2.ubar )
    + V_ub_sq * ( X1.u * X2.b + X1.b * X2.ubar )
    + V_cd_sq * ( X1.c * X2.d + X1.dbar * X2.c )
    + V_cs_sq * ( X1.c * X2.s + X1.s * X2.c )
    + V_cb_sq * ( X1.c * X2.b + X1.b * X2.c ) );
    }
    else if ((fWm==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( X1.u * X2.dbar + X1.dbar * X2.u )
    + V_us_sq * ( X1.u * X2.s + X1.s * X2.u )
    + V_ub_sq * ( X1.u * X2.b + X1.b * X2.u )
    + V_cd_sq * ( X1.c * X2.dbar + X1.dbar * X2.c )
    + V_cs_sq * ( X1.c * X2.s + X1.s * X2.c )
    + V_cb_sq * ( X1.c * X2.b + X1.b * X2.c ) );
    }
    else if ((fWp==1) && (c==pp)) { return pref * (
      V_ud_sq * ( X1.d * X2.u + X1.ubar * X2.dbar )
    + V_us_sq * ( X1.s * X2.u + X1.ubar * X2.s )
    + V_ub_sq * ( X1.b * X2.u + X1.ubar * X2.b )
    + V_cd_sq * ( X1.d * X2.c + X1.c * X2.dbar )
    + V_cs_sq * ( X1.s * X2.c + X1.c * X2.s )
    + V_cb_sq * ( X1.b * X2.c + X1.c * X2.b ) );
    }
    else if ((fWp==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( X1.d * X2.ubar + X1.ubar * X2.d )
    + V_us_sq * ( X1.s * X2.ubar + X1.ubar * X2.s )
    + V_ub_sq * ( X1.b * X2.ubar + X1.ubar * X2.b )
    + V_cd_sq * ( X1.d * X2.c + X1.c * X2.d )
    + V_cs_sq * ( X1.s * X2.c + X1.c * X2.s )
    + V_cb_sq * ( X1.b * X2.c + X1.c * X2.b ) );
    }
    else return 0.;
  }
  else if ((c==pp) || (c==ppbar) || c==piso) { return qq_CE1_lumi(X1,X2,c) ; }
  else return 0.;
}

/* "CF" and luminosity for identical-quark (or identical antiquark)
   exchange contribution. 
   "CF" is same as "CE1" and "CE2" in (Z,gamma) case; but not W^\pm case. */

double qq_CF_lumi(const pdfArray& X1, const pdfArray& X2, collider c){
  if ((f_qq == 0) || (X1.x >= 1.) || (X2.x >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
    if ((fWp==1) && (c==pp)) { return
      uWf * X1.u * X2.u + cWf * X1.c * X2.c
    + dWf * X1.dbar * X2.dbar + sWf * X1.s * X2.s
    + bWf * X1.b * X2.b ;
    }
    if ((fWp==1) && (c==ppbar)) { return 
      uWf * X1.u * X2.ubar + cWf * X1.c * X2.c
    + dWf * X1.dbar * X2.d + sWf * X1.s * X2.s
    + bWf * X1.b * X2.b ;
    }
    if ((fWm==1) && (c==pp)) { return
      uWf * X1.ubar * X2.ubar + cWf * X1.c * X2.c
    + dWf * X1.d * X2.d + sWf * X1.s * X2.s
    + bWf * X1.b * X2.b ;
    }
    if ((fWm==1) && (c==ppbar)) { return
      uWf * X1.ubar * X2.u + cWf * X1.c * X2.c
    + dWf * X1.d * X2.dbar + sWf * X1.s * X2.s
    + bWf * X1.b * X2.b ;
    }
    return 0.;
  }
  else if ((c==pp) || (c==ppbar) || c==piso) { return qq_CE1_lumi(X1,X2,c) ; }
  else return 0.;
}
