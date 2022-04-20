/* =============================================================
   pdf luminosity functions required for Z-gamma^* total cross section
   and rapidity distributions.
   ============================================================= */
#include <iomanip>
#include "NAClasses.h"
#include "pdf.h"  // pdf distributions
#include "Vlumifns.h"
#include "EWparams.h"

// Fudge parameters to weight qqbar, qg and gg contributions
// (all terms, LO, NLO, NNLO!):
int f_qqbar, f_qg, f_gg, f_qq;
// If f_X = 0, then channel X is set to zero. 
// If f_X = 1, then it is included at full strength.
// Thus, normally all should be set to 1.
// For (qqbar only), set f_qqbar = 1, f_qg = f_gg = f_qq = 0, etc.

// procedures to decide which parton channels are being computed:

void compute_all(){ 
 f_qqbar = 1; f_qg = 1; f_gg = 1;  f_qq = 1;
   std::cout << " Computing all initial-state channels " << std::endl << std::endl;
 }

void compute_qqbar(){ 
 f_qqbar = 1; f_qg = 0; f_gg = 0;  f_qq = 0; 
 std::cout << " Computing only q-qbar channel " << std::endl << std::endl;
 }

void compute_qg(){ 
 f_qqbar = 0; f_qg = 1; f_gg = 0;  f_qq = 0; 
 std::cout << " Computing only q-g (qbar-g) channels " << std::endl << std::endl;
 }

void compute_gg(){ 
 f_qqbar = 0; f_qg = 0; f_gg = 1;  f_qq = 0;  
 std::cout << " Computing only g-g channel " << std::endl << std::endl;
 }

void compute_qq(){ 
 f_qqbar = 0; f_qg = 0; f_gg = 0;  f_qq = 1;  
 std::cout << " Computing only q_i q_j (qbar_j) and identical-q channels " 
      << std::endl << std::endl;
 }

void compute_qqbar_plus_qq(){ 
 f_qqbar = 1; f_qg = 0; f_gg = 0;  f_qq = 1;  
 std::cout << " Computing q-qbar, plus q_i q_j (qbar_j), and identical-q channels " 
      << std::endl << std::endl;
 }

// Parameters to weight gamma^* vs. Z vs. W exchange
// (all terms, LO, NLO, NNLO!):
double F_gamma, F_Zgamma, F_Z, F_W;

// flag W+, W-:
int fWp = 0;   int fWm = 0;

// storage of coupling parameters:
double vsqasq_u, vsqasq_d, NF_f, vuvu, vdvd, vuvd; 
double auau, adad, auad, aaf_u, aaf_d;
// these contain CKM angles (for W production):
double sum_W_f, uWf, cWf, dWf, sWf, bWf;

// Procedure to decide which of Z/gamma^*/W exchanges are being computed:

void setV(exchange E, double Q, double alpha, double Nf, int f_quiet){
  if (E == gamma_only){
    F_gamma = 1.;   F_Zgamma = 0.;   F_Z = 0.;   F_W = 0.;
    fWp = 0;   fWm = 0;
    if (f_quiet==0) {
    std::cout << " Computing only gamma^* exchange " << std::endl << std::endl;
    }
  }
  else if (E == Zgamma_interf){
    F_gamma = 0.;   F_Zgamma = N_Zgamma(Q,alpha);   F_Z = 0.;   F_W = 0.;
    if (f_quiet==0) {
    std::cout << " Computing only Z-gamma^* interference " << std::endl << std::endl;
    }
  }
  else if (E == Z_only){
    F_gamma = 0.;   F_Zgamma = 0.;   F_Z = N_Z(Q,alpha);   F_W = 0.;
    fWp = 0;   fWm = 0;
    if (f_quiet==0) {
    std::cout << " Computing only Z exchange " << std::endl << std::endl;
    }
  }
  else if (E == Zgamma){
    F_gamma = 1.;   F_Zgamma = N_Zgamma(Q,alpha);   F_Z = N_Z(Q,alpha);
    F_W = 0.;  
    fWp = 0;   fWm = 0;
    if (f_quiet==0) {
    std::cout << " Computing all terms: gamma + Z + interference " << std::endl << std::endl;
    }
  }
  else if (E == Wplus){
    F_gamma = 0.;   F_Zgamma = 0.;   F_Z = 0.;   F_W = N_W(Q,alpha);
    fWp = 1;   fWm = 0;
    if (f_quiet==0) {
    std::cout << " Computing W^+ production " << std::endl << std::endl;
    }
  }
  else if (E == Wminus){
    F_gamma = 0.;   F_Zgamma = 0.;   F_Z = 0.;   F_W = N_W(Q,alpha);
    fWp = 0;  fWm = 1;
    if (f_quiet==0) {
    std::cout << " Computing W^- production " << std::endl << std::endl;
    }
  }
  else{ std::cout << " Invalid choice for exchange argument to setV " << std::endl; }
// compute all the couplings required for luminosity functions:
  vsqasq_u = VsqAsq_u();  vsqasq_d = VsqAsq_d();  
  NF_f = NFf(int(Nf));  sum_W_f = Sum_W_f(int(Nf));
  uWf = UWf(int(Nf));   cWf = CWf(int(Nf));  
  dWf = DWf(int(Nf));   sWf = SWf(int(Nf));   bWf = BWf(int(Nf)); 
  vuvu = VuVu();   vdvd = VdVd();   vuvd = VuVd(); 
  auau = AuAu();   adad = AdAd();   auad = AuAd(); 
  aaf_u = AAf_u(int(Nf));  aaf_d = AAf_d(int(Nf));  
}


// the gamma-Z interference normalization factor:

double N_Zgamma(double Q, double alpha){
  double Qsq = Q*Q;  double mZsq = mZ*mZ;
  return (1.-4.*sstW)/8./sstW/cstW
    * Qsq*(Qsq-mZsq)/((Qsq-mZsq)*(Qsq-mZsq) + Gamma_Z*Gamma_Z * mZsq);
}

// the pure-Z normalization factor:

double N_Z(double Q, double alpha){
  double Qsq = Q*Q;  double mZsq = mZ*mZ;
  return 3./16./sstW/cstW/alpha * Gamma_Z * Br_l_Z/mZ
    * Qsq*Qsq/((Qsq-mZsq)*(Qsq-mZsq) + Gamma_Z*Gamma_Z * mZsq);
}

// the W normalization factor:

double N_W(double Q, double alpha){
  double Qsq = Q*Q;  double mWsq = mW*mW;
  return 3./4./sstW/alpha * Gamma_W * Br_l_W/mW
    * Qsq*Qsq/((Qsq-mWsq)*(Qsq-mWsq) + Gamma_W*Gamma_W * mWsq);
}

// Computation of coupling parameters:

// The (v^2+a^2) functions appearing in the main q-qbar PLUS qbar-q
//  contribution to Z-gamma^* production (plus several other channels):

double VsqAsq_u(){
 return F_gamma * Q_u*Q_u + F_Zgamma * Q_u * v_u_Z 
      + F_Z * (v_u_Z*v_u_Z + a_u_Z*a_u_Z) ;
}

double VsqAsq_d(){
 return F_gamma * Q_d*Q_d + F_Zgamma * Q_d * v_d_Z 
      + F_Z * (v_d_Z*v_d_Z + a_d_Z*a_d_Z) ;
}

// The v_i*v_j functions appearing in Z-gamma^* production:

double VuVu(){
 return F_gamma * Q_u*Q_u + F_Zgamma * Q_u * v_u_Z + F_Z * v_u_Z*v_u_Z ;
}

double VdVd(){
 return F_gamma * Q_d*Q_d + F_Zgamma * Q_d * v_d_Z + F_Z * v_d_Z*v_d_Z ;
}

double VuVd(){
  return F_gamma * Q_u*Q_d + F_Zgamma * 0.5 * (Q_u * v_d_Z + Q_d * v_u_Z) 
       + F_Z * v_u_Z*v_d_Z ;
}

// The a_i*a_j functions appearing in Z-gamma^* production:

double AAf_u(int Nf){ return F_Z * a_u_Z * Sum_a_f(Nf); }
double AAf_d(int Nf){ return F_Z * a_d_Z * Sum_a_f(Nf); }

double AuAu(){ return F_Z * a_u_Z*a_u_Z ; }
double AdAd(){ return F_Z * a_d_Z*a_d_Z ; }
double AuAd(){ return F_Z * a_u_Z*a_d_Z ; }

//======== Auxiliary functions ============================

// The quantity "NFf" is used for Z-gamma^* emission off a final state 
// q-\bar{q} pair; it returns \sum_{f=1}^{Nf} e_f^2 in the photon case,
// and the (v^2+a^2) generalization in the Z-gamma^* case:

double NFf(int Nf){ 
  if (Nf==1){ return vsqasq_u; }
  else if (Nf==2){ return vsqasq_u + vsqasq_d; }
  else if (Nf==3){ return vsqasq_u + 2.*vsqasq_d; }
  else if (Nf==4){ return 2.*vsqasq_u + 2.*vsqasq_d; }
  else if (Nf==5){ return 2.*vsqasq_u + 3.*vsqasq_d; }
  else if (Nf==6){ return 3.*vsqasq_u + 3.*vsqasq_d; }
  else { return 0.; }
}

// The quantity "Sum_a_f" returns \sum_{f=1}^{Nf} a_f_Z:

double Sum_a_f(int Nf){ 
  if ((Nf==2) || (Nf==4) || (Nf==6)){ return 0.; }
  else if (Nf==1){ return a_u_Z; }
  else if ((Nf==3) || (Nf==5)){ return a_d_Z; }
  else { return 0.; }
}

// The quantity "Sum_W_f" returns
//  F_W * ( v_W^2 + a_W^2) * \sum_{f,f'=1}^{Nf} |V_{ff'}|^2 :

double Sum_W_f(int Nf){ 
  double temp;
  if (Nf==2){ temp = V_ud_sq; }
  else if (Nf==3){ temp = 1. - V_ub_sq; }   // by unitarity....
  else if (Nf==4){ temp = 2. - V_ub_sq - V_cb_sq; }
  else if (Nf==5){ temp = 2.; }
  else if (Nf==6){ temp = 3.; }  // better be well above m_t !!!
  else { temp = 0.; }     // includes Nf=1 case
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

// The quantities "QWf" return
//  F_W * ( v_W^2 + a_W^2) * \sum_{f=1}^{Nf} |V_{Qf}|^2 :

double UWf(int Nf){ 
  double temp;
  if (Nf==1){ temp = 0.; }
  else if (Nf==2){ temp = V_ud_sq; }
  else if ((Nf==3) || (Nf==4)){ temp = 1. - V_ub_sq; }   // by unitarity....
  else temp = 1.;
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

double CWf(int Nf){ 
  double temp;
  if (Nf==1){ temp = 0.; }
  else if (Nf==2){ temp = V_cd_sq; }
  else if ((Nf==3) || (Nf==4)){ temp = 1. - V_cb_sq; }   // by unitarity....
  else temp = 1.;
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

double DWf(int Nf){ 
  double temp;
  if ((Nf==1) || (Nf==2)|| (Nf==3)){ temp = V_ud_sq; }
  else if ((Nf==4) || (Nf==5)){ temp = 1. - V_td_sq; }
  else temp = 1.;
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

double SWf(int Nf){ 
  double temp;
  if ((Nf==1) || (Nf==2)|| (Nf==3)){ temp = V_us_sq; }
  else if ((Nf==4) || (Nf==5)){ temp = 1. - V_ts_sq; }
  else temp = 1.;
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

double BWf(int Nf){ 
  double temp;
  if ((Nf==1) || (Nf==2)|| (Nf==3)){ temp = V_ub_sq; }
  else if ((Nf==4) || (Nf==5)){ temp = V_ub_sq + V_cb_sq; }
  else temp = 1.;
  return F_W * ( v_W*v_W + a_W*a_W ) * temp;
}

#if USE_OLD

//======== Luminosity functions ============================

/* The (v^2+a^2)-weighted parton luminosity function 
   that always appears in main q-qbar PLUS qbar-q 
   contribution to Z-gamma^* production.
   For archaic reasons we use a "process" flag to weight the "charges"
   "Drell-Yan" (for process = 1) -> (v^2+a^2)
   or "gluon" (for "NFf" term) (for process = 0) -> no-charge weighting here.
   Collider can be pp (for collider = 1) or p-pbar (for collider = 2). */

double qqbar_lumi(double x1, double x2, double muF, process p, collider c){
  if ((f_qqbar == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  double pref, pref_u, pref_d;
// The "DY" type W^\pm cases are cumbersome because of CKM elements:
  if ( ((fWp==1) || (fWm==1)) && (p==DY) ){
    pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWp==1) && (c==pp)){ return pref * ( 
       u(x1,muF) * ( V_ud_sq * dbar(x2,muF) 
	           + V_us_sq * str(x2,muF) + V_ub_sq * bot(x2,muF) )
     + chm(x1,muF) * ( V_cd_sq * dbar(x2,muF) 
	           + V_cs_sq * str(x2,muF) + V_cb_sq * bot(x2,muF) )
     + u(x2,muF) * ( V_ud_sq * dbar(x1,muF) 
	           + V_us_sq * str(x1,muF) + V_ub_sq * bot(x1,muF) )
     + chm(x2,muF) * ( V_cd_sq * dbar(x1,muF) 
		     + V_cs_sq * str(x1,muF) + V_cb_sq * bot(x1,muF) ) );
    }
//
    if ((fWp==1) && (c==ppbar)){ return pref * ( 
       u(x1,muF) * ( V_ud_sq * d(x2,muF) 
	           + V_us_sq * str(x2,muF) + V_ub_sq * bot(x2,muF) )
     + chm(x1,muF) * ( V_cd_sq * d(x2,muF) 
	           + V_cs_sq * str(x2,muF) + V_cb_sq * bot(x2,muF) )
     + ubar(x2,muF) * ( V_ud_sq * dbar(x1,muF) 
	           + V_us_sq * str(x1,muF) + V_ub_sq * bot(x1,muF) )
     + chm(x2,muF) * ( V_cd_sq * dbar(x1,muF) 
		     + V_cs_sq * str(x1,muF) + V_cb_sq * bot(x1,muF) ) );
    }
//
    if ((fWm==1) && (c==pp)){ return pref * ( 
       ubar(x1,muF) * ( V_ud_sq * d(x2,muF) 
	           + V_us_sq * str(x2,muF) + V_ub_sq * bot(x2,muF) )
     + chm(x1,muF) * ( V_cd_sq * d(x2,muF) 
	           + V_cs_sq * str(x2,muF) + V_cb_sq * bot(x2,muF) )
     + ubar(x2,muF) * ( V_ud_sq * d(x1,muF) 
	           + V_us_sq * str(x1,muF) + V_ub_sq * bot(x1,muF) )
     + chm(x2,muF) * ( V_cd_sq * d(x1,muF) 
		     + V_cs_sq * str(x1,muF) + V_cb_sq * bot(x1,muF) ) );
    }
//
    if ((fWm==1) && (c==ppbar)){ return pref * ( 
       ubar(x1,muF) * ( V_ud_sq * dbar(x2,muF) 
	           + V_us_sq * str(x2,muF) + V_ub_sq * bot(x2,muF) )
     + chm(x1,muF) * ( V_cd_sq * dbar(x2,muF) 
	           + V_cs_sq * str(x2,muF) + V_cb_sq * bot(x2,muF) )
     + u(x2,muF) * ( V_ud_sq * d(x1,muF) 
	           + V_us_sq * str(x1,muF) + V_ub_sq * bot(x1,muF) )
     + chm(x2,muF) * ( V_cd_sq * d(x1,muF) 
		     + V_cs_sq * str(x1,muF) + V_cb_sq * bot(x1,muF) ) );
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
     pref_u * ( u(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * u(x2,muF) 
		                + 2. * chm(x1,muF) * chm(x2,muF) )
   + pref_d * ( d(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * d(x2,muF) 
	                        + 2. * str(x1,muF) * str(x2,muF) 
		                + 2. * bot(x1,muF) * bot(x2,muF) )
 ; }
  else if (c==ppbar) { return 
     pref_u * ( u(x1,muF) * u(x2,muF) + ubar(x1,muF) * ubar(x2,muF) 
                                + 2. * chm(x1,muF) * chm(x2,muF) )
   + pref_d * ( d(x1,muF) * d(x2,muF) + dbar(x1,muF) * dbar(x2,muF) 
	                        + 2. * str(x1,muF) * str(x2,muF) 
                                + 2. * bot(x1,muF) * bot(x2,muF) ) ; }
  else if (c==piso) { return 
      pref_u * ( u(x1,muF) * (ubar(x2,muF)+dbar(x2,muF))/2. + ubar(x1,muF) * (u(x2,muF)+d(x2,muF))/2.
		 + 2. * chm(x1,muF) * chm(x2,muF) )
      + pref_d * ( d(x1,muF) * (dbar(x2,muF)+ubar(x2,muF))/2. + dbar(x1,muF) * (d(x2,muF)+u(x2,muF))/2.
		   + 2. * str(x1,muF) * str(x2,muF) 
		   + 2. * bot(x1,muF) * bot(x2,muF) ) ; }
  else return 0.;
}

/* For W^\pm production, the "BC" terms get a different type of luminosity
   function.  This function is equal to 2 * qbarq_lumi(DY) for the (Z,gamma)
   case. */

double qqbar_BC_lumi(double x1, double x2, double muF, collider c){
  if ((f_qqbar == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
// The W^+ and W^- cases are the same here:
  if ((fWp==1) || (fWm==1)) { 
    if (c==pp) { return
        uWf * ( u(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * u(x2,muF) ) 
      + 2. * cWf * chm(x1,muF) * chm(x2,muF)
      + dWf * ( d(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * d(x2,muF) )
      + 2. * sWf * str(x1,muF) * str(x2,muF) 
      + 2. * bWf * bot(x1,muF) * bot(x2,muF) ;
    }
    if (c==ppbar) { return
        uWf * ( u(x1,muF) * u(x2,muF) + ubar(x1,muF) * ubar(x2,muF) ) 
      + 2. * cWf * chm(x1,muF) * chm(x2,muF)
      + dWf * ( d(x1,muF) * d(x2,muF) + dbar(x1,muF) * dbar(x2,muF) )
      + 2. * sWf * str(x1,muF) * str(x2,muF) 
      + 2. * bWf * bot(x1,muF) * bot(x2,muF) ;
    }
  }
  else if ((c==pp) || (c==ppbar)) { return 2. * qqbar_lumi(x1,x2,muF,DY,c); }
  return 0.;
}
    
// q-\bar{q} luminosity function required for the "AB axial" terms.
// It should return 0 for W^\pm production.

double qqbar_ax_lumi(double x1, double x2, double muF, collider c){
  if ((f_qqbar == 0) || (fWp==1) || (fWm==1) || (x1 >= 1.) || (x2 >= 1.)) { 
    return 0.; 
  }
  if (c==pp) { return 
     aaf_u * ( u(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * u(x2,muF) 
                                + 2. * chm(x1,muF) * chm(x2,muF) )
   + aaf_d * ( d(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * d(x2,muF) 
                                + 2. * str(x1,muF) * str(x2,muF) 
                                + 2. * bot(x1,muF) * bot(x2,muF) ) ; }
  else if (c==ppbar) { return 
     aaf_u * ( u(x1,muF) * u(x2,muF) + ubar(x1,muF) * ubar(x2,muF) 
                                + 2. * chm(x1,muF) * chm(x2,muF) )
   + aaf_d * ( d(x1,muF) * d(x2,muF) + dbar(x1,muF) * dbar(x2,muF) 
	                        + 2. * str(x1,muF) * str(x2,muF) 
                                + 2. * bot(x1,muF) * bot(x2,muF) ) ; }
  else return 0.;
}

/* The quark-charge weighted parton luminosity function 
   that always appears in the qg PLUS qbar-g (assuming C invariance!) 
   initiated contribution to W+, W-, or gamma^* / Z production
   for pp (for collider = 1) or p-pbar (for collider = 2). */ 

double qg_lumi(double x1, double x2, double muF, collider c){
  if ((f_qg == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if (fWp==1){ return 
      ( uWf * u(x1,muF) + cWf * chm(x1,muF)
      + dWf * dbar(x1,muF) + sWf * str(x1,muF) + bWf * bot(x1,muF) )
      * g(x2,muF);
  }
  else if (fWm==1){ return 
      ( uWf * ubar(x1,muF) + cWf * chm(x1,muF)
      + dWf * d(x1,muF) + sWf * str(x1,muF) + bWf * bot(x1,muF) )
      * g(x2,muF);
  }
  else if ((c==pp) || (c==ppbar)) { return 
     ( vsqasq_u * (u(x1,muF) + ubar(x1,muF) + 2. * chm(x1,muF))
     + vsqasq_d * (d(x1,muF) + dbar(x1,muF) 
                  + 2. * str(x1,muF) + 2. * bot(x1,muF)) )
      * g(x2,muF);
  }
  else return 0.;
}

/* And similarly for gq PLUS g-qbar (not quite so simple in W case;
   can't just use x1 <-> x2 because of p vs. pbar distinction). */

double gq_lumi(double x1, double x2, double muF, collider c){
  if ((f_qg == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if (((fWp==1) && (c==pp)) || ((fWm==1) && (c==ppbar))) { return 
      ( uWf * u(x2,muF) + cWf * chm(x2,muF)
      + dWf * dbar(x2,muF) + sWf * str(x2,muF) + bWf * bot(x2,muF) )
      * g(x1,muF);
  }
  else if (((fWm==1) && (c==pp)) || ((fWp==1) && (c==ppbar))) { return 
      ( uWf * ubar(x2,muF) + cWf * chm(x2,muF)
      + dWf * d(x2,muF) + sWf * str(x2,muF) + bWf * bot(x2,muF) )
      * g(x1,muF);
  }
  else if ((c==pp) || (c==ppbar)) { return 
     ( vsqasq_u * (u(x2,muF) + ubar(x2,muF) + 2. * chm(x2,muF))
     + vsqasq_d * (d(x2,muF) + dbar(x2,muF) 
                  + 2. * str(x2,muF) + 2. * bot(x2,muF)) )
      * g(x1,muF);
  }
  else return 0.;
}

/* The gluon-gluon parton luminosity function (same for pp or p-pbar).
 For (Z,gamma) production it always appears multiplied by "NFf";
 for W^\pm production it is multiplied by sum_W_f.  */

double gg_lumi(double x1, double x2, double muF, collider c){
  if ((f_gg == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  double pref;
  if ((fWp==1) || (fWm==1)) { pref = sum_W_f; }
  else { pref = NF_f; }
  return pref * g(x1,muF) * g(x2,muF) ;
}

/* The (generalized) quark-charge weighted parton luminosity function 
   that appears in the "C^2" (v_1^2 + a_1^2) q_i q_j (qbar_j) PLUS charge-conj
   initiated contributions to Drell-Yan production.
   Here (v_1^2 +a_1^2) means both quark charges are associated with hadron 1.
   The function is the same for both pp (collider = 1) 
   and p-pbar (collider = 2). */

double qq_11_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if (fWp==1){ return 
      ( uWf * u(x1,muF) + cWf * chm(x1,muF)
      + dWf * dbar(x1,muF) + sWf * str(x1,muF) + bWf * bot(x1,muF) )
   * ( u(x2,muF) + ubar(x2,muF) + 2.*chm(x2,muF)
     + d(x2,muF) + dbar(x2,muF) + 2.*str(x2,muF) + 2.*bot(x2,muF) ); 
  }
  else if (fWm==1){ return 
      ( uWf * ubar(x1,muF) + cWf * chm(x1,muF)
      + dWf * d(x1,muF) + sWf * str(x1,muF) + bWf * bot(x1,muF) )
   * ( u(x2,muF) + ubar(x2,muF) + 2.*chm(x2,muF)
     + d(x2,muF) + dbar(x2,muF) + 2.*str(x2,muF) + 2.*bot(x2,muF) ); 
  }
  else if ((c==pp) || (c==ppbar)){ return 
    ( vsqasq_u * ( u(x1,muF) + ubar(x1,muF) + 2.*chm(x1,muF) )
    + vsqasq_d * ( d(x1,muF) + dbar(x1,muF) + 2.*str(x1,muF) + 2.*bot(x1,muF) ) )  
   * ( u(x2,muF) + ubar(x2,muF) + 2.*chm(x2,muF)
     + d(x2,muF) + dbar(x2,muF) + 2.*str(x2,muF) + 2.*bot(x2,muF) );
  }
  else return 0.;
}

/* Same function, except for e_2^2, i.e. both quark charges are associated
   with hadron 2.  Get by letting (x1 <-> x2) in qq_11_lumi, except for
   the case of W^\pm production... */

double qq_22_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if (((fWp==1) && (c==pp)) || ((fWm==1) && (c==ppbar))) { return 
      ( uWf * u(x2,muF) + cWf * chm(x2,muF)
      + dWf * dbar(x2,muF) + sWf * str(x2,muF) + bWf * bot(x2,muF) )
   * ( u(x1,muF) + ubar(x1,muF) + 2.*chm(x1,muF)
     + d(x1,muF) + dbar(x1,muF) + 2.*str(x1,muF) + 2.*bot(x1,muF) ); 
  }
  else if (((fWm==1) && (c==pp)) || ((fWp==1) && (c==ppbar))) { return 
      ( uWf * ubar(x2,muF) + cWf * chm(x2,muF)
      + dWf * d(x2,muF) + sWf * str(x2,muF) + bWf * bot(x2,muF) )
   * ( u(x1,muF) + ubar(x1,muF) + 2.*chm(x1,muF)
     + d(x1,muF) + dbar(x1,muF) + 2.*str(x1,muF) + 2.*bot(x1,muF) ); 
  }
  else if ((c==pp) || (c==ppbar)){ return 
    ( vsqasq_u * ( u(x2,muF) + ubar(x2,muF) + 2.*chm(x2,muF) )
    + vsqasq_d * ( d(x2,muF) + dbar(x2,muF) + 2.*str(x2,muF) + 2.*bot(x2,muF) ) )  
   * ( u(x1,muF) + ubar(x1,muF) + 2.*chm(x1,muF)
     + d(x1,muF) + dbar(x1,muF) + 2.*str(x1,muF) + 2.*bot(x1,muF) );
  }
  else return 0.;
}

/* Same function, except for v_1*v_2, 
   i.e. one quark vector charge is associated with each hadron. 
   Note minus signs under q <-> qbar (vector case).
   Hence the nonvalence quarks do not contribute.
   Overall (-) is so q-qbar enters with plus sign.
   Like qqbar_ax_lumi, this function should return 0 for W^\pm production. */

double qq_12_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (fWp==1) || (fWm==1) || (x1 >= 1.) || (x2 >= 1.)) { 
    return 0.;
  }
  double temp = - ( 
      vuvu * (u(x1,muF) - ubar(x1,muF)) * (u(x2,muF) - ubar(x2,muF))
    + vuvd * ( (u(x1,muF) - ubar(x1,muF)) * (d(x2,muF) - dbar(x2,muF))
	     + (d(x1,muF) - dbar(x1,muF)) * (u(x2,muF) - ubar(x2,muF)) )
    + vdvd * (d(x1,muF) - dbar(x1,muF)) * (d(x2,muF) - dbar(x2,muF)) ) ;
  if (c==pp){ return temp; }
  else if (c==ppbar) { return - temp; }  
  else return 0.;
}

/* Same function as qq_12_lumi, except contains axial vector coupling
  product, a_1*a_2. It also returns 0 for W^\pm production. */

double qq_12_ax_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (fWp==1) || (fWm==1) || (x1 >= 1.) || (x2 >= 1.)) {
     return 0.;
  }
  if ((c==pp) || (c==ppbar)){ return 
      auau * (u(x1,muF) + ubar(x1,muF) + 2. * chm(x1,muF)) 
           * (u(x2,muF) + ubar(x2,muF) + 2. * chm(x2,muF))
    + auad * ( (u(x1,muF) + ubar(x1,muF) + 2. * chm(x1,muF))
             * (d(x2,muF) + dbar(x2,muF) + 2. * str(x2,muF) + 2. * bot(x2,muF))
             + (d(x1,muF) + dbar(x1,muF) + 2. * str(x1,muF) + 2. * bot(x1,muF))
	     * (u(x2,muF) + ubar(x2,muF) + 2. * chm(x2,muF)) )
    + adad * (d(x1,muF) + dbar(x1,muF) + 2. * str(x1,muF) + 2. * bot(x1,muF))
           * (d(x2,muF) + dbar(x2,muF) + 2. * str(x2,muF) + 2. * bot(x2,muF)) ;
  }
  else return 0.;
}

/* "CE1" and "CE2" luminosities for identical-quark (or identical antiquark)
   exchange contribution. 
   "CE2" is same as "CE1" in (Z,gamma) case; but not W^\pm case. */

double qq_CE1_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
  double pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWp==1) && (c==pp)) { return pref * (
      V_ud_sq * ( u(x1,muF) * d(x2,muF) + dbar(x1,muF) * ubar(x2,muF) )
    + V_us_sq * ( u(x1,muF) * str(x2,muF) + str(x1,muF) * ubar(x2,muF) )
    + V_ub_sq * ( u(x1,muF) * bot(x2,muF) + bot(x1,muF) * ubar(x2,muF) )
    + V_cd_sq * ( chm(x1,muF) * d(x2,muF) + dbar(x1,muF) * chm(x2,muF) )
    + V_cs_sq * ( chm(x1,muF) * str(x2,muF) + str(x1,muF) * chm(x2,muF) )
    + V_cb_sq * ( chm(x1,muF) * bot(x2,muF) + bot(x1,muF) * chm(x2,muF) ) );
    }
    if ((fWp==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( u(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * u(x2,muF) )
    + V_us_sq * ( u(x1,muF) * str(x2,muF) + str(x1,muF) * u(x2,muF) )
    + V_ub_sq * ( u(x1,muF) * bot(x2,muF) + bot(x1,muF) * u(x2,muF) )
    + V_cd_sq * ( chm(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * chm(x2,muF) )
    + V_cs_sq * ( chm(x1,muF) * str(x2,muF) + str(x1,muF) * chm(x2,muF) )
    + V_cb_sq * ( chm(x1,muF) * bot(x2,muF) + bot(x1,muF) * chm(x2,muF) ) );
    }
    if ((fWm==1) && (c==pp)) { return pref * (
      V_ud_sq * ( d(x1,muF) * u(x2,muF) + ubar(x1,muF) * dbar(x2,muF) )
    + V_us_sq * ( str(x1,muF) * u(x2,muF) + ubar(x1,muF) * str(x2,muF) )
    + V_ub_sq * ( bot(x1,muF) * u(x2,muF) + ubar(x1,muF) * bot(x2,muF) )
    + V_cd_sq * ( d(x1,muF) * chm(x2,muF) + chm(x1,muF) * dbar(x2,muF) )
    + V_cs_sq * ( str(x1,muF) * chm(x2,muF) + chm(x1,muF) * str(x2,muF) )
    + V_cb_sq * ( bot(x1,muF) * chm(x2,muF) + chm(x1,muF) * bot(x2,muF) ) );
    }
    if ((fWm==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( d(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * d(x2,muF) )
    + V_us_sq * ( str(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * str(x2,muF) )
    + V_ub_sq * ( bot(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * bot(x2,muF) )
    + V_cd_sq * ( d(x1,muF) * chm(x2,muF) + chm(x1,muF) * d(x2,muF) )
    + V_cs_sq * ( str(x1,muF) * chm(x2,muF) + chm(x1,muF) * str(x2,muF) )
    + V_cb_sq * ( bot(x1,muF) * chm(x2,muF) + chm(x1,muF) * bot(x2,muF) ) );
    }
    return 0.;
  }
  else if (c==pp){ return
     vsqasq_u * ( u(x1,muF) * u(x2,muF) + ubar(x1,muF) * ubar(x2,muF) 
                + 2. * chm(x1,muF) * chm(x2,muF) )
   + vsqasq_d * ( d(x1,muF) * d(x2,muF) + dbar(x1,muF) * dbar(x2,muF) 
        + 2. * str(x1,muF) * str(x2,muF) + 2. * bot(x1,muF) * bot(x2,muF) ); }
  else if (c==ppbar) { return
     vsqasq_u * ( u(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * u(x2,muF) 
               + 2. * chm(x1,muF) * chm(x2,muF) )
   + vsqasq_d * ( d(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * d(x2,muF) 
        + 2. * str(x1,muF) * str(x2,muF) + 2. * bot(x1,muF) * bot(x2,muF) ); }
  else return 0.;
}

double qq_CE2_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
  double pref = F_W * ( v_W*v_W + a_W*a_W );
    if ((fWm==1) && (c==pp)) { return pref * (
      V_ud_sq * ( u(x1,muF) * d(x2,muF) + dbar(x1,muF) * ubar(x2,muF) )
    + V_us_sq * ( u(x1,muF) * str(x2,muF) + str(x1,muF) * ubar(x2,muF) )
    + V_ub_sq * ( u(x1,muF) * bot(x2,muF) + bot(x1,muF) * ubar(x2,muF) )
    + V_cd_sq * ( chm(x1,muF) * d(x2,muF) + dbar(x1,muF) * chm(x2,muF) )
    + V_cs_sq * ( chm(x1,muF) * str(x2,muF) + str(x1,muF) * chm(x2,muF) )
    + V_cb_sq * ( chm(x1,muF) * bot(x2,muF) + bot(x1,muF) * chm(x2,muF) ) );
    }
    else if ((fWm==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( u(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * u(x2,muF) )
    + V_us_sq * ( u(x1,muF) * str(x2,muF) + str(x1,muF) * u(x2,muF) )
    + V_ub_sq * ( u(x1,muF) * bot(x2,muF) + bot(x1,muF) * u(x2,muF) )
    + V_cd_sq * ( chm(x1,muF) * dbar(x2,muF) + dbar(x1,muF) * chm(x2,muF) )
    + V_cs_sq * ( chm(x1,muF) * str(x2,muF) + str(x1,muF) * chm(x2,muF) )
    + V_cb_sq * ( chm(x1,muF) * bot(x2,muF) + bot(x1,muF) * chm(x2,muF) ) );
    }
    else if ((fWp==1) && (c==pp)) { return pref * (
      V_ud_sq * ( d(x1,muF) * u(x2,muF) + ubar(x1,muF) * dbar(x2,muF) )
    + V_us_sq * ( str(x1,muF) * u(x2,muF) + ubar(x1,muF) * str(x2,muF) )
    + V_ub_sq * ( bot(x1,muF) * u(x2,muF) + ubar(x1,muF) * bot(x2,muF) )
    + V_cd_sq * ( d(x1,muF) * chm(x2,muF) + chm(x1,muF) * dbar(x2,muF) )
    + V_cs_sq * ( str(x1,muF) * chm(x2,muF) + chm(x1,muF) * str(x2,muF) )
    + V_cb_sq * ( bot(x1,muF) * chm(x2,muF) + chm(x1,muF) * bot(x2,muF) ) );
    }
    else if ((fWp==1) && (c==ppbar)) { return pref * (
      V_ud_sq * ( d(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * d(x2,muF) )
    + V_us_sq * ( str(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * str(x2,muF) )
    + V_ub_sq * ( bot(x1,muF) * ubar(x2,muF) + ubar(x1,muF) * bot(x2,muF) )
    + V_cd_sq * ( d(x1,muF) * chm(x2,muF) + chm(x1,muF) * d(x2,muF) )
    + V_cs_sq * ( str(x1,muF) * chm(x2,muF) + chm(x1,muF) * str(x2,muF) )
    + V_cb_sq * ( bot(x1,muF) * chm(x2,muF) + chm(x1,muF) * bot(x2,muF) ) );
    }
    else return 0.;
  }
  else if ((c==pp) || (c==ppbar)) { return qq_CE1_lumi(x1,x2,muF,c) ; }
  else return 0.;
}

/* "CF" and luminosity for identical-quark (or identical antiquark)
   exchange contribution. 
   "CF" is same as "CE1" and "CE2" in (Z,gamma) case; but not W^\pm case. */

double qq_CF_lumi(double x1, double x2, double muF, collider c){
  if ((f_qq == 0) || (x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  if ((fWp==1) || (fWm==1)) {   // W^\pm case cumbersome due to CKM angles:
    if ((fWp==1) && (c==pp)) { return
      uWf * u(x1,muF) * u(x2,muF) + cWf * chm(x1,muF) * chm(x2,muF)
    + dWf * dbar(x1,muF) * dbar(x2,muF) + sWf * str(x1,muF) * str(x2,muF)
    + bWf * bot(x1,muF) * bot(x2,muF) ;
    }
    if ((fWp==1) && (c==ppbar)) { return 
      uWf * u(x1,muF) * ubar(x2,muF) + cWf * chm(x1,muF) * chm(x2,muF)
    + dWf * dbar(x1,muF) * d(x2,muF) + sWf * str(x1,muF) * str(x2,muF)
    + bWf * bot(x1,muF) * bot(x2,muF) ;
    }
    if ((fWm==1) && (c==pp)) { return
      uWf * ubar(x1,muF) * ubar(x2,muF) + cWf * chm(x1,muF) * chm(x2,muF)
    + dWf * d(x1,muF) * d(x2,muF) + sWf * str(x1,muF) * str(x2,muF)
    + bWf * bot(x1,muF) * bot(x2,muF) ;
    }
    if ((fWm==1) && (c==ppbar)) { return
      uWf * ubar(x1,muF) * u(x2,muF) + cWf * chm(x1,muF) * chm(x2,muF)
    + dWf * d(x1,muF) * dbar(x2,muF) + sWf * str(x1,muF) * str(x2,muF)
    + bWf * bot(x1,muF) * bot(x2,muF) ;
    }
    return 0.;
  }
  else if ((c==pp) || (c==ppbar)) { return qq_CE1_lumi(x1,x2,muF,c) ; }
  else return 0.;
}
#endif
