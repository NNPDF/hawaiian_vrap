#pragma once

#include "settings.h"

/*================================================================
 The full (NLO, NNLO) DY calculation, with cuts, divides up into 
 3 pieces:
 1) A Born kinematics piece (Born), with only DY in final state, and
  proportional to delta(1-z).
 2a) Two boosted Born kinematics pieces (boost), where 
     y_DY = 1/2 * ln(x1/x2) \pm ln(z),
 and
 2b) A fully smeared distribution (real), for the full DY + parton(s)
 final state.
 Here we artificially smear 2a) and add it to 2b) 
=================================================================*/

/*=================================================================
  The overall q qbar -> DY normalization:
  The "Born-level" partonic cross section in pb is sigma_DY_0.
  It omits the quark charges, which will be supplied elsewhere.
  (Now sigma_DY_0 is factored out front of everything.)
  The "2/M" normalization gives d sigma/dM.                     */

double sigma_DY_0(double M, double alpha){
  return GeV2pb * 2./M * 4.*PI*alpha*alpha/3./N_c /M/M ;
}

// Total prefactor:

double DY_prefactor(double M, double alpha){
  return sigma_DY_0(M,alpha) * M*M/E_CM/E_CM;
}

//================================================================
// The Born-type terms in the integrand.
double Born_integrand(double y){
    double asR = alpha_s(muR)/PI;
    double tau = Q*Q/E_CM/E_CM;
    double x1 = sqrt(tau) * exp(y);
    double x2 = sqrt(tau) * exp(-y);

    if ((x1 >= 1.) || (x2 >= 1.)){ return 0.; }

    double LO_term = 1.0;
    double NLO_term = Born_NLO(muF/Q);
    double NNLO_term = Born_NNLO(Nf,muF/Q,muR/Q);
    double Born_term;

    if (order_flag == 0) {
        Born_term = LO_term;
    }
    if (order_flag == 1) {
        Born_term = LO_term + asR * NLO_term;
    }
    if (order_flag == 2) {
        if (f_NNLO_only == 0) {
            Born_term = asR*asR * NNLO_term;
        }
        else {
            Born_term = LO_term + asR * NLO_term + asR*asR * NNLO_term;
        }
    }

    pdfArray X1,X2;

    LHAComputePdf(x1,muF,X1);	
    LHAComputePdf(x2,muF,X2);

    // Now fill in the pineappl grids
    piner.fill_grid(0, &qqbar_lumi, x1, x2, LO_term);
    if (order_flag > 0) piner.fill_grid(1, &qqbar_lumi, x1, x2, NLO_term);
    if (order_flag > 1) piner.fill_grid(4, &qqbar_lumi, x1, x2, NNLO_term);

    // Calls `qqbar_lumi` probably from `Vlumifns_LHApdf.C`
    // but this can change depending on where the headers are coming from...
    return Born_term * qqbar_lumi(X1,X2,DY,coll);
}


//================================================================
// The boost + real-type terms in the integrand.
// NLO and NNLO terms are separate routines.

double int_NLO(double y, double ys, double z){
  double lumiz1, lumiz2, lumiz, lumi0, qg_lumiz1, gq_lumiz2, boost_ans,
         lumi_ys, lumi_0, lumi_1, gq_lumi_ys, qg_lumi_ys, 
         gq_lumi_1, qg_lumi_0, real_ans;
  double tau = Q*Q/E_CM/E_CM;
// Boost terms in integrand:
  // for ys = 0 boost (x2 is the radiator):
  double x1a = sqrt(tau) * exp(y) ;
  double x2a = sqrt(tau) * exp(-y)/z ;

	pdfArray X1a,X2a;
	LHAComputePdf(x1a,muF,X1a);	
	LHAComputePdf(x2a,muF,X2a);	
	lumiz1 = qqbar_lumi(X1a,X2a,DY,coll); 
	qg_lumiz1 = qg_lumi(X1a,X2a,coll);
  // for ys = 1 boost (x1 is the radiator):
  double x1b = sqrt(tau) * exp(y)/z ;
  double x2b = sqrt(tau) * exp(-y) ;
	pdfArray X1b,X2b;
	LHAComputePdf(x1b,muF,X1b);	
	LHAComputePdf(x2b,muF,X2b);	
  lumiz2 = qqbar_lumi(X1b,X2b,DY,coll); 
  gq_lumiz2 = gq_lumi(X1b,X2b,coll);
lumiz = lumiz1 + lumiz2 ;
  // for z = 1 subtraction (q qbar only):
double x1_z1 = sqrt(tau) * exp(y) ;
  double x2_z1 = sqrt(tau) * exp(-y) ;
	pdfArray X1_z1,X2_z1;
	LHAComputePdf(x1_z1,muF,X1_z1);	
	LHAComputePdf(x2_z1,muF,X2_z1);	
	lumi0 = 2. * qqbar_lumi(X1_z1,X2_z1,DY,coll);
  boost_ans = 
  // qqbar
      NLO_qbarq_boost(z,muF/Q) * lumiz/z
    - NLO_qbarq_boost_soft(z,muF/Q) * lumi0
  // gq + qg 
    + NLO_qg_boost(z,muF/Q) * ( gq_lumiz2 + qg_lumiz1 )/z ;
//
// Real terms in integrand:
  double x1 = sqrt(tau) * exp(y) * sqrt((z+(1.-z)*ys)/z/(1.-(1.-z)*ys)) ;
  double x2 = tau/x1/z ;
	pdfArray X1,X2;
	LHAComputePdf(x1,muF,X1);	
	LHAComputePdf(x2,muF,X2);	
	lumi_ys = qqbar_lumi(X1,X2,DY,coll); 
	qg_lumi_ys = qg_lumi(X1,X2,coll);
	gq_lumi_ys = gq_lumi(X1,X2,coll);
  // for ys = 0 subtraction:
  double x1_0 = sqrt(tau) * exp(y) ;
  double x2_0 = tau/x1_0/z ;
	pdfArray X1_0,X2_0;
	LHAComputePdf(x1_0,muF,X1_0);	
	LHAComputePdf(x2_0,muF,X2_0);	
	lumi_0 = qqbar_lumi(X1_0,X2_0,DY,coll); 
	qg_lumi_0 = qg_lumi(X1_0,X2_0,coll);
  // for ys = 1 subtraction:
  double x1_1 = sqrt(tau) * exp(y)/z ;
  double x2_1 = tau/x1_1/z ;
	pdfArray X1_1,X2_1;
	LHAComputePdf(x1_1,muF,X1_1);	
	LHAComputePdf(x2_1,muF,X2_1);	
	lumi_1 = qqbar_lumi(X1_1,X2_1,DY,coll); 
	gq_lumi_1 = gq_lumi(X1_1,X2_1,coll);
  //
  real_ans = 
      // qqbar
        NLO_qbarq_real(ys,z) * lumi_ys
      - NLO_qbarq_real_soft(ys,z) * lumi_0
      - NLO_qbarq_real_soft(1.-ys,z) * lumi_1
     // qg + gq
      + NLO_qg_real(1.-ys,z) * gq_lumi_ys 
      - NLO_qg_real_soft(1.-ys,z) * gq_lumi_1
      + NLO_qg_real(ys,z) * qg_lumi_ys
      - NLO_qg_real_soft(ys,z) * qg_lumi_0 ;
  real_ans = real_ans/z ;
  return boost_ans + real_ans;
}

double int_NNLO(double y, double ys, double z){
  double lumiz1, lumiz2, lumiz, lumi0, qg_lumiz1, gq_lumiz2, 
    qq_11_lumiz1, qq_22_lumiz2, qq_CE1_lumiz1, qq_CE2_lumiz2, boost_ans,
    lumi_ys, lumi_0, lumi_1, lumi_z1, 
    qbarq_BC_lumi_ys, qbarq_NFf_lumi_ys, qbarq_ax_lumi_ys,
    qg_lumi_ys, qg_lumi_0, gq_lumi_ys, gq_lumi_1, gg_lumi_ys,
    qq_11_lumi_ys, qq_11_lumi_0, qq_22_lumi_ys, qq_22_lumi_1,
    qq_12_lumi_ys, qq_12_ax_lumi_ys, 
    qq_CE1_lumi_ys, qq_CE2_lumi_ys, qq_CF_lumi_ys, 
    qq_CE1_lumi_0, qq_CE2_lumi_1, real_ans;
  if ((ys < 1.0e-08) || (ys > 1.0-1.0e-08) || (z > 1.0-1.0e-08)){ return 0.; }
  double tau = Q*Q/E_CM/E_CM;
  if ((tau * exp(2.*y) >= 1.) || (tau * exp(-2.*y) >= 1.)
                                 || (z <= tau)){ return 0.; } 
// Boost terms in integrand:
  // for ys = 0 boost (x2 is the radiator):
  double x1a = sqrt(tau) * exp(y) ;
  double x2a = sqrt(tau) * exp(-y)/z ;
	pdfArray X1a,X2a;
	LHAComputePdf(x1a,muF,X1a);	
	LHAComputePdf(x2a,muF,X2a);	
	lumiz1 = qqbar_lumi(X1a,X2a,DY,coll); 
	qg_lumiz1 = qg_lumi(X1a,X2a,coll); 
	qq_11_lumiz1 = qq_11_lumi(X1a,X2a,coll);
	qq_CE1_lumiz1 = qq_CE1_lumi(X1a,X2a,coll);
  // for ys = 1 boost (x1 is the radiator):
  double x1b = sqrt(tau) * exp(y)/z ;
  double x2b = sqrt(tau) * exp(-y) ;
	pdfArray X1b,X2b;
	LHAComputePdf(x1b,muF,X1b);	
	LHAComputePdf(x2b,muF,X2b);	
	lumiz2 = qqbar_lumi(X1b,X2b,DY,coll); 
	gq_lumiz2 = gq_lumi(X1b,X2b,coll);
	qq_22_lumiz2 = qq_22_lumi(X1b,X2b,coll); 
	qq_CE2_lumiz2 = qq_CE2_lumi(X1b,X2b,coll);
  lumiz = lumiz1 + lumiz2 ;
  // for z = 1 subtraction (q qbar only):
  double x1_z_1 = sqrt(tau) * exp(y) ;
  double x2_z_1 = sqrt(tau) * exp(-y) ;
	pdfArray X1_z_1,X2_z_1;
	LHAComputePdf(x1_z_1,muF,X1_z_1);	
	LHAComputePdf(x2_z_1,muF,X2_z_1);	
	lumi0 = 2. * qqbar_lumi(X1_z_1,X2_z_1,DY,coll);
  boost_ans = 
     ( qbarq_boost_soft(z,Nf,muF/Q,muR/Q) + qbarq_boost_hard(z,Nf,muF/Q,muR/Q)
       // removed to q_iq_j channel:
     - qq11_boost_tot(z,1.) )  * lumiz/z 
    - qbarq_boost_soft(z,Nf,muF/Q,muR/Q) * lumi0 
    + qg_boost_tot(z,Nf,muF/Q,muR/Q) * ( qg_lumiz1/z + gq_lumiz2/z )  
    + qq11_boost_tot(z,muF/Q) * ( qq_11_lumiz1/z + qq_22_lumiz2/z )
    + qq_CE_boost_tot(z,muF/Q) * ( qq_CE1_lumiz1/z + qq_CE2_lumiz2/z )
    ; 
//
// Real terms in integrand:
  double x1 = sqrt(tau) * exp(y) * sqrt((z+(1.-z)*ys)/z/(1.-(1.-z)*ys)) ;
  double x2 = tau/x1/z ;
	pdfArray X1,X2;
	LHAComputePdf(x1,muF,X1);	
	LHAComputePdf(x2,muF,X2);	
	lumi_ys = qqbar_lumi(X1,X2,DY,coll); 
	qbarq_BC_lumi_ys = qqbar_BC_lumi(X1,X2,coll); // BUGFIX, 11/21/2003 (LD)
	qbarq_NFf_lumi_ys = qqbar_lumi(X1,X2,gluon,coll);
	qbarq_ax_lumi_ys = qqbar_ax_lumi(X1,X2,coll);
	qg_lumi_ys = qg_lumi(X1,X2,coll);
	gq_lumi_ys = gq_lumi(X1,X2,coll);
	gg_lumi_ys = gg_lumi(X1,X2,coll);
	qq_11_lumi_ys = qq_11_lumi(X1,X2,coll);
	qq_22_lumi_ys = qq_22_lumi(X1,X2,coll);
	qq_12_lumi_ys = qq_12_lumi(X1,X2,coll);
	qq_12_ax_lumi_ys = qq_12_ax_lumi(X1,X2,coll);
	qq_CE1_lumi_ys = qq_CE1_lumi(X1,X2,coll);
	qq_CE2_lumi_ys = qq_CE2_lumi(X1,X2,coll);
	qq_CF_lumi_ys = qq_CF_lumi(X1,X2,coll);
  // for ys = 0 subtraction:
  double x1_0 = sqrt(tau) * exp(y) ;
  double x2_0 = tau/x1_0/z ;
	pdfArray X1_0,X2_0;
	LHAComputePdf(x1_0,muF,X1_0);	
	LHAComputePdf(x2_0,muF,X2_0);	
	lumi_0 = qqbar_lumi(X1_0,X2_0,DY,coll); 
  qg_lumi_0 = qg_lumi(X1_0,X2_0,coll); 
  qq_11_lumi_0 = qq_11_lumi(X1_0,X2_0,coll); 
  qq_CE1_lumi_0 = qq_CE1_lumi(X1_0,X2_0,coll);
  // for ys = 1 subtraction:
  double x1_1 = sqrt(tau) * exp(y)/z ;
  double x2_1 = tau/x1_1/z ;
	pdfArray X1_1,X2_1;
	LHAComputePdf(x1_1,muF,X1_1);	
	LHAComputePdf(x2_1,muF,X2_1);	
	lumi_1 = qqbar_lumi(X1_1,X2_1,DY,coll); 
  gq_lumi_1 = gq_lumi(X1_1,X2_1,coll); 
  qq_22_lumi_1 = qq_22_lumi(X1_1,X2_1,coll); 
  qq_CE2_lumi_1 = qq_CE2_lumi(X1_1,X2_1,coll);
  // for z = 1 subtraction (q qbar only):
  double x1_z1 = sqrt(tau) * exp(y) ;
  double x2_z1 = sqrt(tau) * exp(-y) ;
	pdfArray X1_z1,X2_z1;
	LHAComputePdf(x1_z1,muF,X1_z1);	
	LHAComputePdf(x2_z1,muF,X2_z1);	
	lumi_z1 = qqbar_lumi(X1_z1,X2_z1,DY,coll);
  //
  real_ans = 
  // q-\bar{q}:
    // 1) nonsingular terms: now includes unsubtracted soft pieces;
    //   BC terms have been removed, now have own luminosity function
    //   C^2 and C*D terms have been removed to q_i q_j channel;
    //   muF- and muR-dependence now included
      qbarq_real_hard(ys,z,Nf,muF/Q,muR/Q) * lumi_ys
    // 2) y -> 0 and y -> 1 singular subtraction terms:
    - ( qbarq_real_soft_0(ys,z,Nf,muF/Q,muR/Q)
       // removed to q_iq_j channel:
       - qq11_real_soft(ys,z,1.) ) * lumi_0
    - ( qbarq_real_soft_1(ys,z,Nf,muF/Q,muR/Q)
       // removed to q_iq_j channel:
       - qq11_real_soft(1.-ys,z,1.) ) * lumi_1
     // 3) z -> 1 singular + subtraction terms (they alone work):
    - qbarq_real_soft_z1(ys,z,Nf,muF/Q) * lumi_z1 * z
  // q\bar{q} BC terms:
    + qbarq_BC_real_hard(ys,z) * qbarq_BC_lumi_ys
  // q\bar{q} NFf terms:
    + qbarq_NFf_real_hard(ys,z) * qbarq_NFf_lumi_ys
    + qq_AB_ax_real_hard(ys,z) * qbarq_ax_lumi_ys  // extra axial term for Z
  // qg:
    + qg_real_hard(ys,z,Nf,muF/Q,muR/Q) * qg_lumi_ys
    + qg_real_hard(1.-ys,z,Nf,muF/Q,muR/Q) * gq_lumi_ys
    - qg_real_soft(ys,z,Nf,muF/Q,muR/Q) * qg_lumi_0
    - qg_real_soft(1.-ys,z,Nf,muF/Q,muR/Q) * gq_lumi_1
  // gg (only "hard" terms):
    + gg_real_hard(ys,z,muF/Q) * gg_lumi_ys
  // q_i q_j (qbar_j)   and   identical quark terms:
    + qq11_real_hard(ys,z,muF/Q) * qq_11_lumi_ys
    + qq11_real_hard(1.-ys,z,muF/Q) * qq_22_lumi_ys
    - qq11_real_soft(ys,z,muF/Q) * qq_11_lumi_0
    - qq11_real_soft(1.-ys,z,muF/Q) * qq_22_lumi_1
    + qq12_real_hard(ys,z) * qq_12_lumi_ys         // a.k.a. CD_V
    + qq_CD_ax_real_hard(ys,z) * qq_12_ax_lumi_ys  // CD_axial term, for Z
    + qq_CE_real_hard(ys,z) * qq_CE1_lumi_ys
    + qq_CE_real_hard(1.-ys,z) * qq_CE2_lumi_ys
    + qq_CF_real_hard(ys,z) * qq_CF_lumi_ys
    - qq_CE_real_soft(ys,z) * qq_CE1_lumi_0
    - qq_CE_real_soft(1.-ys,z) * qq_CE2_lumi_1
    ;
  real_ans = real_ans/z;
  // TEMPORARY DIAGNOSTICS:
  if (std::fabs(boost_ans + real_ans) >= 0.0) { }
  else { std::cout << std::setw(15) << std::setprecision(15)
                << "  y =  " << y
                << "  ;  ys  " << ys 
                << "  ;  z =  " << z << std::endl; 
           std::cout << " boost_ans =  " << boost_ans
                << " ;   real_ans =  " << real_ans << std::endl;
           std::cout << " qq_CD_ax_real_hard(ys,z) =  " 
                << qq_CD_ax_real_hard(ys,z) 
                << " ; qq_12_ax_lumi_ys = "  << qq_12_ax_lumi_ys << std::endl;
           std::cout << " qq_AB_ax_real_hard(ys,z) =  " 
                << qq_AB_ax_real_hard(ys,z) 
                << " ; qbarq_ax_lumi_ys = " << qbarq_ax_lumi_ys << std::endl;
  }
  return boost_ans + real_ans;
}

/*=================================================================
 For integrating over M, we need slightly different integrands:
 We have to reset the lumi-functions for a different M,
 which also requires passing along exch,alphat 
 (we assume Nf stays constant over the integration range).  */

double Born_integrand_M(double y, double M){
  Q = M;    muR = muRrel*Q;    muF = muFrel*Q;
  double alpha_local = alphat;
  // if you want to run alpha_QED here:
  // double alpha_local = alpha_QED(Q); 
  setV(exchM,M,alpha_local,Nf,1);
  return DY_prefactor(Q,alpha_local) * Born_integrand(y);
}

double int_NLO_M(double y, double ys, double z, double M){
  Q = M;    muR = muRrel*Q;    muF = muFrel*Q;
  double alpha_local = alphat;
  // if you want to run alpha_QED here:
  // double alpha_local = alpha_QED(Q); 
  setV(exchM,M,alpha_local,Nf,1);
  return DY_prefactor(Q,alpha_local) * alpha_s(muR)/PI * int_NLO(y,ys,z);
}

double int_NNLO_M(double y, double ys, double z, double M){
  Q = M;    muR = muRrel*Q;    muF = muFrel*Q;
  double alpha_local = alphat;
  // if you want to run alpha_QED here:
  // double alpha_local = alpha_QED(Q); 
  setV(exchM,M,alpha_local,Nf,1);
  return DY_prefactor(Q,alpha_local) * alpha_s(muR)/PI * alpha_s(muR)/PI 
                                                       * int_NNLO(y,ys,z);
}

//============================================================
class surf: public Surface{
public:
  int n;
  surf(int m): n(m){}

  double surface(DVector & x){
// Integrands for evaluation at fixed y:
    if (n==1){  // evaluate NLO integrand at fixed y
      return int_NLO(y, x[1], x[2]) ; }
    if (n==2){  // evaluate NNLO integrand at fixed y
      return int_NNLO(y, x[1], x[2]) ; }
// Integrands for integration over a range in y:
   // introduce change of variables to do integral over y:
   // let xi = (1+tanh(y))/2,  y = ln(x/(1-x))/2.
   // Then let xi = xi_l * (1-xx) + xi_u * xx,
   // xx = (xi-xi_l)/(xi_u-xi_l),   so that xx runs from 0 to 1.
    double xi = xi_l * (1-x[1]) + xi_u * x[1] ;
    double ytemp = 0.5*log(xi/(1-xi));
    double Jy = 0.5*(xi_u-xi_l)/xi/(1.-xi);
    if (n==3) {  // integrate Born_integrand over y
      return Jy * Born_integrand(ytemp) ; }
    if (n==4) {  // integrate NLO integrand over y
      return Jy * int_NLO(ytemp, x[2], x[3]) ; }
    if (n==5) {  // integrate NNLO integrand over y
      return Jy * int_NNLO(ytemp, x[2], x[3]) ; }
// Integrands for integration over a range in M:
    // mapping of M range containing Z pole:
     double Gamma_Z = 2.4952; 
     double u_l = atan((Ml*Ml-m_Z*m_Z)/Gamma_Z/m_Z);
     double u_u = atan((Mu*Mu-m_Z*m_Z)/Gamma_Z/m_Z);
     double u = u_l + (u_u-u_l) * x[1];
     double M = sqrt(m_Z*m_Z + Gamma_Z*m_Z*tan(u));
     double JM = (u_u-u_l) * Gamma_Z*m_Z*(1.+tan(u)*tan(u))/2./M;
   // uniform mapping of M range:
    //    double M = Ml + (Mu-Ml) * x[1];   double JM = Mu-Ml; 
    if (n==6) { // integrate Born integrand over M:
      return JM * Born_integrand_M(y,M); }
    if (n==7) { // integrate NLO integrand over M:
      return JM * int_NLO_M(y, x[2], x[3], M) ; }
    if (n==8) { // integrate NNLO integrand over M:
      return JM * int_NNLO_M(y, x[2], x[3], M) ; }
    else return 0;
  }

};

// "sqint" runs a VEGAS integration for function number n_f,
// over a n-dimensional unit square (hypercube), using surface "surf".
// It returns "integral, sd, chi2" as a DVector:

DVector sqint(int n_f, int n_dim, int n_call, int n_adapt, int n_run, 
              int ranseed){
  double integral, sd, chi2; 
  DVector temp_out(0,2);
  DVector x1(1,n_dim);  DVector x2(1,n_dim);
  for (int i=1; i <= n_dim; i++) {
     x1[i] = 0.;   x2[i] = 1.;
  }
  VegasGrid V(x1,x2,n_call,ranseed);  // last argument is "ranseed"
  surf Sfc(0);
  V.reset(n_call);
  Sfc.n = n_f;
  if (f_quiet==0){ std::cout << "adapting grid"<<std::endl; }
  integral = Sfc.vegas(V,n_adapt,sd,chi2);
  if (f_quiet==0){
  std::cout << std::setw(8) << std::setprecision(8) << integral 
       << "  pm  " << sd << "  chi2 = " << chi2 << std::endl;
  std::cout << "evaluating integral" <<std::endl;
  }  
  integral = Sfc.vegasint(V,n_run,sd,chi2);
  if (f_quiet==0){  
  std::cout << std::setw(8) << std::setprecision(8) << integral 
       << "  pm  " << sd << "  chi2 = " << chi2 << std::endl;
  std::cout << std::endl;
  }
  temp_out.fill(3,integral,sd,chi2);
  return temp_out;
}

// ========== BASIC INTEGRATION ROUTINES FOR \int dy and fixed y =======

// Routine to integrate over y, from y_lower to y_upper,
// or actually, the related variables xi_l to xi_u:

DVector int_y(double xil, double xiu){
  DVector result(0,1);
  DVector Born_int(0,2);   DVector i_NLO(0,2);   DVector i_NNLO(0,2);
  xi_l = xil;   xi_u = xiu;
  double total_int, total_int_error;
  double prefactor = DY_prefactor(Q,alphat);
  double asR = alpha_s(muR)/PI;
  std::cout  << std::endl << "Now computing integral over rapidity" << std::endl << std::endl;
// LO case
  if (order_flag == 0) {
    if (std::fabs(Born_integrand(0.2)) < 1.0e-16){
      std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      result.fill(2,0.,0.); 
    }
    else {
      std::cout << "   Now integrating Born terms " << std::endl;
      Born_int = sqint(3,1,1000,10,20,ranseed);
      total_int = prefactor * Born_int[0];
      total_int_error = prefactor * Born_int[1];
      std::cout << " LO integral =  " << std::setw(8) << std::setprecision(8)
           << prefactor*Born_int[0] << "  pm  " << prefactor*Born_int[1]
           <<  ";    chi^2 =  " << Born_int[2] << std::endl << std::endl;
    }
  }
// NLO case
  if (order_flag == 1) {
    if (std::fabs(Born_integrand(0.2)) < 1.0e-16){
      std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      Born_int.fill(3,0.,0.,0.);
    }
    else {
     std::cout << "   Now integrating Born terms " << std::endl;
     Born_int = sqint(3,1,1000,10,20,ranseed);
     std::cout << " Born-type integral =  " << std::setw(8) << std::setprecision(8)
          << prefactor*Born_int[0] << "  pm  " << prefactor*Born_int[1]
          <<  ";    chi^2 =  " << Born_int[2] << std::endl << std::endl;
    }
    if (std::fabs(int_NLO(.03,.43,.93)) < 1.0e-16) {
      std::cout << "Setting non-Born NLO terms to 0, since integrand ~ 0 " << std::endl;
      i_NLO.fill(3,0.,0.,0.);
    }
    else {
      std::cout << " Now integrating NLO terms " << std::endl;
      i_NLO = sqint(4,3,10000,10,20,ranseed);
      std::cout << " NLO non-Born integral =  " << std::setw(8) << std::setprecision(8)
          << prefactor * asR * i_NLO[0] 
          << "  pm  " << prefactor * asR * i_NLO[1]
          <<  ";    chi^2 =  " << i_NLO[2] << std::endl << std::endl;
    }
    total_int = prefactor * ( Born_int[0] + asR * i_NLO[0] );
    total_int_error = prefactor * sqrt( Born_int[1]*Born_int[1] 
                                      + asR*asR * i_NLO[1]*i_NLO[1] ) ;
    std::cout << " NLO integral =  " << std::setw(8) << std::setprecision(8)
         << total_int << "  pm  " << total_int_error << std::endl << std::endl;
  }
// NNLO case
  if (order_flag == 2) {
    if (std::fabs(Born_integrand(0.2)) < 1.0e-16){
      std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      Born_int.fill(3,0.,0.,0.);
    }
    else {
    std::cout << " Now integrating Born terms " << std::endl;
    Born_int = sqint(3,1,1000,10,20,ranseed);
    std::cout << " Born-type integral =  " << std::setw(8) << std::setprecision(8)
         << prefactor*Born_int[0] << "  pm  " << prefactor*Born_int[1]
         <<  ";    chi^2 =  " << Born_int[2] << std::endl << std::endl;
    }
    if (f_NNLO_only == 1){
      if (std::fabs(int_NLO(.03,.43,.93)) < 1.0e-16) {
        std::cout << "Setting NLO non-Born terms to 0, since integrand ~ 0 " << std::endl;
        i_NLO.fill(3,0.,0.,0.);
      }
      else {
      std::cout << " Now integrating NLO terms " << std::endl;
      i_NLO = sqint(4,3,5000,10,20,ranseed);
      std::cout << " NLO non-Born integral =  " << std::setw(8) << std::setprecision(8)
           << prefactor * asR * i_NLO[0] 
           << "  pm  " << prefactor * asR * i_NLO[1]
           <<  ";    chi^2 =  " << i_NLO[2] << std::endl << std::endl;
      }
    }
    std::cout << "   Now integrating NNLO terms " << std::endl;
    i_NNLO = sqint(5,3,5000,10,30,ranseed);
    std::cout << " NNLO non-Born integral =  " << std::setw(8) << std::setprecision(8)
         << prefactor * asR*asR * i_NNLO[0] 
         << "  pm  " << prefactor * asR*asR * i_NNLO[1]
         <<  ";    chi^2 =  " << i_NNLO[2] << std::endl 
         << std::endl;
    if (f_NNLO_only == 0) {
      total_int = prefactor * ( Born_int[0] + asR*asR * i_NNLO[0] );
      total_int_error = prefactor * sqrt( Born_int[1]*Born_int[1] 
                                   + asR*asR*asR*asR * i_NNLO[1]*i_NNLO[1] ) ;
      std::cout << " NNLO (ONLY) integral =  " << std::setw(8) << std::setprecision(8)
         << total_int << "  pm  " << total_int_error  
         << std::endl << std::endl;
    }
    else {
    total_int = prefactor * ( Born_int[0] 
			    + asR * i_NLO[0] + asR*asR * i_NNLO[0] );
    total_int_error = prefactor * sqrt( Born_int[1]*Born_int[1] 
                            + asR*asR * i_NLO[1]*i_NLO[1]
                   + asR*asR*asR*asR * i_NNLO[1]*i_NNLO[1] ) ;
    std::cout << " NNLO integral =  " << std::setw(8) << std::setprecision(8)
       << total_int << "  pm  " << total_int_error  
       << std::endl << std::endl;
    }
  }
 result.fill(2,total_int,total_int_error);  
 return result;
}

// Routine computes d sigma/dM/dy, plus errors:
DVector rap_y(){
    double total, total_error;
    double prefactor = DY_prefactor(Q,alphat);
    std::cout << "Prefactor = " << prefactor << std::endl;
    double asR = alpha_s(muR)/PI;
    DVector result(0,1);   DVector i_NLO(0,2);  DVector i_NNLO(0,2);
    if (f_quiet==0) {
        std::cout << std::setw(9) << std::setprecision(9) << " working on    y =  " << y << std::endl; 
    }

    // Set the prefactor for the pineappl grid
    // this needs to be done at this stage since it is not included in the weights
    piner.set_prefactor(prefactor);

    // Note that Born_ans is different depending on the `order_flag`
    double Born_ans = Born_integrand(y);

    // LO case
    if (order_flag == 0) { 
        std::cout << "Staring LO calculation: " << std::endl;
        total = prefactor * Born_ans;
        total_error = 0.;
        std::cout << "LO   = " << total << std::endl;
    }

    // NLO case
    else if (order_flag == 1) {
        std::cout << "LO   = " << prefactor * Born_ans << std::endl;

        if (std::fabs(int_NLO(.03,.43,.93)) < 1.0e-16) {
            if (f_quiet==0) std::cout << "Setting NLO non-Born terms to 0, since integrand ~ 0 " << std::endl;

            i_NLO.fill(3,0.,0.,0.);
        } else {
            i_NLO = sqint(1,2,5000,10,20,ranseed); 
        }

        total = prefactor * ( Born_ans + asR * i_NLO[0] );
        total_error = prefactor * asR * i_NLO[1] ;

        std::cout << "NLO  = " << prefactor * (Born_ans + asR * i_NLO[0]) << std::endl;
    }
    // NNLO case
    else {   //  order_flag = 2:
        std::cout << "LO   = " << prefactor * Born_ans << std::endl;
        if (f_NNLO_only == 1) {
            if (f_quiet==0) {
                std::cout << " computing NLO integral for y = " << y << std::endl; 
            }
            if (std::fabs(int_NLO(.03,.43,.93)) < 1.0e-16) {
                if (f_quiet==0) {
                    std::cout << "Setting NLO non-Born terms to 0, since integrand ~ 0 " << std::endl;
                }
                i_NLO.fill(3,0.,0.,0.);
            }
            else {
                i_NLO = sqint(1,2,5000,10,20,ranseed); 
            }
        }

        std::cout << "NLO  = " << prefactor * (Born_ans + asR * i_NLO[0]) << std::endl;
        if (f_quiet==0) {
            std::cout << " computing NNLO integral for y = " << y << std::endl; 
        }

        i_NNLO = sqint(2,2,20000,3,10,ranseed);

        if (f_NNLO_only == 0) {
            total = prefactor * ( Born_ans + asR*asR * i_NNLO[0] );
            total_error = prefactor * asR*asR * i_NNLO[1] ;
        } else {
            total = prefactor * ( Born_ans + asR * i_NLO[0] + asR*asR * i_NNLO[0] );
            total_error = prefactor * asR * sqrt( i_NLO[1]*i_NLO[1]
                                                + asR*asR * i_NNLO[1]*i_NNLO[1] );
            std::cout << "NNLO = " << total << std::endl;
        }
    }

    if (f_quiet==0) {
        std::cout << std::setw(8) << std::setprecision(8) <<  " y = " << y 
            << ";   d^2sigma/dM/dy =  " << total
            << "   pm  " << total_error << std::endl;
    }

    result.fill(2,total,total_error);
    return result;
}

// ============ SCAN IN RAPIDITY y =========================

// Make table of K factors for fixed-target Drell-Yan, M vs. y.
// There are
//       M_steps evenly-spaced values of M, from M_min to M_max,
// vs.   y_steps evenly-spaced values of y, from   0   to y_max.
// For J. Stirling.

void makeKtable(double y_max, int y_steps, 
                double M_min, double M_max, int M_steps){
  double yMAX;   // the physical limit for y
  f_quiet = 1;  // shut off output
  std::cout << std::endl
       << "--------------------------------------------------------" << std::endl; 
  std::cout << std::setw(8) << std::setprecision(8) << "Performing mass/rapidity scan:  " 
       << M_min << "  < M < " << M_max << "  ;   0 < y < " << y_max 
       << std::endl; 
 std::cout << "--------------------------------------------------------" << std::endl; 
 std::cout << "number of mass steps =  " << M_steps 
      << "; number of rapidity steps =  " << y_steps << std::endl << std::endl; 
  double y_step_size = y_max/y_steps;
  double M_step_size = (M_max-M_min)/M_steps;
// print header
  std::cout << " y =         " ;
  for (int j = 0; j <= y_steps; j++) {
    y = 0. + y_step_size * j;
    std::cout << y << "           " ;
  }
  std::cout << std::endl;
  std::cout << "---------------------------------------------" 
       << "---------------------------------------------" << std::endl;
//
  for (int i = 0; i <= M_steps; i++) {
    Q = M_min + M_step_size * i;
    std::cout << " M =  " << Q << " :    " ;
    yMAX = log(E_CM/Q);
    muF = Q;   muR = muF;  // ren. & fact. scales equal to Q here
    for (int j = 0; j <= y_steps; j++) {
      y = 0. + y_step_size * j;
      if (y < yMAX){
// LO MRST:
       order_flag = 0;  alpha_s_Z = 0.130;  
	LHApdfInit(mrst,1);  
       DVector LO_rap_y = rap_y();
// NLO MRST:
       order_flag = 1;  alpha_s_Z = 0.119;  
	LHApdfInit(mrst,2);  

       DVector NLO_rap_y = rap_y();
       std::cout << NLO_rap_y[0]/LO_rap_y[0] << "   " ;
      }
      else { std::cout <<  0.0 << "   " ; }
    }
  std::cout << std::endl;
  }
  std::cout << "------------------------------------------------------- " << std::endl;
}

/* Make table of "C" factors for fixed-target Drell-Yan, M vs. y.
   There are
         M_steps evenly-spaced values of M, from M_min to M_max,
   vs.   y_steps evenly-spaced values of y, from   0   to y_max.
   "C" is defined by J. Stirling as
   sigma(NLO) = sigma (LO, but using NLO pdfs) 
                  * ( 1 + alphas(M)/pi C(M,y) )
			 Now we include errors too.   */
                                                                                     
void makeCtable(double y_max, int y_steps,
                double M_min, double M_max, int M_steps){
  double yMAX;   // the physical limit for y
  f_quiet = 1;  // shut off output
  std::cout << std::endl
       << "--------------------------------------------------------"
       << std::endl << "Running `makeCtable'   " << std::endl;
  std::cout << std::setw(8) << std::setprecision(8) << "Performing mass/rapidity scan:  "
       << M_min << "  < M < " << M_max << "  ;   0 < y < " << y_max
       << std::endl;
  std::cout << "--------------------------------------------------------" <<
                std::endl;
  std::cout << "number of mass steps =  " << M_steps
       << "; number of rapidity steps =  " << y_steps << std::endl << std::endl;
  double y_step_size = y_max/y_steps;
  double M_step_size = (M_max-M_min)/M_steps;
  std::cout <<std::endl<< "=========   Table of C(M,y) values  ===========" <<std::endl<<std::endl;
  for (int i = 0; i <= M_steps; i++) {
    Q = M_min + M_step_size * i;
    std::cout << " M = " << Q << " :" << std::endl 
         << "y    |        C       (   error in C    "  << std::endl ;
    yMAX = log(E_CM/Q);
    muF = Q;   muR = muF;  // ren. & fact. scales equal to Q here
    for (int j = 0; j <= y_steps; j++) {
      y = 0. + y_step_size * j;
      if (y < yMAX){
	// LO, but using NLO [*LO* MRST, NLO Alek.] pdfs:
	//	order_flag = 0;  alpha_s_Z = 0.130;  pdf_init(mrst,2);
        order_flag = 0; 
	LHApdfInit(alekhin,1); 
  alpha_s_Z = Alekhin_alpha_s(m_Z); 
        compute_all();
	DVector LO_rap_y = rap_y();
	// NLO, using SAME NLO [*LO* MRST, NLO Alek.] pdfs::
	//	order_flag = 1;  alpha_s_Z = 0.119; LHApdfInit(mrst,2);
        order_flag = 1; 
	LHApdfInit(alekhin,1); 

	alpha_s_Z = Alekhin_alpha_s(m_Z); 
	// here we can compute only the qqbar, or only the qg, component:
        compute_qqbar();
	DVector NLO_rap_y = rap_y();
        // for qqbar:
	double Cvalue = ( NLO_rap_y[0]/LO_rap_y[0] - 1. )/ (alpha_s(Q)/PI); 
	// For qg use compute_qg() instead of compute_all(), and then this formula:
       	// double Cvalue = NLO_rap_y[0]/LO_rap_y[0]/ (alpha_s(Q)/PI);
	// 
	// error should be same regardless:
        double Cerror = sqrt( NLO_rap_y[0]*NLO_rap_y[0]
                        /LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]
                              * LO_rap_y[1]*LO_rap_y[1]
	       	  + 1./LO_rap_y[0]/LO_rap_y[0] * NLO_rap_y[1]*NLO_rap_y[1] )
                         / (alpha_s(Q)/PI) ;
        std::cout << y << "  |    " << Cvalue  << "   (   " << Cerror << std::endl ;
      }
      else { std::cout << y << "  |    " << 0.0  << "             (   " << 0.0 << std::endl ; }
    }
    std::cout << std::endl<< std::endl ;
  }
  std::cout << "------------------------------------------------------- " << std::endl;
}
  
/* Make table of "D" and "E" factors for fixed-target Drell-Yan, M vs. y.
   There are
         M_steps evenly-spaced values of M, from M_min to M_max,
   vs.   y_steps evenly-spaced values of y, from   0   to y_max.
   "D" and "E" are defined as
   sigma(NNLO) = sigma (LO, but using NNLO pdfs) 
                   * ( 1 + alphas(M)/pi D(M,y) 
                         + (alphas(M)/pi)^2 * E(M,y) )
			 Now we include errors too.   */

void makeDEtable(double y_max, int y_steps,
                double M_min, double M_max, int M_steps){
  double yMAX;   // the physical limit for y
  f_quiet = 1;  // shut off output
  compute_all();
  std::cout << std::endl
       << "--------------------------------------------------------" 
       << std::endl << "Running `makeDEtable'   " << std::endl;
  std::cout << std::setw(8) << std::setprecision(8) << "Performing mass/rapidity scan:  "
       << M_min << "  < M < " << M_max << "  ;   0 < y < " << y_max
       << std::endl;
  std::cout << "--------------------------------------------------------" <<
                std::endl;
  std::cout << "number of mass steps =  " << M_steps
       << "; number of rapidity steps =  " << y_steps << std::endl << std::endl;
  double y_step_size = y_max/y_steps;
  double M_step_size = (M_max-M_min)/M_steps;
  std::cout <<std::endl<< "=========   Table of D(M,y) values  ===========" <<std::endl<<std::endl;
  for (int i = 0; i <= M_steps; i++) {
    Q = M_min + M_step_size * i;
    std::cout << " M = " << Q << " :" << std::endl 
         << "y    |        D       (   error in D    "  << std::endl ;
    yMAX = log(E_CM/Q);
    muF = Q;   muR = muF;  // ren. & fact. scales equal to Q here
    for (int j = 0; j <= y_steps; j++) {
      y = 0. + y_step_size * j;
      if (y < yMAX){
	// LO, but using MRST NNLO average pdfs:
        order_flag = 0;  alpha_s_Z = 0.1155;  
	 LHApdfInit(mrst,6);
        // compute_all();
	DVector LO_rap_y = rap_y();
	// NLO, using SAME NNLO pdfs:
        order_flag = 1;  alpha_s_Z = 0.1155;  
	 LHApdfInit(mrst,6);
	// We can also compute only the qqbar, or only the qg, component.
        // For qqbar, use compute_qqbar() and same Dvalue formula.
        // compute_all();
	DVector NLO_rap_y = rap_y();
	double Dvalue = ( NLO_rap_y[0]/LO_rap_y[0] - 1. )/ (alpha_s(Q)/PI);
	// For qg use compute_qg() instead of compute_all(), and then this formula:
       	// double Dvalue = NLO_rap_y[0]/LO_rap_y[0]/ (alpha_s(Q)/PI);
	// 
	// error should be same regardless:
        double Derror = sqrt( NLO_rap_y[0]*NLO_rap_y[0]
                        /LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]
                              * LO_rap_y[1]*LO_rap_y[1]
	       	  + 1./LO_rap_y[0]/LO_rap_y[0] * NLO_rap_y[1]*NLO_rap_y[1] )
                         / (alpha_s(Q)/PI) ;
        std::cout << y << "  |    " << Dvalue  << "   (   " << Derror << std::endl ;
      }
      else { std::cout << y << "  |    " << 0.0  << "             (   " << 0.0 << std::endl ; }
    }
    std::cout << std::endl<< std::endl ;
  }
// print E table header
  std::cout <<std::endl<< "=========   Table of E(M,y) values  ===========" <<std::endl<<std::endl;
  for (int i = 0; i <= M_steps; i++) {
    Q = M_min + M_step_size * i;
    std::cout << " M = " << Q << " :" << std::endl 
         << "y    |        E       (   error in E    "  << std::endl ;
    yMAX = log(E_CM/Q);
    muF = Q;   muR = muF;  // ren. & fact. scales equal to Q here
    for (int j = 0; j <= y_steps; j++) {
      y = 0. + y_step_size * j;
      if (y < yMAX){
	// LO, but using MRST NNLO average pdfs:
        order_flag = 0;  alpha_s_Z = 0.1155;  
	 LHApdfInit(mrst,6);
        // compute_all();
	DVector LO_rap_y = rap_y();
	// NNLO, using SAME NNLO pdfs:
        order_flag = 2;  alpha_s_Z = 0.1155;  
	 LHApdfInit(mrst,6);
	// We can also compute only the qqbar, or only the qg, component.
        // For qqbar, use compute_qqbar() and same Dvalue formula.
        f_NNLO_only = 0;
        // compute_all();
	DVector NNLO_rap_y = rap_y();
	double Evalue = NNLO_rap_y[0]/LO_rap_y[0]/(alpha_s(Q)/PI)/(alpha_s(Q)/PI); 
        double Eerror = sqrt( NNLO_rap_y[0]*NNLO_rap_y[0]
                        /LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]/LO_rap_y[0]
                              * LO_rap_y[1]*LO_rap_y[1]
	       	  + 1./LO_rap_y[0]/LO_rap_y[0] * NNLO_rap_y[1]*NNLO_rap_y[1] )
                         /(alpha_s(Q)/PI)/(alpha_s(Q)/PI) ;
        std::cout << y << "  |    " << Evalue  << "   (   " << Eerror << std::endl ;
      }
      else { std::cout << y << "  |    " << 0.0  << "             (   " << 0.0 << std::endl ; }
    }
  }
  std::cout << "----------------------------------------------------- " << std::endl;
}
  
// Routine uses "rap_y" to compute un-normalized d sigma/dM/dy, plus errors,
// in the SYMMETRIC case, just scans from y = 0 to y = ymax,
// with a number of points equal to n_points:
// ymax is computed from E_CM and Q; then we rescale to omit the last point.
// For n_points = 19, this agrees with Frank's choice of points.
// forwrev = 1  -> print output from -ymax to +ymax
// forwrev = -1  -> print output from +ymax to -ymax
// output_format = 0  -> print for topdraw:  y,    d sigma/dy  ( error
// output_format = 1  -> just list d sigma/dy values.

DMatrix sym_scan_rap_y(int n_points, int forwrev, int output_format){
  double ymax = log(E_CM/Q) * n_points/(n_points+1.);  
 std::cout << std::endl
      << "--------------------------------------------------------" << std::endl; 
 std::cout << "Performing rapidity scan - symmetric case (0 < y < ymax)" << std::endl; 
 std::cout << "--------------------------------------------------------" << std::endl; 
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
  double step_size = ymax/n_points;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
  // y = -ymax + step_size * j;
     // full result:
     DVector temp_rap_y = rap_y();
     resultMatrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
   if (output_format==0) {
// print results for topdraw:
    std::cout << std::endl
         << "   y               |    d^2sigma/dM/dy   |   (  error   " << std::endl;
    std::cout << "------------------------------------------------------- " << std::endl;
    for (int j = n_points; j > 0; j--) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << -forwrev*resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl;
    }
    for (int j = 0; j <= n_points; j++) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << forwrev*resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl;
    }
    std::cout << std::endl ; }
   if (output_format==1) {
     std::cout << std::endl;
    for (int j = 0; j <= n_points; j++) {
      std::cout << std::setw(8) << std::setprecision(8) << resultMatrix[j][1] << "," ;
	}
    std::cout << std::endl << std::endl; }
//
 return resultMatrix;
}

// Same as "sym_scan_rap_y", but for the ASYMMETRIC case:
// scans from y = -ymax to y = ymax.
// For n_points = 38, this agrees with Frank's choice of points.

DMatrix asym_scan_rap_y(int n_points, int forwrev, int output_format){
  double ymax = log(E_CM/Q) * n_points/(n_points+2.);  
 std::cout << std::endl
      << "-------------------------------------------------------------" 
      << std::endl; 
 std::cout << "Performing rapidity scan - asymmetric case (-ymax < y < ymax)" 
      << std::endl; 
 std::cout << "-------------------------------------------------------------"
      << std::endl; 
  double step_size = 2.*ymax/n_points;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     y = - ymax + step_size * j;
     // full result:
     DVector temp_rap_y = rap_y();
     resultMatrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
  }
   if (output_format==0) {
// print results for topdraw:
   std::cout << std::endl
        << "   y               |    d^2sigma/dM/dy   |   (  error   " << std::endl;
   std::cout << "------------------------------------------------------- " << std::endl;
   if (forwrev==1) {
     for (int j = 0; j <= n_points; j++) {
       std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
            << "       " << resultMatrix[j][1] << " "
            << "   (   " << resultMatrix[j][2] << std::endl;
     }
   }
   if (forwrev==-1) {
     for (int j = n_points; j >=0; j--) {
       std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
            << "       " << resultMatrix[j][1] << " "
            << "   (   " << resultMatrix[j][2] << std::endl;
     }
   }
   std::cout << std::endl ; }
   if (output_format==1) {
     for (int j = 0; j <= n_points; j++) {
       std::cout << std::setw(8) << std::setprecision(8) << resultMatrix[j][1] << "," ;
	}
    std::cout << std::endl << std::endl; }
//
 return resultMatrix;
}

// Same as "sym_scan_rap_y", but does all orders -- LO, NLO, NNLO --
// for muR = muF = mu1*Q, and for muR = muF = mu2*Q,
// using MRST or Alekhin distributions.
// Prints results in format convenient for topdraw.
// pdftype must be set to either mrst or alekhin.
// For n_points = 19, agrees with Frank's choice of points.

void sym_scan_all(pdf_type pdftype, int n_points, double mu1, double mu2){
 double ymax = log(E_CM/Q) * n_points/(n_points+1.);  
 std::cout << std::endl
      << "-----------------------------------------------------------" << std::endl;
 std::cout << "LO, NLO, NNLO rapidity scan - symmetric case (0 < y < ymax)" << std::endl;
 std::cout << "-----------------------------------------------------------" << std::endl;
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
 std::cout << " muF_lower = muR_lower = " << mu1 << " * Q " << std::endl;
 std::cout << " muF_upper = muR_upper = " << mu2 << " * Q " << std::endl << std::endl;
 std::cout << " pdf type = " << pdftype << " (2 = mrst, 3 = alekhin) " 
      << std::endl << std::endl;
  double step_size = ymax/n_points; 
  f_NNLO_only = 1;   compute_all();
//
// Do both LO curves first:
  order_flag = 0;  
  if (pdftype==mrst) { alpha_s_Z = 0.130;  
	 LHApdfInit(mrst,1);
  }
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);   
	alpha_s_Z = Alekhin_alpha_s(m_Z);
  }
  DMatrix LO1Matrix(0,n_points,0,2);
  DMatrix LO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     LO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     LO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// Next do both NLO curves:
  order_flag = 1;  
  if (pdftype==mrst) { 
	alpha_s_Z = 0.119; 
	LHApdfInit(mrst,2);
	  
	}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
      
	alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NLO1Matrix(0,n_points,0,2);
  DMatrix NLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     NLO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     NLO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// Finally do both NNLO curves:
  order_flag = 2;  
  if (pdftype==mrst) { alpha_s_Z = 0.1155;  
	LHApdfInit(mrst,6);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
	alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NNLO1Matrix(0,n_points,0,2);
  DMatrix NNLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     NNLO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     NNLO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
//
// print results: 
// print LO:
  std::cout << std::endl << "set color blue" << std::endl << std::endl;
   std::cout << "title 0.05 10.5 data 'LO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -LO2Matrix[j][0]
            << "       " << LO2Matrix[j][1] << " "
            << "   (   " << LO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO2Matrix[j][0]
            << "       " << LO2Matrix[j][1] << " "
            << "   (   " << LO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO1Matrix[j][0]
            << "       " << LO1Matrix[j][1] << " "
            << "   (   " << LO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -LO1Matrix[j][0]
            << "       " << LO1Matrix[j][1] << " "
            << "   (   " << LO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// print NLO:
  std::cout << std::endl << "set color green" << std::endl << std::endl;
   std::cout << "title 0.05 20.5 data 'NLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NLO2Matrix[j][0]
            << "       " << NLO2Matrix[j][1] << " "
            << "   (   " << NLO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO2Matrix[j][0]
            << "       " << NLO2Matrix[j][1] << " "
            << "   (   " << NLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO1Matrix[j][0]
            << "       " << NLO1Matrix[j][1] << " "
            << "   (   " << NLO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NLO1Matrix[j][0]
            << "       " << NLO1Matrix[j][1] << " "
            << "   (   " << NLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// print NNLO:
  std::cout << std::endl << "set color red" << std::endl << std::endl;
   std::cout << "title 0.05 25.5 data 'NNLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NNLO2Matrix[j][0]
            << "       " << NNLO2Matrix[j][1] << " "
            << "   (   " << NNLO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO2Matrix[j][0]
            << "       " << NNLO2Matrix[j][1] << " "
            << "   (   " << NNLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO1Matrix[j][0]
            << "       " << NNLO1Matrix[j][1] << " "
            << "   (   " << NNLO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NNLO1Matrix[j][0]
            << "       " << NNLO1Matrix[j][1] << " "
            << "   (   " << NNLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
}

// Same as "sym_scan_all", but does W+/W- ratio,
// all orders -- LO, NLO, NNLO --
// for muR = muF = mu1*Q, and for muR = muF = mu2*Q,
// using MRST or Alekhin distributions.
// Prints results in format convenient for topdraw.
// pdftype must be set to either mrst or alekhin.
// For n_points = 19, agrees with Frank's choice of points.

void Wpm_ratio_scan_all(pdf_type pdftype, int n_points, double mu1, double mu2){
 double ymax = log(E_CM/Q) * n_points/(n_points+1.);  
 double temp_ratio, temp_err;
 Q = m_W;  // assumes we are on-shell...
 f_quiet = 1;  // turn off most output
 std::cout << std::endl
      << "-----------------------------------------------------------" << std::endl;
 std::cout << "LO, NLO, NNLO W+/W- rap. scan - sym. case (0 < y < ymax)" << std::endl;
 std::cout << "-----------------------------------------------------------" << std::endl;
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
 std::cout << " muF_lower = muR_lower = " << mu1 << " * Q " << std::endl;
 std::cout << " muF_upper = muR_upper = " << mu2 << " * Q " << std::endl << std::endl;
 std::cout << " pdf type = " << pdftype << " (2 = mrst, 3 = alekhin) " 
      << std::endl << std::endl;
  double step_size = ymax/n_points; 
  f_NNLO_only = 1;   compute_all();
//
// Do both LO curves first:
  order_flag = 0;  
  if (pdftype==mrst) { alpha_s_Z = 0.130; 
	LHApdfInit(mrst,1);
  
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
	alpha_s_Z = Alekhin_alpha_s(m_Z);
  }
  DMatrix LO1Matrix(0,n_points,0,2);
  DMatrix LO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     DVector plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     DVector minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     LO1Matrix.fillRow(j,3,y,temp_ratio,temp_err);
//  
     muF = mu2*Q;   muR = mu2*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     LO2Matrix.fillRow(j,3,y,temp_ratio,temp_err);
   }
// print LO:
  std::cout << std::endl << "set color blue" << std::endl << std::endl;
   std::cout << "title 0.05 10.5 data 'LO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -LO2Matrix[j][0]
            << "       " << LO2Matrix[j][1] << " "
            << "   (   " << LO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO2Matrix[j][0]
            << "       " << LO2Matrix[j][1] << " "
            << "   (   " << LO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO1Matrix[j][0]
            << "       " << LO1Matrix[j][1] << " "
            << "   (   " << LO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -LO1Matrix[j][0]
            << "       " << LO1Matrix[j][1] << " "
            << "   (   " << LO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// Next do both NLO curves:
  order_flag = 1;  
  if (pdftype==mrst) { alpha_s_Z = 0.119; 
	LHApdfInit(mrst,2);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
	alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NLO1Matrix(0,n_points,0,2);
  DMatrix NLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     DVector plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     DVector minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NLO1Matrix.fillRow(j,3,y,temp_ratio,temp_err);
//  
     muF = mu2*Q;   muR = mu2*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NLO2Matrix.fillRow(j,3,y,temp_ratio,temp_err);
   }
// print NLO:
  std::cout << std::endl << "set color green" << std::endl << std::endl;
   std::cout << "title 0.05 20.5 data 'NLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NLO2Matrix[j][0]
            << "       " << NLO2Matrix[j][1] << " "
            << "   (   " << NLO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO2Matrix[j][0]
            << "       " << NLO2Matrix[j][1] << " "
            << "   (   " << NLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO1Matrix[j][0]
            << "       " << NLO1Matrix[j][1] << " "
            << "   (   " << NLO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NLO1Matrix[j][0]
            << "       " << NLO1Matrix[j][1] << " "
            << "   (   " << NLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// Finally do both NNLO curves:
  order_flag = 2;  
  if (pdftype==mrst) { alpha_s_Z = 0.1155;  
	LHApdfInit(mrst,6);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NNLO1Matrix(0,n_points,0,2);
  DMatrix NNLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     DVector plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     DVector minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NNLO1Matrix.fillRow(j,3,y,temp_ratio,temp_err);
//  
     muF = mu2*Q;   muR = mu2*Q;
     setV(Wplus,Q,alphat,Nf,f_quiet);
     plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NNLO2Matrix.fillRow(j,3,y,temp_ratio,temp_err);
   }
//
// print NNLO:
  std::cout << std::endl << "set color red" << std::endl << std::endl;
   std::cout << "title 0.05 25.5 data 'NNLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NNLO2Matrix[j][0]
            << "       " << NNLO2Matrix[j][1] << " "
            << "   (   " << NNLO2Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO2Matrix[j][0]
            << "       " << NNLO2Matrix[j][1] << " "
            << "   (   " << NNLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j > 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO1Matrix[j][0]
            << "       " << NNLO1Matrix[j][1] << " "
            << "   (   " << NNLO1Matrix[j][2] << std::endl;
   }
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << -NNLO1Matrix[j][0]
            << "       " << NNLO1Matrix[j][1] << " "
            << "   (   " << NNLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
}

// Same as "sym_scan_all", but handles case of an asymmetric 
// rapidity distribution, e.g. W production at Tevatron.
// For n_points = 38, agrees with Frank's choice of points.

void asym_scan_all(pdf_type pdftype, int n_points, double mu1, double mu2){
 double ymax = log(E_CM/Q) * n_points/(n_points+2.);  
 std::cout << std::endl
      << "-----------------------------------------------------------" << std::endl;
 std::cout << "LO, NLO, NNLO rapidity scan - asymmetric case (|y| < ymax)" << std::endl;
 std::cout << "-----------------------------------------------------------" << std::endl;
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
 std::cout << " muF_lower = muR_lower = " << mu1 << " * Q " << std::endl;
 std::cout << " muF_upper = muR_upper = " << mu2 << " * Q " << std::endl << std::endl;
 std::cout << " pdf type = " << pdftype << " (2 = mrst, 3 = alekhin) " 
      << std::endl << std::endl;
  double step_size = 2.*ymax/n_points; 
  f_NNLO_only = 1;   compute_all();
//
// Do both LO curves first:
  order_flag = 0;  
  if (pdftype==mrst) { alpha_s_Z = 0.130;  
	LHApdfInit(mrst,1);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
  alpha_s_Z = Alekhin_alpha_s(m_Z);
  }
  DMatrix LO1Matrix(0,n_points,0,2);
  DMatrix LO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = - ymax + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     LO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     LO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// Next do both NLO curves:
  order_flag = 1;  
  if (pdftype==mrst) { alpha_s_Z = 0.119; 
	LHApdfInit(mrst,2);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
 alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NLO1Matrix(0,n_points,0,2);
  DMatrix NLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = -ymax + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     NLO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     NLO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// Finally do both NNLO curves:
  order_flag = 2;  
  if (pdftype==mrst) { alpha_s_Z = 0.1155;  
	LHApdfInit(mrst,6);
}
  else if (pdftype==alekhin) {
	LHApdfInit(alekhin,order_flag);
      alpha_s_Z = Alekhin_alpha_s(m_Z); 
  }
  DMatrix NNLO1Matrix(0,n_points,0,2);
  DMatrix NNLO2Matrix(0,n_points,0,2);
  for (int j = 0; j <= n_points; j++) {
     y = -ymax + step_size * j;
     muF = mu1*Q;   muR = mu1*Q;
     DVector temp_rap_y = rap_y();
     NNLO1Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
//  
     muF = mu2*Q;   muR = mu2*Q;
     temp_rap_y = rap_y();
     NNLO2Matrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
//
// print results: 
// print LO:
  std::cout << std::endl << "set color blue" << std::endl << std::endl;
   std::cout << "title 0.05 10.5 data 'LO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO2Matrix[j][0]
            << "       " << LO2Matrix[j][1] << " "
            << "   (   " << LO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j >=0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LO1Matrix[j][0]
            << "       " << LO1Matrix[j][1] << " "
            << "   (   " << LO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// print NLO:
  std::cout << std::endl << "set color green" << std::endl << std::endl;
   std::cout << "title 0.05 20.5 data 'NLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO2Matrix[j][0]
            << "       " << NLO2Matrix[j][1] << " "
            << "   (   " << NLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j >= 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLO1Matrix[j][0]
            << "       " << NLO1Matrix[j][1] << " "
            << "   (   " << NLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
// print NNLO:
  std::cout << std::endl << "set color red" << std::endl << std::endl;
   std::cout << "title 0.05 25.5 data 'NNLO' SIZE 2.5" << std::endl << std::endl;
   std::cout << "( mu/Q = " << mu2 << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO2Matrix[j][0]
            << "       " << NNLO2Matrix[j][1] << " "
            << "   (   " << NNLO2Matrix[j][2] << std::endl;
   }
   std::cout << "( mu/Q = " << mu1 << std::endl;
   for (int j = n_points; j >= 0; j--) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLO1Matrix[j][0]
            << "       " << NNLO1Matrix[j][1] << " "
            << "   (   " << NNLO1Matrix[j][2] << std::endl;
   }
   std::cout << "join fill" << std::endl;
}

// =========== SCAN IN REN/FACT. SCALES mu ============================

// At fixed y, scan muF = muR from mu_r_lower*Q to mu_r_upper*Q,
// (or simple variations on this, like muF only, muR only, etc.)
// with mu values uniformly spaced in log(mu), and print results.
// X is the variable on the horizontal axis for plotting.
// Do LO, NLO and NNLO MRST curves at once.

void scan_mu(double yy, double mu_r_lower, double mu_r_upper,
                int n_points){
 std::cout << std::endl
      << "-------------------------------------------------------------" 
      << std::endl; 
 std::cout << "Performing mu scan from " << mu_r_lower << " to " 
      << mu_r_upper << " in " << n_points 
      << " steps;  at rapidity y = " << yy << std::endl; 
 std::cout << "-------------------------------------------------------------"
      << std::endl; 
  y = yy;
  double n_points_real = n_points;
  double step_size = log(mu_r_upper/mu_r_lower)/n_points_real;
  double X;
// matrices to store results:
  DMatrix LOMatrix(0,n_points,0,2); 
  DMatrix NLOMatrix(0,n_points,0,2);  
  DMatrix NNLOMatrix(0,n_points,0,2);  
  f_NNLO_only = 1;   compute_all();
  for (int j = 0; j <= n_points; j++) {
     muF = Q * mu_r_lower * exp(step_size * j);
     muR = Q*Q/muF;
     X = muF/Q;
// LO MRST:
     order_flag = 0;  alpha_s_Z = 0.130;  
	LHApdfInit(mrst,1);
     DVector temp_rap_y = rap_y();
     LOMatrix.fillRow(j,3,X,temp_rap_y[0],temp_rap_y[1]);
// NLO MRST:
     order_flag = 1;  alpha_s_Z = 0.119;  
	LHApdfInit(mrst,2);
     temp_rap_y = rap_y();
     NLOMatrix.fillRow(j,3,X,temp_rap_y[0],temp_rap_y[1]);
// NNLO MRST:
     order_flag = 2;  alpha_s_Z = 0.1155;  
	LHApdfInit(mrst,6);
     temp_rap_y = rap_y();
     NNLOMatrix.fillRow(j,3,X,temp_rap_y[0],temp_rap_y[1]);
  }
// print results: 
// print LO:
  std::cout << std::endl << "set color blue" << std::endl << std::endl;
   std::cout << "title 0.05 10.5 data 'LO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LOMatrix[j][0]
            << "       " << LOMatrix[j][1] << " "
            << "   (   " << LOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
// print NLO:
  std::cout << std::endl << "set color green" << std::endl << std::endl;
   std::cout << "title 0.05 20.5 data 'NLO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLOMatrix[j][0]
            << "       " << NLOMatrix[j][1] << " "
            << "   (   " << NLOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
// print NNLO:
  std::cout << std::endl << "set color red" << std::endl << std::endl;
   std::cout << "title 0.05 25.5 data 'NNLO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLOMatrix[j][0]
            << "       " << NNLOMatrix[j][1] << " "
            << "   (   " << NNLOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
   std::cout << std::endl ;
}

// Same as "scan mu", but does the W+/W- ratio.
// X is the variable on the horizontal axis for plotting.

void Wpm_ratio_scan_mu(double yy, double mu_r_lower, double mu_r_upper,
                int n_points){
 double temp_ratio, temp_err;
 Q = m_W;  // assumes we are on-shell...
 f_quiet = 1;  // turn off most output
 std::cout << std::endl
      << "-------------------------------------------------------------" 
      << std::endl; 
 std::cout << "Performing W+/W- ratio mu scan from " << mu_r_lower << " to " 
      << mu_r_upper << " in " << n_points 
      << " steps;  at rapidity y = " << yy << std::endl; 
 std::cout << "-------------------------------------------------------------"
      << std::endl; 
  y = yy;
  double n_points_real = n_points;
  double step_size = log(mu_r_upper/mu_r_lower)/n_points_real;
  double X;
// matrices to store results:
  DMatrix LOMatrix(0,n_points,0,2); 
  DMatrix NLOMatrix(0,n_points,0,2);  
  DMatrix NNLOMatrix(0,n_points,0,2);  
  f_NNLO_only = 1;   compute_all();
  for (int j = 0; j <= n_points; j++) {
     std::cout << "Now working on point " << j << " of " << n_points << std::endl; 
     muR = Q * mu_r_lower * exp(step_size * j);
     muF = Q;
     X = muR/Q;
// LO MRST:
     order_flag = 0;  alpha_s_Z = 0.130; 
	LHApdfInit(mrst,1);
     setV(Wplus,Q,alphat,Nf,f_quiet);
     DVector plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     DVector minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     LOMatrix.fillRow(j,3,X,temp_ratio,temp_err);
// NLO MRST:
     order_flag = 1;  alpha_s_Z = 0.119;  
	LHApdfInit(mrst,2);
     setV(Wplus,Q,alphat,Nf,f_quiet);
     plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NLOMatrix.fillRow(j,3,X,temp_ratio,temp_err);
// NNLO MRST:
     order_flag = 2;  alpha_s_Z = 0.1155;  
	LHApdfInit(mrst,6);
     setV(Wplus,Q,alphat,Nf,f_quiet);
     plus_rap_y = rap_y();
     setV(Wminus,Q,alphat,Nf,f_quiet);
     minus_rap_y = rap_y();
     temp_ratio = plus_rap_y[0]/minus_rap_y[0];
     temp_err = 1/minus_rap_y[0] * sqrt( 
                  plus_rap_y[1]*plus_rap_y[1] 
               + plus_rap_y[0]*plus_rap_y[0]/minus_rap_y[0]/minus_rap_y[0]
		  * minus_rap_y[1]*minus_rap_y[1] );
     NNLOMatrix.fillRow(j,3,X,temp_ratio,temp_err);
  }
// print results: 
// print LO:
  std::cout << std::endl << "set color blue" << std::endl << std::endl;
   std::cout << "title 0.5 1.1 data 'LO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << LOMatrix[j][0]
            << "       " << LOMatrix[j][1] << " "
            << "   (   " << LOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
// print NLO:
  std::cout << std::endl << "set color green" << std::endl << std::endl;
   std::cout << "title 0.5 1.2 data 'NLO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NLOMatrix[j][0]
            << "       " << NLOMatrix[j][1] << " "
            << "   (   " << NLOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
// print NNLO:
  std::cout << std::endl << "set color red" << std::endl << std::endl;
   std::cout << "title 0.5 1.3 data 'NNLO' SIZE 2.5" << std::endl << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << NNLOMatrix[j][0]
            << "       " << NNLOMatrix[j][1] << " "
            << "   (   " << NNLOMatrix[j][2] << std::endl;
   }
   std::cout << "join" << std::endl;
   std::cout << std::endl ;
}

// Like "scan_mu" above, but for the total cross section, computed
// "the hard way":   FIX !!!!!!!

DMatrix scan_mu_tot(double mu_r_lower, double mu_r_upper,
                int n_points){
  double n_points_real = n_points;
  double step_size = log(mu_r_upper/mu_r_lower)/n_points_real;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     muF = Q * mu_r_lower * exp(step_size * j);
     muR = muF;
     // full result:
     DVector temp_int_y = int_y(0.,1.);
     resultMatrix.fillRow(j,3,muF/Q,temp_int_y[0],temp_int_y[1]);
   }

   // print absolute results again, in form convenient for copying: 
   std::cout << "   (muF=muR)/Q    |    dsigma/dM  |   (  error   " << std::endl;
   std::cout << "------------------------------------------------------- " << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
            << "       " << resultMatrix[j][1] << " "
            << "   (   " << resultMatrix[j][2]   // include errors
            << std::endl;
   }
   std::cout << std::endl ;
//
 return resultMatrix;
}

// =========== INTEGRATE OVER INVARIANT MASS M ==================

// "int_M_sym_scan_y" 
// integrates over invariant masses M, from M_lower to M_upper,
// just using Simpson's rule with m_points in the M direction,
// and also scans in rapidity (symmetric version).  
// Very time-consuming, so probably cannot
// be run at NNLO for more than a couple of rapidity points at a time.
// If last argument, "containsZ" is 1, we remap the interval to take
// into account the Z peak. 
// At the moment, though, use "containsZ = 0"; the remapping seems
// to perform worse than the original at the same number of points???

void int_M_sym_scan_y(exchange exch, double M_lower, double M_upper, 
   int m_points, int n_points, int containsZ){
 double u_l, u_u, u, temp_prefactor, tempres, temperr;
 double Gamma_Z = 2.4952;
 double ymax = log(E_CM/Q) * n_points/(n_points+1.);  
 std::cout << std::endl
      << "--------------------------------------------------------" << std::endl; 
 std::cout << "Performing M-integration rapidity scan - "
      << "symmetric case (0 < y < ymax)" << std::endl; 
 std::cout << "--------------------------------------------------------" << std::endl; 
 std::cout << " M_lower =  " << M_lower << ";  M_upper =  " << M_upper 
      << ";    number of mass steps (must be even!) =  " 
      << m_points << std::endl << std::endl;
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
 setV(exch,M_lower,alphat,Nf,0);
 f_quiet = 1;
 double step_size = ymax/n_points;
 double M_step_size = (M_upper-M_lower)/m_points;
 if (containsZ==1) { 
   Gamma_Z = 2.4952;
   u_l = atan((M_lower*M_lower-m_Z*m_Z)/Gamma_Z/m_Z);
   u_u = atan((M_upper*M_upper-m_Z*m_Z)/Gamma_Z/m_Z);
   M_step_size = (u_u-u_l)/m_points;
 }
//
 DMatrix resultMatrix(0,n_points,0,2);  // matrix to store results
 for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
  // y = -ymax + step_size * j;
     // full result/error:
     double temp_sum_y = 0.;    double temp_err_y = 0.;
      for (int k = 0; k <= m_points; k++) {
        Q = M_lower + M_step_size * k;
        if (containsZ==1) {
	  u = u_l + M_step_size * k;
          Q = sqrt(m_Z*m_Z + Gamma_Z*m_Z*tan(u));
	}
  //      muR = muRrel*Q;    muF = muFrel*Q;
	setV(exch,Q,alphat,Nf,1);
        DVector temp_rap_y = rap_y();
        temp_prefactor = 1.;
        if (containsZ==1) { 
          temp_prefactor = Gamma_Z*m_Z*(1.+tan(u)*tan(u))/2./Q;
	}     
        tempres = temp_prefactor * temp_rap_y[0];
        temperr = temp_prefactor * temp_rap_y[1];
        if ((k==0) || (k==m_points)) {   // endpoints in Simpson's rule:
	  temp_sum_y = temp_sum_y + M_step_size/3. * tempres;
	  temp_err_y = temp_err_y + M_step_size*M_step_size/9. 
                                       * temperr*temperr;
	}
        else if (k/2.==INT(k/2.)) {     // even points:
	  temp_sum_y = temp_sum_y + M_step_size * 2./3. * tempres;
	  temp_err_y = temp_err_y + M_step_size*M_step_size * 4./9. 
                                       * temperr*temperr;
        }
        else {                          // odd points:
	  temp_sum_y = temp_sum_y + M_step_size * 4./3. * tempres;
	  temp_err_y = temp_err_y + M_step_size*M_step_size * 16./9. 
                                       * temperr*temperr;
        }
        temp_err_y = sqrt(temp_err_y);
        resultMatrix.fillRow(j,3,y,temp_sum_y,temp_err_y);
      }
 }
// print results: 
   std::cout << std::endl
        << "   y               |    d^2sigma/dM/dy   |   (  error   " << std::endl;
   std::cout << "------------------------------------------------------- " << std::endl;
   for (int j = 0; j <= n_points; j++) {
     std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
            << "       " << resultMatrix[j][1] << " "
            << "   (   " << resultMatrix[j][2]   // include errors
            << std::endl;
   }
   std::cout << std::endl ;
}

// Routine computes \int_Ml^Mu dM * d sigma/dM/dy, plus errors:

DVector rap_y_M(){
  double total, total_error;
//  Set relative ren/fact. factors (they are held fixed as we
// integrate over M (Q):
  DVector result(0,1);   
  DVector i_Born(0,2);   DVector i_NLO(0,2);  DVector i_NNLO(0,2);
  if (f_quiet==0) {
    std::cout << std::setw(9) << std::setprecision(9) << " working on    y =  " << y << std::endl; 
  } 
// LO case
  if (order_flag == 0) { 
    if (std::fabs(Born_integrand_M(0.2,20.)) < 1.0e-16){
      std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      result.fill(2,0.,0.); 
    }
    else {
      std::cout << "   Now integrating Born terms " << std::endl;
      i_Born = sqint(6,2,10000,10,30,ranseed); 
      result.fill(2,i_Born[0],i_Born[1]) ;
    }
  }
// NLO case
  else if (order_flag == 1) {
    if (std::fabs(Born_integrand_M(0.2,20.)) < 1.0e-16){
      std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      i_Born.fill(3,0.,0.,0.);
    }
    else {
      std::cout << "   Now integrating Born terms " << std::endl;
      i_Born = sqint(6,2,5000,10,20,ranseed); 
    }
    if (std::fabs(int_NLO_M(.03,.43,.93,20.)) < 1.0e-16) {
       if (f_quiet==0) {
       std::cout << "Setting NLO non-Born terms to 0, since integrand ~ 0 " << std::endl;
       }
       i_NLO.fill(3,0.,0.,0.);
    }
    else { i_NLO = sqint(7,3,5000,10,30,ranseed); }
    total = i_Born[0] + i_NLO[0];
    total_error = sqrt(i_Born[1]*i_Born[1] + i_NLO[1]*i_NLO[1]);
    result.fill(2,total,total_error);
  }
// NNLO case
  else {   //  order_flag = 2:
    if (f_NNLO_only == 1) {
      if (std::fabs(Born_integrand_M(0.2,20.)) < 1.0e-16){
        std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
      i_Born.fill(3,0.,0.,0.);
      }
      else {
        std::cout << "   Now integrating Born terms " << std::endl;
        i_Born = sqint(6,2,5000,10,20,ranseed); 
      }
      if (f_quiet==0) {
        std::cout << " computing NLO integral for y = " << y << std::endl; 
      }
      if (std::fabs(int_NLO_M(.03,.43,.93,20.)) < 1.0e-16) {
        if (f_quiet==0) {
        std::cout <<"Setting NLO non-Born terms to 0, since integrand ~ 0 "<< std::endl;
        }
        i_NLO.fill(3,0.,0.,0.);
      }
      else { i_NLO = sqint(7,3,5000,10,30,ranseed); }
    }
    if (f_quiet==0) {
    std::cout << " computing NNLO integral for y = " << y << std::endl; 
    }
    i_NNLO = sqint(8,3,20000,5,15,ranseed);
    if (f_NNLO_only == 0) {
      total = i_Born[0] + i_NNLO[0];
      total_error = sqrt(i_Born[1]*i_Born[1] + i_NNLO[1]*i_NNLO[1]);
    }
    else {
      total = i_Born[0] + i_NLO[0] + i_NNLO[0];
      total_error = sqrt( i_Born[1]*i_Born[1] + i_NLO[1]*i_NLO[1]
                                             + i_NNLO[1]*i_NNLO[1] );
    }
    result.fill(2,total,total_error);
  }
  if (f_quiet==0) {
  std::cout << std::setw(8) << std::setprecision(8) <<  " y = " << y 
       << ";  int_ml^Mu dM d^2sigma/dM/dy =  " << result[0]
       << "   pm  " << result[1] << std::endl;
  }
  return result;
}

// Routine uses "rap_y_M" to compute un-normalized 
// integral_Mu^Ml dM * d sigma/dM/dy, plus errors,
// in the SYMMETRIC case, just scans from y = 0 to y = ymax,
// with a number of points equal to n_points:
// ymax is computed from E_CM and Q; then we rescale to omit the last point.
// For n_points = 19, this agrees with Frank's choice of points.
// forwrev = 1  -> print output from -ymax to +ymax
// forwrev = -1  -> print output from +ymax to -ymax

DMatrix sym_scan_rap_y_M(exchange exch, double M_l, double M_u, 
            int n_points, int forwrev){
  Ml = M_l;  Mu = M_u;   exchM = exch;
  double ymax = log(E_CM/Q) * n_points/(n_points+1.);  
 std::cout << std::endl
      << "--------------------------------------------------------" << std::endl; 
 std::cout << "Performing rapidity scan - symmetric case (0 < y < ymax)" << std::endl; 
 std::cout << "Includes integration over M from " 
      << Ml << "  to " << Mu << " GeV " << std::endl; 
 std::cout << "--------------------------------------------------------" << std::endl; 
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
  double step_size = ymax/n_points;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     y = 0. + step_size * j;
  // y = -ymax + step_size * j;
     // full result:
     DVector temp_rap_y = rap_y_M();
     resultMatrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// print results for topdraw:
    std::cout << std::endl
         << "   y       |  int_Ml^Mu dM d^2sigma/dM/dy  |  ( error   " << std::endl;
    std::cout << "------------------------------------------------------- " << std::endl;
    for (int j = n_points; j > 0; j--) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << -forwrev*resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl;
    }
    for (int j = 0; j <= n_points; j++) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << forwrev*resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl;
    }
   std::cout << std::endl ;
//
 return resultMatrix;
}

// Same as "sym_scan_rap_y_M", except for the ASYMMETRIC case, 
// e.g. W (W') production at the Tevatron. 
// It scans from y = -ymax to y = ymax, with a number of points equal to n_points:
// ymax is computed from E_CM and Q; then we rescale to omit the last
// point on each end.  Try e.g. n_points = 38.
// forwrev = 1  -> print output from -ymax to +ymax
// forwrev = -1  -> print output from +ymax to -ymax

DMatrix asym_scan_rap_y_M(exchange exch, double M_l, double M_u, 
            int n_points, int forwrev){
  Ml = M_l;  Mu = M_u;   exchM = exch;
  double ymax = log(E_CM/Q) * n_points/(n_points+2.);  
 std::cout << std::endl
      << "--------------------------------------------------------" << std::endl; 
 std::cout << "Performing rapidity scan - asymmetric case (|y| < ymax)" << std::endl; 
 std::cout << "Includes integration over M from " 
      << Ml << "  to " << Mu << " GeV " << std::endl; 
 std::cout << "--------------------------------------------------------" << std::endl; 
 std::cout << " ymax =  " << ymax
      << ";    number of rapidity steps =  " << n_points << std::endl << std::endl; 
  double step_size = 2.*ymax/n_points;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     y = - ymax + step_size * j;
     // full result:
     DVector temp_rap_y = rap_y_M();
     resultMatrix.fillRow(j,3,y,temp_rap_y[0],temp_rap_y[1]);
   }
// print results for topdraw:
    std::cout << std::endl
         << "   y       |  int_Ml^Mu dM d^2sigma/dM/dy  |  ( error   " << std::endl;
    std::cout << "------------------------------------------------------- " << std::endl;
    if (forwrev==1) {
    for (int j = 0; j <= n_points; j++) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl; }
    }
    if (forwrev==-1) {
    for (int j = n_points; j >= 0; j--) {
      std::cout << std::setw(8) << std::setprecision(8) <<  " " << resultMatrix[j][0]
           << "       " << resultMatrix[j][1] << " "
           << "   (   " << resultMatrix[j][2] << std::endl; }
    }
   std::cout << std::endl ;
//
 return resultMatrix;
}

#if USE_OLD

//====================== TESTS ======================================

// Test the luminosity functions:

void lumi_test(){
 double x_1, x_2, mu_temp;
 std::cout << "type x1"  << std::endl;
 cin >> x_1;
 std::cout << "type x2"  << std::endl;
 cin >> x_2;
 std::cout << "type mu"  << std::endl;
 cin >> mu_temp;
 std::cout << " x1 = " << x_1 << "  ;   x2 = " << x_2 
      << "  ;   mu = " << mu_temp << std::endl;
 std::cout << std::setw(8) << std::setprecision(8) << " qqbarNS(x1,x2,mu)  =   " 
      << qqbar_lumi(x_1,x_2,mu_temp,DY,coll) << std::endl;
 std::cout << " qqbarB2(x1,x2,mu)  =   " 
      << qqbar_lumi(x_1,x_2,mu_temp,gluon,coll) << std::endl;
 std::cout << " qqbarBC(x1,x2,mu)  = [exactly 2 * NS for (gamma,Z) case] " 
      << qqbar_BC_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " qqbarABv(x1,x2,mu) not required; hard cross section vanishes  " 
      << std::endl;
 std::cout << " qqbarABax(x1,x2,mu)  =   " 
      << qqbar_ax_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " qg_lumi(x1,x2,mu) = " << qg_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " gq_lumi(x1,x2,mu) = " << gq_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiC2(x1,x2,mu)  =   " 
      << qq_11_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiD2(x1,x2,mu)  =   " 
      << qq_22_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiCD_V(x1,x2,mu)  =   " 
      << qq_12_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiCD_ax(x1,x2,mu)  =   " 
      << qq_12_ax_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiCE1(x1,x2,mu)  =   " 
      << qq_CE1_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiCE2(x1,x2,mu)  =   " 
      << qq_CE2_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " lumiCF(x1,x2,mu)  = [same as CE for (gamma,Z) case] " 
      << qq_CF_lumi(x_1,x_2,mu_temp,coll) << std::endl;
 std::cout << " gg_lumi(x1,x2,mu) = " << gg_lumi(x_1,x_2,mu_temp,coll) << std::endl;
}
#endif
