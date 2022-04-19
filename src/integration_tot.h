
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
  double tau = Q*Q/E_CM/E_CM;
  double asR = alpha_s(muR)/PI;
  double x1 = sqrt(tau) * exp(y);
  double x2 = sqrt(tau) * exp(-y);
  if ((x1 >= 1.) || (x2 >= 1.)){ return 0.; }
  double LO_term = 1.;
  double NLO_term = SV_Born_NLO(muF/Q);
  double NNLO_term = SV_Born_NNLO(Nf,muF/Q,muR/Q);
  double Born_term;
  if (order_flag == 0) {
    Born_term = LO_term; }
  if (order_flag == 1) {
    Born_term = LO_term + asR * NLO_term; }
  if (order_flag == 2) {
    if (f_NNLO_only == 0) { Born_term = asR*asR * NNLO_term; }
    else { Born_term = LO_term + asR * NLO_term + asR*asR * NNLO_term; }
  }
	pdfArray X1,X2;
	X1.x=x1; X2.x=x2;
	LHAComputePdf(x1,muF,X1);	
	LHAComputePdf(x2,muF,X2);	
	return Born_term * qqbar_lumi(X1,X2,DY,coll);

}

//================================================================
// The boost-type terms in the integrand.  
// Here NLO and NNLO are combined.
// Form after change of variables to x1h, zh (much better convergence):

double boost_integrand(double x1h, double zh){
  double asR, tau, c, x1, d, z, x2, Jx1hzh, lumiz, lumi0,  
         qbarq_BC_lumiz, qbarq_NFf_lumiz, qbarq_ax_lumiz,
         gq_lumiz, qg_lumiz, gg_lumiz,
         qq_11_lumiz, qq_22_lumiz, qq_12_lumiz, qq_12_ax_lumiz,
         qq_CE1_lumiz, qq_CF_lumiz, main_ans_NLO, main_ans_NNLO;
  asR = alpha_s(muR)/PI; 
  tau = Q*Q/E_CM/E_CM;     c = -log(tau);    x1 = exp(-c*x1h) ;
  d = c + log(x1);   z = exp(-d*zh);  x2 = tau/x1/z;
  Jx1hzh = c * d * z ;
 
	pdfArray X1,X2;
	X1.x=x1; X2.x=x2;
	LHAComputePdf(x1,muF,X1);	
	LHAComputePdf(x2,muF,X2);	
 lumiz = qqbar_lumi(X1,X2,DY,coll); 
  qbarq_BC_lumiz = qqbar_BC_lumi(X1,X2,coll); 
  qbarq_NFf_lumiz = qqbar_lumi(X1,X2,gluon,coll);
  qbarq_ax_lumiz = qqbar_ax_lumi(X1,X2,coll);
  qg_lumiz = qg_lumi(X1,X2,coll);
  gq_lumiz = gq_lumi(X1,X2,coll);
  gg_lumiz = gg_lumi(X1,X2,coll);
  qq_11_lumiz = qq_11_lumi(X1,X2,coll);
  qq_22_lumiz = qq_22_lumi(X1,X2,coll);
  qq_12_lumiz = qq_12_lumi(X1,X2,coll);
  qq_12_ax_lumiz = qq_12_ax_lumi(X1,X2,coll);
  qq_CE1_lumiz = qq_CE1_lumi(X1,X2,coll);
  qq_CF_lumiz = qq_CF_lumi(X1,X2,coll);
  // for z = 1 subtraction (q qbar only):
  x2 = tau/x1;
	X2.x=x2; 
	LHAComputePdf(x2,muF,X2);	
  lumi0 = qqbar_lumi(X1,X2,DY,coll);

  main_ans_NLO = Jx1hzh * (
  // qqbar
     ( SV_boost_NLO(z,muF/Q) + qbarq_hard_NLO(z,muF/Q) ) * lumiz/z
   - ( SV_boost_NLO(z,muF/Q) + SV_boost_NLOe(x2,muF/Q) ) * lumi0
  // qg
   + qg_NLO(z,muF/Q) * (gq_lumiz + qg_lumiz)/z 
                          ) ;
  if (order_flag == 1) { return asR * main_ans_NLO ; }
  else if (order_flag == 2) {
    // SOFT TEST VERSION:
    /*
      main_ans_NNLO =
      SV_boost_NNLO_TEST(z,Nf,muF/Q,muR/Q) * z * lumiz/z
    - ( SV_boost_NNLO_TEST(z,Nf,muF/Q,muR/Q) 
	+ SV_boost_NNLOe_TEST(x2,Nf,muF/Q,muR/Q) ) * lumi0 ;
    */
    // REAL VERSION:
    main_ans_NNLO =
       // q-\bar{q}:  Note that 
      // 1) BC terms now have their own luminosity function
      // 2) CC and CD_V terms are removed to "qq".
    ( SV_boost_NNLO(z,Nf,muF/Q,muR/Q) 
    + qbarq_NS_hard_NNLO(z,Nf,muF/Q,muR/Q) ) * lumiz/z
    - ( SV_boost_NNLO(z,Nf,muF/Q,muR/Q) 
      + SV_boost_NNLOe(x2,Nf,muF/Q,muR/Q) ) * lumi0
       // q-\bar{q} BC terms:
    + qbarq_BC(z) * qbarq_BC_lumiz/z
       // q-\bar{q} NFf terms:
    + qbarq_NFf(z) * qbarq_NFf_lumiz/z
    + AB_ax(z) * qbarq_ax_lumiz/z  // extra axial terms for Z
        // qg:
    + ( qg_C_A(z,muF/Q,muR/Q) + qg_C_F(z,muF/Q,muR/Q) 
      + Nf * qg_NF(z,muF/Q,muR/Q) ) // FIXED missing Nf, 11/21/2003 (LD)
         * (gq_lumiz + qg_lumiz)/z
        // gg:
    + (gg_C_A(z) + gg_C_F(z,muF/Q)) * gg_lumiz/z
        // q_i \neq q_j (or \bar{q}_j)   plus   identical quark terms:
    + qbarq_CC(z,muF/Q) * (qq_11_lumiz + qq_22_lumiz)/z 
    + qbarq_CD_V(z) * qq_12_lumiz/z
    + CD_ax(z) * qq_12_ax_lumiz/z
    + 2. * qq_CE_tot(z,muF/Q) * qq_CE1_lumiz/z
    + qq_CF(z) * qq_CF_lumiz/z
    ;
    main_ans_NNLO = Jx1hzh * main_ans_NNLO;
    if (f_NNLO_only == 0) { return asR*asR * main_ans_NNLO; }
    else { return asR * main_ans_NLO + asR*asR * main_ans_NNLO; }
  }
  else { return 0.; }
}

//============================================================
class surf: public Surface{
public:
  int n;
  surf(int m): n(m){}

  double surface(DVector & x){
    double xi = x[1] ;
    double ytemp = 0.5*log(xi/(1-xi));
    double Jy = 0.5/xi/(1.-xi);
// Integrands for total cross section:
   // introduce change of variables to do integral over y:
   // let xi = (1+tanh(y))/2,  y = ln(x/(1-x))/2.
   // xi = x[1] already runs from 0 to 1.
    if (n==1) {  // integrate Born_integrand over y
      return Jy * Born_integrand(ytemp) ; }
    if (n==2) { 
                   // OLD way of doing y integration:
      //     return Jy * boost_integrand_OLD(ytemp, x[2]) ; 
    // integrate boost integrand over y 
    // (but we actually use x[1] = x1h, x[2] = zh \in [0,1]):
      return boost_integrand(x[1], x[2]) ; 
    }
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
  std::cout << "adapting grid"<<std::endl;
  integral = Sfc.vegas(V,n_adapt,sd,chi2);
  std::cout << std::setw(8) << std::setprecision(8) << integral 
       << "  pm  " << sd << "  chi2 = " << chi2 << std::endl;
  std::cout << "evaluating integral" <<std::endl;
  integral = Sfc.vegasint(V,n_run,sd,chi2);
  std::cout << std::setw(8) << std::setprecision(8) << integral 
       << "  pm  " << sd << "  chi2 = " << chi2 << std::endl;
  std::cout << std::endl;
  temp_out.fill(3,integral,sd,chi2);
  return temp_out;
}

// Routine to integrate over y, from -infinity to +infinity.

DVector int_y(){
  DVector result(0,1);
  DVector Born_int(0,2);
  DVector boost_int(0,2);
  double prefactor = DY_prefactor(Q,alphat);
  if (order_flag == 0) {
    if (std::fabs(Born_integrand(0.2)) < 1.0e-16){
    std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
    Born_int.fill(3,0.,0.,0.);
    }
    else {
    std::cout << "   Now integrating Born terms " << std::endl;
    Born_int = sqint(1,1,1000,10,20,ranseed);
    std::cout << " LO integral =  " << std::setw(8) << std::setprecision(8)
         << prefactor*Born_int[0] << "  pm  " << prefactor*Born_int[1]
         <<  ";    chi^2 =  " << Born_int[2] << std::endl << std::endl;
    result.fill(2,prefactor*Born_int[0],prefactor*Born_int[1]); 
    } 
  }
  if ((order_flag == 1) || (order_flag == 2)) {
    if (std::fabs(Born_integrand(0.2)) < 1.0e-16){
    std::cout << "   Setting Born terms to 0, since f_qqbar = 0  " << std::endl;
    Born_int.fill(3,0.,0.,0.);
    }
    else {
    std::cout << "   Now integrating Born terms " << std::endl;
    Born_int = sqint(1,1,1000,10,20,ranseed);
    std::cout << " Born-type integral =  " << std::setw(8) << std::setprecision(8)
         << prefactor*Born_int[0] << "  pm  " << prefactor*Born_int[1]
         <<  ";    chi^2 =  " << Born_int[2] << std::endl << std::endl;
    }
    std::cout << "   Now integrating boost terms " << std::endl;
    boost_int = sqint(2,2,3000,10,40,ranseed);
    std::cout << " boost-type integral =  " << std::setw(8) << std::setprecision(8)
         << prefactor*boost_int[0] 
         << "  pm  " << prefactor*boost_int[1]
         <<  ";    chi^2 =  " << boost_int[2] << std::endl 
         << std::endl;
    double total_int = prefactor * ( Born_int[0] + boost_int[0] );
    double total_int_error = prefactor * sqrt( Born_int[1]*Born_int[1] 
                            + boost_int[1]*boost_int[1] ) ;
   if (order_flag == 1) {
     std::cout << " NLO integral =  " << std::setw(8) << std::setprecision(8)
        << total_int << "  pm  " << total_int_error  
        << std::endl << std::endl;
   }
   else {
     std::cout << " NNLO integral =  " << std::setw(8) << std::setprecision(8)
        << total_int << "  pm  " << total_int_error  
        << std::endl << std::endl;
   }
  result.fill(2,total_int,total_int_error);  
  }
 return result;
}

// Scan muF = muR from mu_r_lower*Q to mu_r_upper*Q,
// with mu values uniformly spaced in log(mu), and print results:

DMatrix scan_mu(double mu_r_lower, double mu_r_upper, int n_points){
  double n_points_real = n_points;
  double step_size = log(mu_r_upper/mu_r_lower)/n_points_real;
  DMatrix resultMatrix(0,n_points,0,2);  // matrix to store unnorm. results
  for (int j = 0; j <= n_points; j++) {
     muF = Q * mu_r_lower * exp(step_size * j);
     muR = muF;
     // full result:
     DVector temp_int_y = int_y();
     resultMatrix.fillRow(j,3,muF/Q,temp_int_y[0],temp_int_y[1]);
   }

   // print absolute results again, in form convenient for copying: 
   std::cout << "   (muF=muR)/Q    |    d sigma/dM   |   (  error   " << std::endl;
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
