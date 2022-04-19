/* Basic constants and functions for QCD calculations.
  Inherited from "kinematics/kinematics.C(h)" 
  but without maple/kinematics stuff.
  Also, uses different way of flagging order in alpha_s (order_flag), 
  and alpha_s incorporates 3-loop case.  Includes L_i and Ls_i functions */

#include "dilog.h"
#include "QCDbasics.h"


int order_flag;

// 1-loop running of alpha_QED, starting from 0, valid up to about 2 * m_mu
// ("1" includes e only):
double alpha_QED_low(double mu){
  double alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3. );
  return 1./alpha_inv;
}

// 1-loop running of alpha_QED, starting from m_Z, 
// valid between 2*m_b and 2*m_W ("20/3" includes e,mu,tau,d,s,b,u,c):
double alpha_QED_OLD(double mu){
  double alpha_inv = 1/alpha_QED_Z - 2./3./PI * (20./3.) * log(mu/m_Z);
  return 1./alpha_inv;
}

// 1-loop running of alpha_QED, starting from 0, 
// ~ valid between 0 and 2*m_W ("20/3" includes e,mu,tau,d,s,b,u,c).
// artificial quark masses mock up actual hadron vacuum polarization
// (See hep-ph/0306234, Ciccolini, Dittmaier & Kramer.)

double alpha_QED(double mu){
  double m_u_fake = 0.066;
  double m_s_fake = 0.150;
  double alpha_inv;
  if (std::fabs(mu) < 2.*m_e) { alpha_inv = 1/alpha0; }
  else if (std::fabs(mu) < m_u_fake) { // e only in vac. pol.:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3. ); }
  else if (std::fabs(mu) < 2.*m_mu) { // add u and d quarks together here:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. )  ); }
  else if (std::fabs(mu) < 2.*m_s) { // add mu:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3. ); }
  else if (std::fabs(mu) < 2.*m_c) { // add s quark:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) ); }
  else if (std::fabs(mu) < 2.*m_tau) { // add c quark:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) 
                 + 4./3. * ( 2. * log(mu/m_c) - 5./3. ) ); }
  else if (std::fabs(mu) < 2.*m_b) { // add tau:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) 
                 + 4./3. * ( 2. * log(mu/m_c) - 5./3. )
                 + 2. * log(mu/m_tau) - 5./3. ); }
  else if (std::fabs(mu) < 2.*m_W) { // add b:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) 
                 + 4./3. * ( 2. * log(mu/m_c) - 5./3. )
                 + 2. * log(mu/m_tau) - 5./3. 
                 + 1./3. * ( 2. * log(mu/m_b) - 5./3. ) ); }
  else if (std::fabs(mu) < 2.*m_t) { // add W:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) 
                 + 4./3. * ( 2. * log(mu/m_c) - 5./3. )
                 + 2. * log(mu/m_tau) - 5./3. 
                 + 1./3. * ( 2. * log(mu/m_b) - 5./3. )
 // used "Kaplunovsky formula" (not sure about "-5/3") for W case!
                 - 21./2. * log(mu/m_W) ); }
  else { // add t quark:
    alpha_inv = 1/alpha0 - 1./3./PI * ( 2. * log(mu/m_e) - 5./3.
                 + 5./3. * ( 2. * log(mu/m_u_fake) - 5./3. ) 
                 + 2. * log(mu/m_mu) - 5./3.
                 + 1./3. * ( 2. * log(mu/m_s_fake) - 5./3. ) 
                 + 4./3. * ( 2. * log(mu/m_c) - 5./3. )
                 + 2. * log(mu/m_tau) - 5./3. 
                 + 1./3. * ( 2. * log(mu/m_b) - 5./3. ) 
 // used "Kaplunovsky formula" (not sure about "-5/3") for W case!
                 - 21./2. * log(mu/m_W)
                 + 4./3. * ( 2. * log(mu/m_t) - 5./3. ) ); }
  return 1./alpha_inv;
}

double alpha_s_Z;

// General running alpha_s formula:

double alpha_s(double mu){
  if (order_flag==0) { return PI*as1(mu); }
  else if (order_flag==1) { return PI*as2(mu); }
  else if (order_flag==2) { return PI*as3(mu); }
  else return 0.;
}

// 1-loop is only valid for mu > m_b (5 flavors):
double as1(double mu){
  double N_f = 5.;
  double b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
  double Lambda = m_Z * exp(-PI/b0/alpha_s_Z);
  return 1./b0/log(mu/Lambda) ;
}

// 2-loop (from ALEPH Phys. Lett. B284:163 (1992)):
double as2(double mu){
  double w, as, b0, b1, N_f;
  if (mu >= m_b) {
    N_f = 5.;
    b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
    b1 = 0.5 * (17./3.*C_A*C_A-(5./3.*C_A+C_F)*N_f) ;
    as = alpha_s_Z/PI;
    w = 1. - b0 * as * log(m_Z/mu) ;
    return as/w * ( 1. - b1/2./b0 * as * log(w)/w ) ; }
  else {
  // run with 4 flavors below mu = m_b; match with formula
  // alpha_{N_f-1}(M) = alpha_{N_f}(M) [ + O([alpha_{N_f}(M)]^3) ]:
    N_f = 4.;
    b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
    b1 = 0.5 * (17./3.*C_A*C_A-(5./3.*C_A+C_F)*N_f) ;
    as = as2(m_b);
    w = 1. - b0 * as * log(m_b/mu) ;
    return as/w * ( 1. - b1/2./b0 * as * log(w)/w ) ; }
}

/* 3-loop PDG formula for alpha_s(mu)/PI:
   For mu > m_b we run with 5 flavors from alpha_s(m_Z) = alpha_s_Z.
   For mu < m_b we run with 4 flavors and match at m_b using
   alpha_{N_f-1}(M) = alpha_{N_f}(M) + 11/72/Pi^2 * [alpha_{N_f}(M)]^3 */

double as3(double mu){
  double lnx, Lambda_5, Lambda_4, L0, L4, l0, l4,
        b0, b1, b2, alpha_s_5_mb, alpha_s_4_mb, N_f;
  if (mu >= m_b){ 
    N_f = 5.;
    b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
    b1 = 0.5 * (17./3.*C_A*C_A-(5./3.*C_A+C_F)*N_f) ;
    b2 = 0.0625 * ( 2857./27.*C_A*C_A*C_A 
       + (-1415./27.*C_A*C_A-205./9.*C_A*C_F+2.*C_F*C_F)*N_f 
       + (79./27.*C_A+22./9.*C_F)*N_f*N_f ) ;
   // first calculate Lambda_5 from alpha_s_Z.
   // L0 is first approximation to ln(m_Z^2/Lambda_5^2), L4 is 4th:
    L0 = 2.*PI/b0/alpha_s_Z ;
    L4 = iterlnx3(iterlnx3(iterlnx3(iterlnx3(L0,alpha_s_Z,5),
                  alpha_s_Z,5),alpha_s_Z,5),alpha_s_Z,5);
    Lambda_5 = m_Z * exp(-L4/2.) ;
    lnx = 2. * log(mu/Lambda_5);
    return 2./b0/lnx * ( 1. - b1/b0/b0 * log(lnx)/lnx
          + b1*b1/b0/b0/b0/b0/lnx/lnx 
         * ( (log(lnx) - 0.5)*(log(lnx) - 0.5) + b2*b0/b1/b1 - 1.25 ) );
  }

  // run with 4 flavors below mu = m_b; match at m_b with formula
  // alpha_{N_f-1}(M) = alpha_{N_f}(M) + 11/72/Pi^2 * [alpha_{N_f}(M)]^3
  alpha_s_5_mb = PI*as3(m_b);
  alpha_s_4_mb = alpha_s_5_mb 
                   * ( 1. + 11./72./PI/PI * alpha_s_5_mb * alpha_s_5_mb );
  // now calculate Lambda_4 from alpha_s_4_mb.
  // l0 is first approximation to ln(m_b^2/Lambda_4^2), l4 is 4th:
  N_f = 4.;
  b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
  b1 = 0.5 * (17./3.*C_A*C_A-(5./3.*C_A+C_F)*N_f) ;
  b2 = 0.0625 * ( 2857./27.*C_A*C_A*C_A 
     + (-1415./27.*C_A*C_A-205./9.*C_A*C_F+2.*C_F*C_F)*N_f 
              	+ (79./27.*C_A+22./9.*C_F)*N_f*N_f ) ;
  l0 = 2.*PI/b0/alpha_s_4_mb ;
  l4 = iterlnx3(iterlnx3(iterlnx3(iterlnx3(l0,alpha_s_4_mb,4),
               alpha_s_4_mb,4),alpha_s_4_mb,4),alpha_s_4_mb,4);
  Lambda_4 = m_b * exp(-l4/2.) ;
  lnx = 2. * log(mu/Lambda_4);
  return 2./b0/lnx * ( 1. - b1/b0/b0 * log(lnx)/lnx
         + b1*b1/b0/b0/b0/b0/lnx/lnx 
        * ( (log(lnx) - 0.5)*(log(lnx) - 0.5) + b2*b0/b1/b1 - 1.25 ) );
}

// Iterate the 3rd order solution to solve for lnx = ln(m^2/Lambda^2) 
// given alpha_s(m), and for N_f flavors:

double iterlnx3(double lnx, double alpha_s, double N_f){
 double b0 = 0.5 * (11./3.*C_A-2./3.*N_f) ;
 double b1 = 0.5 * (17./3.*C_A*C_A-(5./3.*C_A+C_F)*N_f) ;
 double b2 = 0.0625 * ( 2857./27.*C_A*C_A*C_A 
     + (-1415./27.*C_A*C_A-205./9.*C_A*C_F+2.*C_F*C_F)*N_f 
			+ (79./27.*C_A+22./9.*C_F)*N_f*N_f ) ;
 return 2.*PI/b0/alpha_s * ( 1 - b1/b0/b0 * log(lnx)/lnx
    + b1*b1/b0/b0/b0/b0/lnx/lnx * ( (log(lnx)-.5)*(log(lnx)-.5) 
				    + b2*b0/b1/b1 - 1.25 ) ) ;
}

/* Three-loop PDG formula for alpha_s(mu)/PI:
   For mu > m_b we run with 5 flavors from alpha_s(m_Z) = 0.118.
   For mu < m_b we run with 4 flavors and match at m_b using
   alpha_{N_f-1}(M) = alpha_{N_f}(M) + 11/72/Pi^2 * [alpha_{N_f}(M)]^3 */

double as3_fixed(double mu){
 double b0, b1, b2, N_f, lnx, Lambda;
 if (mu > m_b) { Lambda = .20837;     N_f = 5.; }
 else { Lambda = .28846;    N_f= 4.; }
 b0 = 0.5 * ( 11. - 2./3. * N_f );
 b1 = 0.5 * ( 51. - 19./3. * N_f );
 b2 = 0.0625 * ( 2857. - 5033./9. * N_f + 325./27. * N_f * N_f );
#
 lnx = 2. * log(mu/Lambda);
 return 2./b0/lnx * ( 1. - b1/b0/b0 * log(lnx)/lnx
       + b1*b1/b0/b0/b0/b0/lnx/lnx 
         * ( (log(lnx) - .5) * (log(lnx) - .5) + b2*b0/b1/b1 - 1.25 ) );
}

/* Formulae to run b and c quark MS-bar masses up to higher scales */

// Intermediate procedure used e.g. in HDECAY:

double cq(double x, int Nf){
 double b0, A1, A2;
 b0 = 0.5 * (11.-2./3.*Nf);
 A1 = 499./24. - 7./6. * Nf + 5./54. * Nf*Nf;
 A2 = 6375961./4608. + ( - 775333./3456. - 6655./48. * ZETA3 ) * Nf
     + ( 31573./2304. + 605./24. * ZETA3 ) * Nf*Nf
     + ( 527./2592. - 55./36. * ZETA3 ) * Nf*Nf*Nf
     + ( - 2009./46656. + 5./162. * ZETA3 ) * Nf*Nf*Nf*Nf 
     + 35./34992. * Nf*Nf*Nf*Nf*Nf;
 return pow(x/b0, 2./b0) * ( 1. + x * A1/b0/b0 + x*x * A2/b0/b0/b0/b0 ) ;
}

double M_q(double m_q, double mu1, double mu2, int Nf){
  return m_q * cq(as3(mu2),Nf)/cq(as3(mu1),Nf) ;
}

// M_t(mu) runs the t-quark mass from m_t = m_t(m_t) up to mu,
// OR DOWN TO mu.

double M_t(double mu){
  if ((mu >= m_b) && (mu <= m_t)) { return M_q(m_t,m_t,mu,5) ; }
  else { return M_q(m_t,m_t,mu,6) ; }
}

// M_b(mu) runs the b-quark mass from m_b = m_b(m_b) up to mu.
// Need mu >= m_b.

double M_b(double mu){
  if ((mu >= m_b) && (mu <= m_t)) { return M_q(m_b,m_b,mu,5) ; }
  else { return M_q(M_b(m_t),m_t,mu,6) ; }
}

// M_c(mu) runs the c-quark mass from m_c = m_c(m_c) up to mu:
// Need mu >= m_c.

double M_c(double mu){
  if ((mu >= m_c) && (mu <= m_b)) { return M_q(m_c,m_c,mu,4) ; }
  else if ((mu > m_b) && (mu <= m_t)) { return M_q(M_c(m_b),m_b,mu,5); }
  else { return M_q(M_c(m_t),m_t,mu,6); }
}

/* From Vegas.C: It is necessary to supply a random number generator.
  ran() gives the next random number and ran(n) initializes the sequence 
  when n is a negative integer.  */


inline int signum(double a){return (a < 0 ? -1: 1);}

//==========================================================
// Basic functions of s_{ij}:

// The L_i and Ls_i functions.

complex Clog(double s1, double s2){  // complex log of a ratio, s1/s2
 complex temp = log(ABS(s1/s2));
 if (s1 > 0) { temp = temp - I * PI; }
 if (s2 > 0) { temp = temp + I * PI; }
 return temp;
}

complex Clog1(double s){   // complex log of a single s
 complex temp = log(ABS(s));
 if (s > 0) { temp = temp - I * PI; }
 return temp;
}

complex L0(double s1, double s2){ 
  double r = s1/s2;
  return Clog(s1,s2)/(1.-r);
}

complex L1(double s1, double s2){ 
  double r = s1/s2;
  return (Clog(s1,s2)+1.-r)/(1.-r)/(1.-r);
}

complex L2(double s1, double s2){ 
  double r = s1/s2;
  return (Clog(s1,s2)-(r-1./r)/2.)/(1.-r)/(1.-r)/(1.-r);
}

// Li_2(1-r), for r = s1/s2:

complex CLi2r(double s1, double s2){
  double r = s1/s2;  
  return li2(1.-r) - log(ABS(1.-r)) * I * imag(Clog(s1,s2));
}

complex Lsm1(double s1, double s2, double s3, double s4){
 return CLi2r(s1,s2) + CLi2r(s3,s4) + Clog(s1,s2)*Clog(s3,s4) - PISQ6;
}

complex Ls0(double s1, double s2, double s3, double s4){
 double r1 = s1/s2;
 double r2 = s3/s4; 
 return Lsm1(s1,s2,s3,s4)/(1.-r1-r2);
}

complex Ls1(double s1, double s2, double s3, double s4){
 double r1 = s1/s2;
 double r2 = s3/s4; 
 return ( Ls0(s1,s2,s3,s4) + L0(s1,s2) + L0(s3,s4) )/(1.-r1-r2);
}

complex Ls2(double s1, double s2, double s3, double s4){
 double r1 = s1/s2;
 double r2 = s3/s4; 
 return ( Ls1(s1,s2,s3,s4) + 0.5 * ( L1(s1,s2) + L1(s3,s4) ) )/(1.-r1-r2);
}

complex Ls3(double s1, double s2, double s3, double s4){
 double r1 = s1/s2;
 double r2 = s3/s4; 
 return ( Ls2(s1,s2,s3,s4) + 1./3. * ( L2(s1,s2) + L2(s3,s4) ) 
          - 1./6. * (1./r1+1./r2) )/(1.-r1-r2);
}
