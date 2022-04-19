/* =============================================================
   The functions required for the Born kinematics and
   NLO contributions to the DY rapidity distributions.
   ============================================================= */

#include <iomanip>
#include "dilog.h"

// The NLO Born kinematics terms:
double Born_NLO(double muFQ){
  double C_F = 4./3.;  double lnQF = -2.*log(muFQ);
  return C_F * ( 1.5 * lnQF + 2. * ZETA2 - 4. );
}

// The NNLO Born kinematics terms:
double Born_NNLO(double Nf, double muFQ, double muRQ){
  double Z2 = ZETA2;  double Z3 = ZETA3;
  double C_F = 4./3.;  double C_A = 3.;  
  double lnQF = -2.*log(muFQ);  double lnRF = 2.*log(muRQ/muFQ);
  return 1./16. * ( 
    C_A * C_F * (
    - 11. * lnQF*lnQF + (193./3. - 24 * Z3) * lnQF
    - 12./5. * Z2*Z2 + 28. * Z3 + 592./9. * Z2 - 1535./12. )
  + C_F * C_F * ( ( 18. - 32. * Z2 ) * lnQF*lnQF
    + (24. * Z2 + 176. * Z3 - 93.) * lnQF
		+ 8./5. * Z2*Z2 - 70. * Z2 - 60. * Z3 + 511./4. )
  + Nf * C_F * ( 2. * lnQF*lnQF - 34./3. * lnQF
		 + 8. * Z3 - 112./9. * Z2 + 127./6. )
// terms for muR <> muF:
  + C_F * (11./3. * C_A - 2./3. * Nf) * lnRF 
        * ( 6. * lnQF + 8. * Z2 - 16. ) );
}

// The NLO boost terms, and "soft" (z=1) subtraction (for q-\bar{q]):

double NLO_qbarq_boost(double z, double muFQ){
  double C_F = 4./3.;
  return  C_F/2. * ( (1.+z*z)/(1.-z) * ( 2.*log((1.-z)/muFQ) - log(z) ) 
		   + 1. - z ) ;
}
 
double NLO_qbarq_boost_soft(double z, double muFQ){
  double C_F = 4./3.;
  return C_F/2. * 4. * log((1.-z)/muFQ) /(1.-z) ;
}

double NLO_qg_boost(double z, double muFQ){
  double T_R = 0.5;
  return T_R * ( 0.5 * (z*z+(1.-z)*(1.-z)) * ( 2.*log((1.-z)/muFQ) - log(z) )
	       + z - z*z ) ;
}

// The NLO real terms, and y=0, y=1, subtractions:

double NLO_qbarq_real(double y, double z){
  double C_F = 4./3.;
  return C_F/2. * ( (1.+z*z)/(1.-z) * ( 1./y + 1./(1.-y) ) 
		  + 2. * (-1.+z) );
}

double NLO_qbarq_real_soft(double y, double z){
  double C_F = 4./3.;
  return C_F/2. * (1.+z*z)/(1.-z) /y ;
}

double NLO_qg_real(double y, double z){
  double T_R = 0.5;
  return T_R/2. * ( (z*z+(1.-z)*(1.-z)) /y 
		    + 2.*z*(1.-z) + (1.-z)*(1.-z)*y ) ;
}

double NLO_qg_real_soft(double y, double z){ 
  double T_R = 0.5;
  return T_R/2. * (z*z+(1.-z)*(1.-z)) /y ;
}
