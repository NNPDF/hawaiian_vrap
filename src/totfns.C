/*==============================================================
 Functions for the DY total cross section; 
 formulae from Hamberg, van Neerven, & Matsuura
==============================================================*/

#include <iomanip>
#include "dilog.h"
#include "complex.h"

// "plus" distributions (without the subtraction):

double D1(double n, double y){ return pow(log(y),n)/y ; }

// Functions appearing in extra term required by the change of 
// variables to x1h and z1h.  To be evaluated at z = x2:

double D1e(double n, double y){ return -1./(n+1.)/y * pow(log(y),n+1.) ; }

// NNLO functions ================================

// the NLO delta(1-z) terms:

double SV_Born_NLO(double muFQ){
 double C_F,Z2,lnQF;
 C_F = 4./3.;    Z2 = ZETA2;    lnQF = -2.*log(muFQ);
 return 1./4. * C_F * ( 6. * lnQF + 8. * Z2 - 16. );
}

// the NLO plus distribution terms:

double SV_boost_NLO(double z, double muFQ){
 double C_F, lnQF;
 C_F = 4./3.;    lnQF = -2.*log(muFQ);
 return 1./4. * C_F * ( 8. * D1(0.,1.-z) * lnQF + 16. * D1(1.,1.-z) );
}

// Extra term required by the change of variables to x1h and z1h.
// To be evaluated at z = x2:

double SV_boost_NLOe(double z, double muFQ){
 double C_F, lnQF;
 C_F = 4./3.;    lnQF = -2.*log(muFQ);
 return 1./4. * C_F * ( 8. * D1e(0.,1.-z) * lnQF + 16. * D1e(1.,1.-z) );
}

// the rest of the NLO q-\bar{q} terms:

double qbarq_hard_NLO(double z, double muFQ){
 double C_F, lnQF;
 C_F = 4./3.;    lnQF = -2.*log(muFQ);
 return 1./4. * C_F * ( - 4. * (1.+z) * lnQF - 8. * (1.+z) * log(1.-z)
                        - 4. * (1.+z*z)/(1.-z) * log(z) );
}

// NLO qg terms:

double qg_NLO(double z, double muFQ){
 double T_F, lnQF;
 T_F = 1./2.;   lnQF = -2.*log(muFQ);
 return 1./4. * T_F * ( 2. * (1. + 2.*z*z - 2.*z) *
                              ( 2. * log(1.-z) - log(z) + lnQF )
                      + 1. - 7. * z*z + 6. * z );
}

// NNLO functions ================================

// the NNLO delta(1-z) terms:

double SV_Born_NNLO(double NF, double muFQ, double muRQ){
 double C_A,C_F,Z2,Z3,lnQF,lnRF;
 C_A = 3.;   C_F = 4./3.;   Z2 = ZETA2;   Z3 = ZETA3;
 lnQF = -2.*log(muFQ);   lnRF = 2.*log(muRQ) + lnQF;
 return 1./16. * ( 
  C_A * C_F * (
    - 11. * lnQF*lnQF + ( 193./3. - 24. * Z3 ) * lnQF 
    - 12./5. * Z2*Z2 + 592./9. * Z2 + 28. * Z3 - 1535./12. )
+ C_F * C_F * ( 
      ( 18. - 32. * Z2 ) * lnQF*lnQF
    + ( 24. * Z2 + 176. * Z3 - 93. ) * lnQF
      + 8./5. * Z2*Z2 - 70. * Z2 - 60. * Z3 + 511./4. )
+ NF * C_F * ( 2. * lnQF*lnQF - 34./3. * lnQF 
             + 8. * Z3 - 112./9. * Z2 + 127./6. )
// terms for muR <> muF:
  + C_F * (11./3. * C_A - 2./3. * NF) * lnRF 
  * ( 6. * lnQF + 8. * ZETA2 - 16. ) );
}

// the NNLO plus distribution terms:

double SV_boost_NNLO(double z, double NF, double muFQ, double muRQ){
 double C_A,C_F,Z2,Z3,lnQF;
 C_A = 3.;   C_F = 4./3.;  Z2 = ZETA2;  Z3 = ZETA3;
 lnQF = -2.*log(muFQ);
 return 1./16. * (
 C_A * C_F * ( 
  - 44./3. * D1(0.,1.-z) * lnQF*lnQF
  + ( ( 536./9. - 16. * Z2 ) * D1(0.,1.-z) - 176./3. * D1(1.,1.-z) )
         * lnQF
  - 176./3. * D1(2.,1.-z) + ( 1072./9. - 32. * Z2 ) * D1(1.,1.-z)
  + ( 56. * Z3 + 176./3. * Z2 - 1616./27. ) * D1(0.,1.-z) )
+ C_F * C_F * ( 
    ( 64. * D1(1.,1.-z) + 48. * D1(0.,1.-z) ) * lnQF*lnQF
  + ( 192. * D1(2.,1.-z) + 96. * D1(1.,1.-z) 
  - (128. + 64. * Z2) * D1(0.,1.-z) )
         * lnQF
  + 128. * D1(3.,1-z) - (128. * Z2 + 256.) * D1(1.,1.-z)
    + 256. * Z3 * D1(0.,1.-z) )
 + NF * C_F * ( 8./3. * D1(0.,1.-z) * lnQF*lnQF 
		+ ( 32./3. * D1(1.,1.-z) - 80./9. * D1(0.,1.-z) ) * lnQF
                + 32./3. * D1(2.,1.-z) - 160./9. * D1(1.,1.-z)
                + ( 224./27. - 32./3. * Z2 ) * D1(0.,1.-z) ) )
// extra terms for muR <> muF:
 + 1/4. * (11/3. * C_A - 2/3. * NF) * 2. * log(muRQ/muFQ) 
        * SV_boost_NLO(z,muFQ);
}

// Extra term required by the change of variables to x1h and z1h.
// To be evaluated at z = x2:

double SV_boost_NNLOe(double z, double NF, double muFQ, double muRQ){
 double C_A,C_F,Z2,Z3,lnQF;
 C_A = 3.;   C_F = 4./3.;  Z2 = ZETA2;  Z3 = ZETA3;
 lnQF = -2.*log(muFQ);
 return 1./16. * (
 C_A * C_F * ( 
  - 44./3. * D1e(0.,1.-z) * lnQF*lnQF
  + ( ( 536./9. - 16. * Z2 ) * D1e(0.,1.-z) - 176./3. * D1e(1.,1.-z) )
         * lnQF
  - 176./3. * D1e(2.,1.-z) + ( 1072./9. - 32. * Z2 ) * D1e(1.,1.-z)
  + ( 56. * Z3 + 176./3. * Z2 - 1616./27. ) * D1e(0.,1.-z) )
+ C_F * C_F * ( 
    ( 64. * D1e(1.,1.-z) + 48. * D1e(0.,1.-z) ) * lnQF*lnQF
  + ( 192. * D1e(2.,1.-z) + 96. * D1e(1.,1.-z) 
  - (128. + 64. * Z2) * D1e(0.,1.-z) )
         * lnQF
  + 128. * D1e(3.,1-z) - (128. * Z2 + 256.) * D1e(1.,1.-z)
    + 256. * Z3 * D1e(0.,1.-z) )
 + NF * C_F * ( 8./3. * D1e(0.,1.-z) * lnQF*lnQF 
		+ ( 32./3. * D1e(1.,1.-z) - 80./9. * D1e(0.,1.-z) ) * lnQF
                + 32./3. * D1e(2.,1.-z) - 160./9. * D1e(1.,1.-z)
                + ( 224./27. - 32./3. * Z2 ) * D1e(0.,1.-z) ) )
// extra terms for muR <> muF:
 + 1/4. * (11/3. * C_A - 2/3. * NF) * 2. * log(muRQ/muFQ) 
        * SV_boost_NLOe(z,muFQ);
}

// SOFT TEST VERSIONS of the NNLO plus distribution terms, corresponding
// to integral over y of just the qbarq_real_soft_z1 terms (for muF=muR=Q):

double SV_boost_NNLO_TEST(double z, double NF, double muFQ, double muRQ){
  double Z2,Z3, boost_term, real_soft_term;
 Z2 = ZETA2;  Z3 = ZETA3;
 boost_term = 2. * ( 64./9. * D1(3.,1.-z)
   + ( - 22./3. + 4./9. * NF ) * D1(2.,1.-z)
   + ( - 4./9. * Z2 + 2./3. - 20./27. * NF ) * D1(1.,1.-z)
   + ( 95./9. * Z3 + 11./3. * Z2 - 202./27. 
       + ( - 2./9. * Z2 + 28./81. ) * NF ) * D1(0.,1.-z) );
 real_soft_term = - 64./3. * Z2 * D1(1.,1.-z) 
   + ( 64./3. * Z3 
     + 2. * Z2 * ( 11./3. - 2./9. * NF ) ) * D1(0.,1.-z); 
 return 0.*boost_term + real_soft_term;
}

double SV_boost_NNLOe_TEST(double z, double NF, double muFQ, double muRQ){
 double Z2,Z3, boost_term, real_soft_term;
 Z2 = ZETA2;  Z3 = ZETA3;
 boost_term = 2. * ( 64./9. * D1e(3.,1.-z)
   + ( - 22./3. + 4./9. * NF ) * D1e(2.,1.-z)
   + ( - 4./9. * Z2 + 2./3. - 20./27. * NF ) * D1e(1.,1.-z)
   + ( 95./9. * Z3 + 11./3. * Z2 - 202./27. 
       + ( - 2./9. * Z2 + 28./81. ) * NF ) * D1e(0.,1.-z) );
 real_soft_term = - 64./3. * Z2 * D1e(1.,1.-z) 
   + ( 64./3. * Z3 
     + 2. * Z2 * ( 11./3. - 2./9. * NF ) ) * D1e(0.,1.-z); 
 return 0.*boost_term + real_soft_term;
}

// the q-\bar{q} "NFf" terms:

double qbarq_NFf(double Z){
 double C_F, LOGZ, LOG1PZ, DIMZ;
 C_F = 4./3.;  LOGZ = log(Z);   LOG1PZ = log(1.+Z);   DIMZ = li2(-Z);
 return 1./16. * C_F * (
        LOGZ*LOGZ * ( 4./3. + 4./3.*Z*Z + 8./3.*Z )
      + LOGZ*LOG1PZ * (  - 16./3. - 16./3.*Z*Z
                         - 32./3.*Z )
      + LOGZ * ( 4. + 4.*Z*Z + 16./3.*Z )
      + DIMZ * (  - 16./3. - 16./3.*Z*Z - 32./3.*Z )
      + ZETA2 * (  - 8./3. - 8./3.*Z*Z - 16./3.*Z )
	+ 20./3 - 20./3.*Z*Z );
}

// the q-\bar{q} "NF" terms:

double qbarq_NF(double Z, double muFQ, double muRQ){
 double C_F, ZMIN, LOGZ, LOG1MZ, DI1MZ, FACLM, RENLR;
 C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);   DI1MZ = li2(1.-Z);
 FACLM = -2.*log(muFQ);   RENLR = -2.*log(muRQ);
 return 1./16. * C_F * (  ( 1. + Z ) * ( 16./3.*ZETA2 + 8./3.*DI1MZ -
                  10./3.*LOGZ*LOGZ + 32./3.* LOG1MZ*LOGZ -
                  16./3.*LOG1MZ*LOG1MZ ) +
          ( -64./3.*LOG1MZ*LOGZ - 8./3.*DI1MZ +
            8.*LOGZ*LOGZ + 40./3.*LOGZ ) / ZMIN +
          LOG1MZ * ( - 16./9. + 176./9.*Z ) +
          LOGZ * ( - 4./3 - 44./3.*Z ) +
          188./27. - 412./27.*Z
     +    FACLM*FACLM*( 4./3. + 4./3.*Z ) +
          FACLM*RENLR * (  - 8./3. - 8./3.*Z ) +
          FACLM   *( LOGZ * ( 8./3. + 8./3.*Z
                           - 16./3./ZMIN )
                   - 8./9. + 88./9.*Z                  ) +
          RENLR   *( LOGZ * ( 8./3. + 8./3.*Z
                           - 16./3./ZMIN )
		     + LOG1MZ * (  - 16./3. - 16./3.*Z ) ) );
}

// TEST function:  Integral of just the "hard" part of the NF terms.

double diffNF(double z){ 
  double lnz, ln1mz, lnhf1pz, z2, z3, z4;
  lnz = log(z);  ln1mz = log(1.-z);  lnhf1pz = log((1.+z)/2.);
  z2 = z*z;  z3 = z2*z;  z4 = z2*z2;
  return 1./9./z * (-1.-4.*z-2.*z3+5.*z4+4.*z2)/(-1.+z)/(z+1.) * lnz*lnz
    - 2./9. * (z3+1.+2.*z4*z-2.*z2)/(-1.+z)/(z+1.)/z2 * lnz*ln1mz
    + 2./27. * lnz*(2.*z2-9.*z+11.)/(z+1.) + 4./9. * (-1.+z) * ln1mz
    - 1./9. * (z2+1.)/(-1.+z) * lnhf1pz*lnhf1pz + 4./9. * lnhf1pz
    + 1./27. * (z4+3.*z3-2.)/z/(-1.+z)/(z+1.) * PISQ
    - 2./9. * (z2+1.)/(-1.+z) * li2((1.-z)/2.)
    - 2./9. * (-z3+2.*z4+1.)/(-1.+z)/(z+1.)/z2 * li2(z) + 28./27. * (1.-z);
}

// the q-\bar{q} C_A terms:

double qbarq_C_A(double Z, double muFQ, double muRQ){
 double C_A, C_F, ZMIN, LOGZ, LOG1MZ, DI1MZ, S1MZ, TRI1MZ, FACLM, RENLR;
 C_A = 3.;   C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);   
 S1MZ = S12(1.-Z);   DI1MZ = li2(1.-Z);   TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);   RENLR = -2.*log(muRQ);
 return 1./16. * C_A * C_F * ( ( 1. + Z ) *
               ( 20.*S1MZ - 28.*ZETA3 +
                 16.*LOG1MZ*DI1MZ - 8.*LOGZ*DI1MZ +
                 16.*LOG1MZ*ZETA2 - 8.*LOGZ*ZETA2 +
                 88./3.*LOG1MZ*LOG1MZ ) +
          (-8.*S1MZ - 24.*TRI1MZ
           -16.*LOG1MZ*DI1MZ + 16.*LOGZ*DI1MZ +
            16.*LOGZ*ZETA2 + 8./3.*DI1MZ +
            280./3.*LOG1MZ*LOGZ - 29.*LOGZ*LOGZ
           -208./3.*LOGZ )/ZMIN
           -DI1MZ*( 32./3. + 56./3.*Z )
           -ZETA2*( 76./3. + 100./3.*Z )
           -LOG1MZ*LOGZ*( 176./3. + 176./3.*Z ) +
            LOGZ*LOGZ*( 55./3 + 55./3.*Z )
           -LOG1MZ*( 152./9. + 956./9.*Z ) +
            LOGZ*( 52./3. + 218./3.*Z )
           -446./27. + 2278./27.*Z
      +
          FACLM*FACLM*(  - 22./3. - 22./3.*Z ) +
          FACLM*RENLR * ( 44./3. + 44./3.*Z ) +
          FACLM   *( LOGZ * (  - 44./3. - 44./3.*Z +
                                 52./3./ZMIN )
                   + DI1MZ * ( 8. + 8.*Z - 16./ZMIN )
                   + ZETA2 * ( 8. + 8.*Z )
                   - 76./9 - 496./9.*Z                   ) +
          RENLR   *( LOGZ * (  - 44./3 - 44./3.*Z +
                                 88./3./ZMIN )
                   + LOG1MZ * ( 88./3. + 88./3.*Z ) ) );
}

// q-\bar{q} C_A terms prior to Harlnader-Kilgore correction:

double qbarq_C_A_preHK(double Z, double muFQ, double muRQ){
 double C_A, C_F, LOGZ, LOG1MZ, DI1MZ;
 C_A = 3.;   C_F = 4./3.;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  DI1MZ = li2(1.-Z);
 return qbarq_C_A(Z,muFQ,muRQ) - 1./16. * C_A * C_F * ( 
      - 16 * Z * LOGZ * LOG1MZ + 8 * Z * LOGZ*LOGZ - 16 * Z * DI1MZ );
}

// the q-\bar{q} C_F terms:

double qbarq_C_F(double Z, double muFQ){
 double C_F, ZMIN, LOGZ, LOG1MZ, DI1MZ, S1MZ, TRI1MZ, FACLM;
 C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);   
 S1MZ = S12(1.-Z);   DI1MZ = li2(1.-Z);   TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);
 return 1./16. * C_F * C_F * ( ( 1. + Z ) *
               ( 48.*S1MZ - 32.*TRI1MZ - 128.*ZETA3 +
                 24.*LOG1MZ*DI1MZ + 24.*LOGZ*DI1MZ +
                 64.*LOG1MZ*ZETA2 -96.*LOGZ*ZETA2
                -64.*LOG1MZ*LOG1MZ*LOG1MZ + 156.*LOG1MZ*LOG1MZ*LOGZ
                -96.*LOG1MZ*LOGZ*LOGZ + 50./3.*LOGZ*LOGZ*LOGZ ) +
          (-64.*S1MZ - 16.*TRI1MZ +
            48.*LOG1MZ*DI1MZ - 48.*LOGZ*DI1MZ +
            128.*LOGZ*ZETA2
           -248.*LOG1MZ*LOG1MZ*LOGZ + 144.*LOG1MZ*LOGZ*LOGZ
           -24.*LOGZ*LOGZ*LOGZ + 112.*LOGZ )/ZMIN +
          DI1MZ*( -24. - 32.*Z ) + ZETA2*( 64. - 64.*Z ) +
          LOG1MZ*LOG1MZ*( - 64. + 64.*Z ) +
          LOG1MZ*LOGZ*( 64. - 112.*Z ) +
          LOGZ*LOGZ*( - 8. + 24.*Z ) +
          LOG1MZ*( 256. + 12.*Z ) +
          LOGZ*( - 104. + 48.*Z ) - 72. + 48.*Z
       +
         FACLM*FACLM*(  LOGZ * ( 24. + 24.*Z - 32./ZMIN )
                  + LOG1MZ * (  - 32. - 32.*Z )
                  - 40. - 8.*Z                          ) +
         FACLM   *( LOGZ*LOGZ * (  - 36. - 36.*Z + 48./ZMIN )
                  + LOGZ*LOG1MZ * ( 144. + 144.*Z - 224./ZMIN )
                   + LOGZ * ( 56. - 24.*Z - 48./ZMIN )
                   + LOG1MZ*LOG1MZ * (  - 96. - 96.*Z )
                   + LOG1MZ * (  - 112. + 16.*Z )
                   + DI1MZ * ( 16. + 16.*Z + 32./ZMIN )
                   + ZETA2 * ( 32. + 32.*Z )
                   + 120. + 16.*Z ) );
}

// q-\bar{q} C_F terms prior to Harlnader-Kilgore correction:

double qbarq_C_F_preHK(double Z, double muFQ){
 double C_F, LOGZ, LOG1MZ, DI1MZ;
 C_F = 4./3.;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  DI1MZ = li2(1.-Z);
 return qbarq_C_F(Z,muFQ) - 1./16. * C_F * C_F * ( 
   - (48+16*Z) * LOGZ * LOG1MZ - 16 * LOGZ + (24+8*Z) * LOGZ*LOGZ
   - (48+16*Z) * DI1MZ ) ;
}

// q-\bar{q} AC interference terms:

double qbarq_AC(double Z, double muFQ){
 double C_A, C_F, ZMIN, LOGZ, LOG1MZ, DI1MZ, S1MZ, TRI1MZ, FACLM;
 C_A = 3.;   C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);   
 S1MZ = S12(1.-Z);  DI1MZ = li2(1.-Z);  TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);  
 return 1./16. * C_F * ( C_F - C_A/2. ) * (
      ( LOGZ*LOGZ*LOGZ * 16./3.  - LOGZ*LOGZ*LOG1MZ * 16.
      + LOGZ*LOGZ * 15.       - LOGZ*LOG1MZ * 24.
      - LOGZ*DI1MZ * 24.   + LOGZ * 24.
      - LOG1MZ*DI1MZ * 32. - DI1MZ * 12.
      + TRI1MZ * 32.      - S1MZ * 72.         )/ZMIN
      + LOGZ*LOGZ*LOGZ * (  - 2. - 2.*Z )
      + LOGZ*LOGZ*LOG1MZ * ( 8. + 8.*Z )
      + LOGZ*LOGZ * ( 4. + 4.*Z )
      + LOGZ*LOG1MZ * (  - 16. - 16.*Z )
      + LOGZ*DI1MZ * ( 16. + 16.*Z )
      + LOGZ * ( 32. - 30.*Z )
      + LOG1MZ*DI1MZ * ( 16. + 16.*Z )
      + LOG1MZ * (  - 64. + 56.*Z )
      + DI1MZ * (  - 20. - 20.*Z )
      + TRI1MZ * (  - 24. - 24.*Z )
      + S1MZ * ( 36. + 36.*Z )
      + 94. - 78.*Z )
    + 1./16. * C_F * ( C_F - C_A/2. ) * FACLM * (
                LOGZ*LOGZ * ( 4. + 4.*Z - 8./ZMIN )
              + LOGZ * (  - 8. - 8.*Z - 12./ZMIN )
              + DI1MZ * ( 8. + 8.*Z - 16./ZMIN )
		- 32. + 28.*Z ) ;
}

// q-\bar{q} BC interference terms:

double qbarq_BC(double Z){
 double C_A, C_F, ZMIN, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ, SMZ, S1MZ, TRIMZ;
 C_A = 3.;   C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z);   
 DIMZ = li2(-Z);  DI1MZ = li2(1.-Z);  SMZ = S12(-Z);  S1MZ = S12(1.-Z);  
 TRIMZ = li3(-Z); 
 return 1./16. * C_F * ( C_F - C_A/2. ) * (
        LOGZ*LOGZ*LOGZ * ( 4./3. + 4./3.*Z*Z + 16./3.*Z )
      + LOGZ*LOGZ*LOG1PZ * ( 20. + 20.*Z*Z + 40.*Z )
      + LOGZ*LOGZ * ( 12. - 30.*Z*Z - 16.*Z )
      + LOGZ*LOG1PZ*LOG1PZ * (  - 24. - 24.*Z*Z - 48.*Z )
      + LOGZ*LOG1PZ * ( 24. + 24.*Z*Z + 48.*Z )
      + LOGZ*DI1MZ * ( 16. + 16.*Z*Z + 48.*Z )
      + LOGZ*DIMZ * ( 24. + 24.*Z*Z + 48.*Z )
      + LOGZ*ZETA2 * ( 8. + 8.*Z*Z + 16.*Z )
      + LOGZ * ( 36. + 44.*Z )
      + LOG1PZ*DIMZ * (  - 48. - 48.*Z*Z - 96.*Z )
      + LOG1PZ*ZETA2 * (  - 24. - 24.*Z*Z - 48.*Z )
      + DI1MZ * ( 36. - 36.*Z*Z )
      + DIMZ * ( 24. + 24.*Z*Z + 48.*Z )
      + ZETA2 * ( 12. + 12.*Z*Z + 24.*Z )
      + TRIMZ * (  - 8. - 8.*Z*Z - 16.*Z )
      + S1MZ * ( 32. + 32.*Z*Z + 96.*Z )
      + SMZ * (  - 48. - 48.*Z*Z - 96.*Z )
      + 54. - 26.*Z*Z - 28.*Z );
}

// q-\bar{q} (or qq) C^2 terms:
// This is "CC_tot" in HvNM's NPB paper
// (or CC_tot_F, in FORTRAN program) -- see HvNMresults.

double qbarq_CC(double Z, double muFQ){ 
  double C_A, C_F, ZMIN, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ, SMZ, S1MZ,
    TRI1MZ, FACLM, PART1, PART2, FRTERM;
 C_A = 3.;   C_F = 4./3.;
 ZMIN = 1.-Z;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z);   
 DIMZ = li2(-Z);  DI1MZ = li2(1.-Z);  SMZ = S12(-Z);  S1MZ = S12(1.-Z);  
 TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);  
 PART1 = ( 1. + Z ) * (
             -8.*TRI1MZ + 12.*S1MZ + 8.*LOG1MZ*DI1MZ +
            2.*LOGZ*DI1MZ - 4.*LOGZ*ZETA2 + 3./2.*LOGZ*LOGZ*LOGZ
           -4.*LOGZ*LOGZ*LOG1MZ + 4.*LOGZ*LOG1MZ*LOG1MZ );
 PART2 = DI1MZ * ( 13. + 8./3.*Z*Z + 5.*Z +
                     16./3./Z )
         + ZETA2 * ( - 2. + 8./3.*Z*Z + 2.*Z
                     - 8./3./Z )
         + LOG1MZ*LOG1MZ * ( 2. - 8./3.*Z*Z - 2.*Z +
                         8./3./Z )
         + LOGZ*LOG1MZ * ( 6. + 8.*Z*Z + 12.*Z )
         + LOGZ*LOGZ * ( - 5./4. - 10./3.*Z*Z
                       - 25./4.*Z )
         + LOG1MZ * ( - 26./3. - 44./9.*Z*Z +
                        26./3.*Z + 44./9./Z )
         + LOGZ * ( 115./6. + 10./9.*Z*Z
                  - 8./3.*Z )
         + 593./36. + 703./108.*Z*Z
         - 433./18.*Z + 29./27./Z ;
 FRTERM = FACLM*FACLM * ( LOGZ * ( 1. + Z ) +
                           1./2. - 2./3.*Z*Z - 1./2.*Z +
                           2./3./Z ) +
              FACLM * ( DI1MZ * ( 4. + 4.*Z ) +
                        LOGZ*LOGZ * ( - 2. - 2.*Z ) +
                        LOGZ*LOG1MZ * ( 4. + 4.*Z ) +
                        LOGZ * ( 3. + 4.*Z*Z + 6.*Z ) +
                        LOG1MZ * (  2. - 8./3.*Z*Z
                                  - 2.*Z + 8./3./Z )
                       -13./3. - 22./9.*Z*Z +
                        13./3.*Z + 22./9./Z ) ;
 return 1./16. * ( C_A*C_A - 1. )/ C_A * ( PART1 + PART2 + FRTERM );
}
 
// q-\bar{q} (or qq) CD (vector) interference terms.
// This is "CD_V_tot" in HvNM's NPB paper
// (or -2*CD_V_tot_F, in FORTRAN program) -- see HvNMresults.

double qbarq_CD_V(double Z){ 
 double C_F, T_F, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ, SMZ, S1MZ, TRIMZ, TRI1MZ;
 C_F = 4./3.;  T_F = 1./2.;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z);   
 DIMZ = li2(-Z);  DI1MZ = li2(1.-Z);  SMZ = S12(-Z);  S1MZ = S12(1.-Z);  
 TRIMZ = li3(-Z); TRI1MZ = li3(1.-Z);
 return 1./16. * C_F * T_F * (
  (2. + Z + 2./Z) * ( 32. * S1MZ - 96. * SMZ - 48. * LOG1PZ*LOG1PZ * LOGZ
        - 48. * ZETA2 * LOG1PZ + 40. * LOGZ*LOGZ * LOG1PZ
        - 96. * DIMZ * LOG1PZ )
 + (1.+Z) * ( 80. * DIMZ + 80. * LOGZ * LOG1PZ + 40. * ZETA2 )
 + 8. * (-6. + 3. * Z + 4./Z) * TRI1MZ 
 - 16. * (-10. + 3. * Z + 10./Z) * TRIMZ 
 - 24. * (-6. + Z + 4./Z) * ZETA3 + 8. * (10.-Z) * DI1MZ * LOGZ
 - 16./3. * Z * LOGZ*LOGZ*LOGZ
 + 32. * (2.*Z + 5./Z) * DIMZ * LOGZ + 8. * (10.+Z) * ZETA2 * LOGZ
 + 8. * (5.-4.*Z) * DI1MZ - 52. * Z * LOGZ*LOGZ - 16. * (5.+4.*Z) * LOGZ
  - 160. * (1.-Z) );
}

// The total q-\bar{q} "hard" cross section, corresponding to "resqq_total":

double qbarq_hard_NNLO(double z, double NF, double muFQ, double muRQ){
 return NF * qbarq_NF(z,muFQ,muRQ) 
      + qbarq_C_A(z,muFQ,muRQ) + qbarq_C_F(z,muFQ) 
      + 2. * qbarq_AC(z,muFQ) + 2. * qbarq_BC(z) 
      + 2. * qbarq_CC(z,muFQ) + qbarq_CD_V(z);
}

// In HvNM's notation, just the nonsinglet (NS) hard terms:

double qbarq_NS_hard_NNLO(double z, double NF, double muFQ, double muRQ){
 return NF * qbarq_NF(z,muFQ,muRQ) 
     + qbarq_C_A(z,muFQ,muRQ) + qbarq_C_F(z,muFQ) + 2. * qbarq_AC(z,muFQ);
}

// In HvNM's notation, just the nonsinglet (NS) hard terms,
// PRIOR TO HK CORRECTION.

double qbarq_NS_hard_NNLO_preHK(double z, double NF, double muFQ, double muRQ){
 return NF * qbarq_NF(z,muFQ,muRQ) 
     + qbarq_C_A_preHK(z,muFQ,muRQ) + qbarq_C_F_preHK(z,muFQ) 
     + 2. * qbarq_AC(z,muFQ);
}

// TEST function:  Integral of just the "hard" part of the NO-NF terms.

double diffnoNF(double z){ 
  double ln2, lnz, lnz2, ln1mz, ln1mz2, lnhf1pz, lnhf1pz2, 
         z2, z3, z4, z5, z6;
  ln2 = log(2.);  
  lnz = log(z);  lnz2 = lnz*lnz;  
  ln1mz = log(1.-z);  ln1mz2 = ln1mz*ln1mz; 
  lnhf1pz = log((1.+z)/2.);  lnhf1pz2 = lnhf1pz*lnhf1pz; 
  z2 = z*z;  z3 = z2*z;  z4 = z2*z2;  z5 = z3*z2;  z6 = z3*z3;
  return 
-1./27*(280*z4-57*z3+97*z2+13*z+z5+24)*lnz*lnz2/(z+1)/z/(-1+z)
+1./9*(-36*z+169*z5+6*z4+115*z3-24-6*z2)*lnz2*ln1mz/z2/(-1+z)/(z+1)
-28./9*(z+1)*ln1mz*li2(-z)
-2./9*(39*z4-7*z2-48*z+48)*ln1mz2*lnz/z2/(-1+z)
-1./27*(7*z2+73*z-15*z4+24*z3+3*z5+24)*6.*ZETA2*lnz/(z+1)/z/(-1+z)
-2./9*lnz*li2(-z)*(3*z4-37*z3+23*z2-33*z+30)/z/(-1+z)
+16./9*lnz*li2((1.+z)/2.)*(z2+1)/(-1+z)
-1./9*(-12*z+4*z5+13*z4-162*z3+24+49*z2)*lnz*li2(z)/(z+1)/z/(-1+z)
-1./9*(z-16*z3+3*z4-12*z2+18)/z/(-1+z)*ln1mz*lnhf1pz2
+1./27*(-62*z3-31*z+3*z4-12*z2+18)/z/(-1+z)*lnhf1pz*lnhf1pz2
+8./3*(z2+1)/(-1+z)*lnz*lnhf1pz2
-1./9*(5*z5-89*z4-38*z3-39*z2+31*z+30)/(z+1)/z/(-1+z)*lnz2*lnhf1pz
+1./27*(42*z+18+3*z5-15*z4+2*z3-20*z2)/(z+1)/z/(-1+z)*6.*ZETA2*lnhf1pz
+4./9*(8*z+15*z3-8-z2)/(-1+z)/(z+1)*li2(z)*lnhf1pz
-4./9*(35*z3+19*z2+21*z+5)/(-1+z)/(z+1)*lnz*ln1mz*lnhf1pz
+16./9*(z2+1)/(-1+z)*ln2*lnz*lnhf1pz
-2.*(z2+1)/(-1+z)*li2((1.+z)/2.)*lnhf1pz
+64./9*(z2+1)/(-1+z)*ln1mz2*lnhf1pz
-28./9*z2/(-1+z)*li2(-z)*lnhf1pz
+8./9*ln1mz*6.*ZETA2*(z3+6*z2+5*z+4)/z/(z+1)
+2./9*ln1mz*li2(z)*(17*z4-62*z3+17*z2+48*z+48)/(z+1)/z2
+4./9*(-33*z5-54*z4-z3+24*z2+18*z-12+2*z6)*li3(z)/z2/(-1+z)/(z+1)
-2./9*(3*z4-9*z3-12*z2+8*z+18)*li3((1.-z)/2.)/z/(-1+z)
-2./9*(-18*z3-z-12*z2+3*z4+18)*li3((1.-z)/(z+1.))/z/(-1+z)
-28./9*(z+1)*ln1mz*ln2*lnz
-1./9*lnz2*ln2*(5*z4-46*z3+8*z2+z+30)/z/(-1+z)
+2./9*(2*z3-51*z+z4+38*z2+30)*li3(-z)/z/(-1+z)
-2./9*(-4*z3-z+3*z4-12*z2+18)*li3((1.+z)/2.)/z/(-1+z)
+1./9*(z-16*z3+3*z4-12*z2+18)*ln2*lnhf1pz2/z/(-1+z)
-1./9*(3*z4-11*z3-12*z2-8*z+18)*li3((1.-z)*(z+1.))/z/(-1+z)
+1./9*(-55*z5-238*z4+69*z3+252*z2+12*z6+60*z-144)
  *li3(1.-z)/z2/(-1+z)/(z+1)
+4./9*(8+15*z2)*li3(2.*z/(z+1.))/(-1+z)
-1./18*(-288-402*z+43*z5-176*z4+339*z3+z6+251*z2)*ZETA3/z2/(-1+z)/(z+1)
//
-1./9.*(15.*z2-11.*z3-21.*z-18.+2.*z4)/z/(-1.+z)*lnhf1pz2
-4./9./z*(3.*z2+12.*z+10.)*ln1mz*lnhf1pz
-2./9.*lnhf1pz*lnz*(3.*z5-19.*z4-47.*z3+13.*z2+58.*z+16.)/z/(-1.+z)/(z+1.)
+1./54.*(10.*z5-90.*z4-28.*z3+69.*z2-97.*z+46.)*6.*ZETA2/z/(-1.+z)/(z+1.)
+1./18.*lnz2*(43.*z5-185.*z4+35.*z3-118.*z2+126.*z+33.)/z/(-1.+z)/(z+1.)
+64./9.*(-1.+z)*ln1mz2
-1./9.*li2(z)*(9.*z6-50.*z5-39.*z4+95.*z3-40.*z2-8.*z-33.)/z2/(-1.+z)/(z+1.)
-2./9.*li2((1.-z)/2.)*(22.*z2-22.*z3-27.*z-8.+2.*z4)/z/(-1.+z)
-1./9.*ln1mz*lnz*(17.*z6-22.*z5-57.*z4-73.*z3+110.*z2-8.*z-33.)
    /z2/(-1.+z)/(z+1.)
-2./9.*(3.*z4-11.*z3-22.*z2+26.*z+16.)*(ln2*lnz+li2(-z))/(-1.+z)/z
//
-1./54*(220*z4-734*z3+79*z2+754*z-559)*lnz/(-1+z)/(z+1)
+4./27*(-1+z)*(11*z2-118*z+14)*ln1mz/z
-2./27*(27*z3+33*z2-10*z-20)*lnhf1pz/(-1+z)/z
//
+1./108/z*(403*z3+1688*z2-2071*z-140) ;
}

// TEST function:  Integral of the "lonely" 1/(1-z)_+ parts of the NF terms.

double lonely(double z, double NF){
  return 2. * ( 
      - 32./9. * 6.*ZETA2 * D1(1.,1.-z)
      + ( 32. * ZETA3 - ( -22./3. + 4./9. * NF ) * ZETA2 ) * D1(0.,1.-z) );
}

// to be evaluated at z = x2:

double lonely_e(double z, double NF){
  return 2. * ( 
      - 32./9. * 6.*ZETA2 * D1e(1.,1.-z)
      + ( 32. * ZETA3 - ( -22./3. + 4./9. * NF ) * ZETA2 ) * D1e(0.,1.-z) );
}

//=========== quark-gluon terms ============================

double qg_C_A(double Z, double muFQ, double muRQ){ 
 double C_A, ZCOM, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ,
    S1MZ, SMZ, TRIMZ, TRI1MZ, TRIPCO, TRIMCO, FACLM, RENLR, DELCA, HELPCA;
 C_A = 3.; LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z); 
 S1MZ = S12(1.-Z);  SMZ = S12(-Z); 
 DIMZ = li2(-Z);   DI1MZ = li2(1.-Z);  
 TRIMZ = li3(-Z);  TRI1MZ = li3(1.-Z);
 ZCOM = (1.-Z)/(1.+Z);
 TRIPCO = li3(ZCOM);  TRIMCO = li3(-ZCOM);
 FACLM = -2.*log(muFQ);   RENLR = -2.*log(muRQ);
//
   DELCA =
        LOGZ*LOGZ*LOGZ * (  - 3 - 20/3.*Z )
      + LOGZ*LOGZ*LOG1MZ * ( 6 - 4*Z*Z + 28*Z )
      + LOGZ*LOGZ*LOG1PZ * (  - 6 - 12*Z*Z - 12*Z )
      + LOGZ*LOGZ * ( 5/2. + 173/3.*Z*Z + 2*Z )
      + LOGZ*LOG1MZ*LOG1MZ * (  - 2 + 12*Z*Z - 44*Z )
      + LOGZ*LOG1MZ*LOG1PZ * ( 8 + 16*Z*Z + 16*Z )
      + LOGZ*LOG1MZ * (  - 10 - 130*Z*Z + 16*Z )
      + LOGZ*LOG1PZ * (  - 4 - 16*Z*Z - 20*Z )
      + LOGZ*DI1MZ * ( 8*Z*Z - 28*Z )
      + LOGZ*DIMZ * (  - 8 - 16*Z*Z - 16*Z )
      + LOGZ*ZETA2 * (  - 16*Z*Z + 40*Z )
      + LOGZ * (  - 118/3. + 457/9.*Z*Z + 16/3.*Z )
      + LOG1MZ*LOG1MZ*LOG1MZ * (  - 13/3. - 26/3.*Z*Z + 26/3.*Z )
      + LOG1MZ*LOG1MZ * (  - 4 + 154/3.*Z*Z - 42*Z
                       - 16/3./Z )
      + LOG1MZ*DI1MZ * (  - 28 - 20*Z*Z - 40*Z )
      + LOG1MZ*DIMZ * ( 8 + 16*Z*Z + 16*Z )
      + LOG1MZ*ZETA2 * ( 16 + 32*Z*Z - 16*Z )
      + LOG1MZ * ( 70/3. - 74/9.*Z*Z
                  - 25/3.*Z - 88/9./Z )
      + DI1MZ * (  - 22 - 88/3.*Z*Z - 64*Z - 32/3./Z )
      + DIMZ * (  - 4 - 16*Z*Z - 20*Z )
      + ZETA2 * (  - 10 - 214/3.*Z*Z + 56*Z
                   + 16/3./Z )
      + TRI1MZ * ( 30 + 24*Z*Z + 68*Z )
      + TRIMZ * ( 4 + 8*Z*Z + 8*Z )
      + TRIMCO * ( 8 + 16*Z*Z + 16*Z )
      + TRIPCO * (  - 8. - 16.*Z*Z - 16.*Z )
      + S1MZ * (  - 36 - 16*Z*Z - 64*Z )
      + ZETA3 * ( 2 + 4*Z*Z + 8*Z )
      - 539/18. - 1837/54.*Z*Z + 613/9.*Z
      - 58/27./Z ;
//
 HELPCA = 
          FACLM*FACLM*( LOGZ * ( 2 + 8*Z )
                  + LOG1MZ * ( 2 + 4*Z*Z - 4*Z )
                  + ( 14 - 9*Z*Z + 2*Z + 4./Z ) / 3. ) +
         FACLM*RENLR*( -11/3. - 22/3.*Z*Z + 22/3.*Z ) +
         FACLM   *( LOGZ*LOGZ * (  - 4 - 12*Z )
                  + LOGZ*LOG1MZ * ( 4 - 8*Z*Z + 40*Z )
                  + LOGZ*LOG1PZ * (  - 4 - 8*Z*Z - 8*Z )
                  + LOGZ * ( 7 + 146*Z*Z + 10*Z ) / 3.
                  + LOG1MZ*LOG1MZ * ( 6 + 12*Z*Z - 12*Z )
                  + LOG1MZ * ( 40/3. - 98/3.*Z*Z +
                               64/3.*Z + 16/3./Z )
                  + DI1MZ * ( 12 + 8*Z*Z + 24*Z )
                  + DIMZ * (  - 4 - 8*Z*Z - 8*Z )
                  + ZETA2 * (  - 8 - 16*Z*Z + 8*Z )
                  - 47/6. - 85/18.*Z*Z + 29/3.*Z
                  + 44/9./Z  ) +
         RENLR   *( LOGZ * ( 11 + 22*Z*Z - 22*Z ) / 3.
                  + LOG1MZ * ( -22 - 44*Z*Z + 44*Z ) / 3.
                  - 11/6. + 77/6.*Z*Z - 11*Z        ) ;
//
 return 1./16. * ( - C_A * DELCA + C_A * HELPCA );
}

// qg C_A terms prior to Harlnader-Kilgore correction:

double qg_C_A_preHK(double Z, double muFQ, double muRQ){
 double C_A, T_F, LOGZ, LOG1MZ, DI1MZ;
 C_A = 3.;   T_F = 1./2.;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  DI1MZ = li2(1.-Z);
 return qg_C_A(Z,muFQ,muRQ) - 1./16. * C_A * T_F * ( 
      8 * Z * LOGZ * LOG1MZ - 8 * Z * LOGZ 
    - 4 * Z * LOGZ*LOGZ + 8 * Z * DI1MZ );
}

double qg_C_F(double Z, double muFQ, double muRQ){ 
 double C_F, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ,
    S1MZ, SMZ, TRIMZ, TRI1MZ, FACLM, DELCF, HELPCF;
 C_F = 4./3.;  LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z); 
 S1MZ = S12(1.-Z);  SMZ = S12(-Z); 
 DIMZ = li2(-Z);   DI1MZ = li2(1.-Z);  
 TRIMZ = li3(-Z);  TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);
//
    DELCF =
        LOGZ*LOGZ*LOGZ * ( 17/6. + 26/3.*Z*Z - 17/3.*Z )
      + LOGZ*LOGZ*LOG1MZ * (  - 12 - 40*Z*Z + 24*Z )
      + LOGZ*LOGZ * ( 11/4. + Z*Z - 15*Z )
      + LOGZ*LOG1MZ*LOG1MZ * ( 21 + 66*Z*Z - 42*Z )
      + LOGZ*LOG1MZ * (  - 14 - 96*Z*Z + 96*Z )
      + LOGZ*LOG1PZ * (  - 8 + 24*Z*Z + 16*Z )
      + LOGZ*DI1MZ * (  - 2 + 4*Z )
      + LOGZ*DIMZ * ( 8 + 16*Z*Z - 16*Z )
      + LOGZ*ZETA2 * (  - 12 - 48*Z*Z + 24*Z )
      + LOGZ * ( 31/2. + 87*Z*Z - 201/2.*Z )
      + LOG1MZ*LOG1MZ*LOG1MZ * (  - 35/3. - 70/3.*Z*Z + 70/3.*Z )
      + LOG1MZ*LOG1MZ * ( 23 + 63*Z*Z - 80*Z )
      + LOG1MZ*DI1MZ * ( 6 + 52*Z*Z - 12*Z )
      + LOG1MZ*ZETA2 * ( 8 + 16*Z*Z - 16*Z )
      + LOG1MZ * (  - 26 - 88*Z*Z + 135*Z )
      + DI1MZ * (  9 - 40*Z*Z + 24*Z )
      + DIMZ * (  - 8 + 24*Z*Z + 16*Z )
      + ZETA2 * (  - 10 + 24*Z*Z - 4*Z )
      + TRI1MZ * ( 2 - 36*Z*Z - 4*Z )
      + TRIMZ * (  - 16 - 32*Z*Z + 32*Z )
      + S1MZ * ( 22 + 68*Z*Z - 44*Z )
      + ZETA3 * (  - 50 - 100*Z*Z + 100*Z )
      + 157/4. + 305/4.*Z*Z - 221/2.*Z ;
//
      HELPCF =
         FACLM*FACLM*( LOGZ * (  - 3 - 12*Z*Z + 6*Z )
                  + LOG1MZ * ( 6 + 12*Z*Z - 12*Z )
                  - 3/2. + 6*Z                            ) +
         FACLM   *( LOGZ*LOGZ * ( 4 + 16*Z*Z - 8*Z )
                  + LOGZ*LOG1MZ * (  - 20 - 64*Z*Z + 40*Z )
                  + LOGZ * ( 5 + 46*Z*Z - 40*Z )
                  + LOG1MZ*LOG1MZ * ( 18 + 36*Z*Z - 36*Z )
                  + LOG1MZ * (  - 16 - 46*Z*Z + 68*Z )
                  + DI1MZ * (  - 24*Z*Z )
                  + ZETA2 * (  - 4 - 8*Z*Z + 8*Z )
                  + 12 + 11*Z*Z - 34*Z                  );
//
 return 1./16. * ( -C_F * DELCF + C_F * HELPCF );
}

// qg C_F terms prior to Harlnader-Kilgore correction:

double qg_C_F_preHK(double Z, double muFQ, double muRQ){
 double C_F, T_F, LOGZ, LOG1MZ, DI1MZ;
 C_F = 4./3.;   T_F = 1./2.;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  DI1MZ = li2(1.-Z);
 return qg_C_F(Z,muFQ,muRQ) - 1./16. * C_F * T_F * ( 
    - (24-8*Z) * LOGZ * LOG1MZ - 24 * (1-Z) * LOG1MZ
    + (28-44*Z) * LOGZ + (12-4*Z) * LOGZ*LOGZ
    - (24-8*Z) * DI1MZ + 12 * (1-Z) );
}

double qg_NF(double Z, double muFQ, double muRQ){ 
 double LOGZ, LOG1MZ, FACLM, RENLR;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);
 FACLM = -2.*log(muFQ);   RENLR = -2.*log(muRQ/muFQ);
 return 1/16. * RENLR * 2/3. * ( 
        ( FACLM - LOGZ + 2.*LOG1MZ ) * ( 1. + 2.*Z*Z - 2.*Z )
	+ 1/2. * ( 1. - 7.*Z*Z + 6.*Z ) );
}

//=========== gluon-gluon terms ============================

double gg_C_A(double Z){ 
 double C_A, LOGZ, LOG1PZ, DIMZ, S1MZ, SMZ, TRIMZ, TRI1MZ;
 C_A = 3.;   LOGZ = log(Z);  LOG1PZ = log(1.+Z); 
 S1MZ = S12(1.-Z);  SMZ = S12(-Z); 
 DIMZ = li2(-Z);  TRIMZ = li3(-Z);  TRI1MZ = li3(1-Z);
 return 1./16. * C_A*C_A / ( C_A*C_A - 1. ) *
           ( S1MZ * ( -8*Z*Z + 16*Z - 8 ) +
             SMZ * ( 16*Z*Z + 32*Z + 16 ) +
             TRIMZ * ( 24*Z*Z + 48*Z + 24 ) +
             ZETA3 * ( 16*Z*Z + 32*Z + 16 ) +
             DIMZ * LOGZ * ( -24*Z*Z - 48*Z - 24 ) +
             DIMZ * LOG1PZ * ( 16*Z*Z + 32*Z + 16 ) +
             DIMZ * ( 16*Z*Z + 32*Z +16 )/3. +
             ZETA2 * LOG1PZ * ( 8*Z*Z + 16*Z + 8 ) +
             ZETA2 * ( 8*Z*Z + 16*Z +8 )/3. +
             LOGZ*LOGZ * LOG1PZ * ( -12*Z*Z - 24*Z - 12 ) +
             LOGZ*LOGZ * ( 50*Z*Z + 4*Z - 4 )/3. +
             LOGZ * LOG1PZ*LOG1PZ * ( 8*Z*Z + 16*Z + 8 ) +
             LOGZ * LOG1PZ * ( 16*Z*Z + 32*Z + 16)/3. +
             LOGZ * ( -50*Z*Z - 76/3.*Z - 4 ) +
             191/3.*Z*Z - 48*Z - 47/3. );
}

double gg_C_F(double Z, double muFQ){ 
 double LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ,
    S1MZ, SMZ, TRIMZ, TRI1MZ, FACLM, DELCF1, DELCF2, FRTERM;
 LOGZ = log(Z);   LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z); 
 S1MZ = S12(1.-Z);  SMZ = S12(-Z); 
 DIMZ = li2(-Z);   DI1MZ = li2(1.-Z);  
 TRIMZ = li3(-Z);  TRI1MZ = li3(1.-Z);
 FACLM = -2.*log(muFQ);
//
 DELCF1 = TRI1MZ * ( 16 + 64*Z + 64*Z*Z ) +
            TRIMZ * ( - 8 - 16*Z + 8*Z*Z ) +
            S1MZ * ( - 8 - 80*Z - 56*Z*Z ) +
            SMZ * ( - 16 - 32*Z - 16*Z*Z ) +
            ZETA3 * ( - 4 - 8*Z + 8*Z*Z ) +
            LOG1MZ*DI1MZ * ( - 16 - 64*Z - 64*Z*Z ) +
            LOGZ*DI1MZ * ( - 4 - 16*Z - 16*Z*Z ) +
            LOG1PZ*DIMZ * ( - 16 - 32*Z - 16*Z*Z ) +
            LOGZ*DIMZ * ( 16 + 32*Z + 8*Z*Z ) +
            LOGZ*ZETA2 * ( 12 + 40*Z + 40*Z*Z ) +
            LOG1PZ*ZETA2 * ( - 8 - 16*Z - 8*Z*Z ) +
            LOG1MZ*LOG1MZ*LOGZ * ( - 8 - 32*Z - 32*Z*Z ) +
            LOG1MZ*LOGZ*LOGZ * ( 4 + 16*Z + 16*Z*Z ) +
            LOGZ*LOGZ*LOGZ * ( - 2 - 16/3.*Z
                        - 16/3.*Z*Z ) +
            LOGZ*LOGZ*LOG1PZ * ( 12 + 24*Z + 12*Z*Z ) +
            LOGZ*LOG1PZ*LOG1PZ * ( - 8 - 16*Z - 8*Z*Z ) ;
//
   DELCF2 = DI1MZ * ( - 20 - 16*Z + 56*Z*Z ) +
            DIMZ * ( 8 + 8*Z ) +
            ZETA2 * ( 20 + 36*Z - 48*Z*Z ) +
            LOG1MZ*LOG1MZ * ( - 16 - 32*Z + 48*Z*Z ) +
            LOG1MZ*LOGZ * ( 4 + 32*Z - 16*Z*Z ) +
            LOGZ*LOGZ * ( - 6 - 14*Z - 8*Z*Z ) +
            LOGZ*LOG1PZ * ( 8 + 8*Z ) +
            LOG1MZ * ( 14 + 120*Z - 134*Z*Z ) +
            LOGZ * ( - 23 - 64*Z + 105*Z*Z )
            - 32 - 66*Z + 98*Z*Z ;
//
   FRTERM = ( LOGZ * ( - 2 - 8*Z - 8*Z*Z )
                - 4 - 8*Z + 12*Z*Z          ) * FACLM*FACLM +
              ( DI1MZ * ( - 8 - 32*Z - 32*Z*Z ) +
                LOG1MZ*LOGZ * ( - 8 - 32*Z - 32*Z*Z ) +
                LOGZ*LOGZ * ( 2 + 8*Z + 8*Z*Z ) +
                LOG1MZ * ( - 16 - 32*Z + 48*Z*Z ) +
                LOGZ * ( 2 + 16*Z - 8*Z*Z ) +
                7 + 60*Z - 67*Z*Z                ) * FACLM ;
//
 return 1./16. * ( DELCF1 + DELCF2 + FRTERM );
}

//=========== quark-quark terms ============================

// Interference exchange terms, for identical final state quarks.

// This is qq_CE_tot = qq_CE_tot_F(z): 

double qq_CE_tot(double Z, double muFQ){
 double C_A, C_F, ZPLUS, ZCOM, LOGZ, LOG1MZ, LOG1PZ, DIMZ, DI1MZ,
    S1MZ, SMZ, TRIMZ, TRI1MZ, TRIPCO, TRIMCO, FACLM, DELCAF, FRTERM;
 C_A = 3.;   C_F = 4./3.;  
 ZPLUS = 1.+Z;   ZCOM = (1.-Z)/(1.+Z);
 LOGZ = log(Z);  LOG1MZ = log(1.-Z);  LOG1PZ = log(1.+Z); 
 S1MZ = S12(1.-Z);  SMZ = S12(-Z); 
 DI1MZ = li2(1.-Z);  DIMZ = li2(-Z);
 TRIMZ = li3(-Z);  TRI1MZ = li3(1.-Z);
 TRIPCO = li3(ZCOM);  TRIMCO = li3(-ZCOM);
 FACLM = -2.*log(muFQ);  
 DELCAF = 1./16. * C_F * ( C_F - C_A/2. ) * ( 
   LOGZ*LOGZ*LOGZ * ( 2 - 2*Z - 16./3./ZPLUS )
      + LOGZ*LOGZ*LOG1MZ * (  - 8 + 8*Z + 16/ZPLUS )
      + LOGZ*LOGZ*LOG1PZ * (  - 24 + 24*Z + 56/ZPLUS )
      + LOGZ*LOGZ * (  - 4 - 12*Z )
      + LOGZ*LOG1MZ*LOG1PZ * ( 32 - 32*Z - 64/ZPLUS )
      + LOGZ*LOG1MZ * ( 16 + 16*Z )
      + LOGZ*LOG1PZ*LOG1PZ * (  - 16/ZPLUS )
      + LOGZ*LOG1PZ * ( 8 + 8*Z )
      + LOGZ*DI1MZ * (  - 24 + 24*Z + 48/ZPLUS )
      + LOGZ*DIMZ * (  - 32 + 32*Z + 64/ZPLUS )
      + LOGZ*ZETA2 * (  - 8 + 8*Z + 24/ZPLUS )
      + LOGZ * (  - 18 + 14*Z )
      + LOG1MZ*DIMZ * ( 32 - 32*Z - 64/ZPLUS )
      + LOG1MZ*ZETA2 * ( 16 - 16*Z - 32/ZPLUS )
      + LOG1MZ * ( 32 - 32*Z )
      + LOG1PZ*DIMZ * (  - 32/ZPLUS )
      + LOG1PZ*ZETA2 * (  - 16/ZPLUS )
//
      + DI1MZ * ( 24 + 8*Z )
      + DIMZ * ( 8 + 8*Z )
      + ZETA2 * ( 4 + 4*Z )
      + TRI1MZ * ( 32 - 32*Z - 64/ZPLUS )
      + TRIMZ * ( 16 - 16*Z - 16/ZPLUS )
      + TRIMCO * ( 32 - 32*Z - 64/ZPLUS )
      + TRIPCO * (  - 32 + 32*Z + 64/ZPLUS )
      + S1MZ * (  - 32 + 32*Z + 64/ZPLUS )
      + SMZ * (  - 32/ZPLUS )
      + ZETA3 * ( 12 - 12*Z - 8/ZPLUS )
      - 34 + 34*Z ) ;
//
   FRTERM = 1./16. * C_F * ( C_F - C_A/2. ) * FACLM * (
              LOGZ*LOGZ * (  - 4 + 4*Z + 8/ZPLUS )
            + LOGZ*LOG1PZ * ( 16 - 16*Z - 32/ZPLUS )
            + LOGZ * ( 8 + 8*Z )
            + DIMZ * ( 16 - 16*Z - 32/ZPLUS )
            + ZETA2 * ( 8 - 8*Z - 16/ZPLUS )
            + 16 - 16*Z                 ) ;
 return DELCAF + FRTERM;
}

// This is qq_CF_tot = 2 * qq_CF_tot_F(z): 

double qq_CF(double Z){
 double C_A, C_F, LOGZ, DI1MZ, S1MZ, TRI1MZ;
 C_A = 3.;   C_F = 4./3.;  
 LOGZ = log(Z);   S1MZ = S12(1.-Z); 
 DI1MZ = li2(1.-Z);  TRI1MZ = li3(1.-Z);
 return 2. * 1./16. * C_F * ( C_F - C_A/2. ) * (
        LOGZ*LOGZ*LOGZ * (  - 4./3. - 4./3.*Z*Z + 8./3.*Z )
      + LOGZ*LOGZ * (  - 6 - 6*Z*Z + 12*Z )
      + LOGZ*DI1MZ * (  - 8 - 8*Z*Z + 16*Z )
      + LOGZ * (  - 14 + 12*Z )
      + DI1MZ * (  - 12 - 12*Z*Z + 24*Z )
      + TRI1MZ * ( 8 + 8*Z*Z - 16*Z )
      + S1MZ * (  - 8 - 8*Z*Z + 16*Z )
      - 15 - 13*Z*Z + 28*Z ) ;
}

// Axial vector contributions:

// For AB, we insert a T_F with respect to the HvNM NPB formula !!!
// This is AB_ax_tot = AB_ax_tot_F(z):

double AB_ax(double Z){
 double C_F, LOGZ, ZMIN;
 C_F = 4./3.;  
 LOGZ = log(Z);  ZMIN = 1.-Z; 
 return 1./16. * C_F * ( LOGZ * (  - 8 + 8*Z + 16/ZMIN )
      + 24 - 8*Z ) ;
}

// This is CD_ax_tot = 2 * CD_ax_tot_F(z):

double CD_ax(double Z){
 double C_F, LOGZ, LOG1PZ, DIMZ, DI1MZ, SMZ, S1MZ, TRIMZ, TRI1MZ, temp;
 C_F = 4./3.;  
 LOGZ = log(Z);   LOG1PZ = log(1.+Z);  
 SMZ = S12(-Z);  S1MZ = S12(1.-Z);
 DIMZ = li2(-Z);  DI1MZ = li2(1.-Z);  
 TRIMZ = li3(-Z);  TRI1MZ = li3(1.-Z);
 temp = LOGZ*LOGZ*LOGZ * (  - 4./3.*Z )
      + LOGZ*LOGZ*LOG1PZ * ( 20 + 10*Z )
      + LOGZ*LOGZ * (  - Z )
      + LOGZ*LOG1PZ*LOG1PZ * (  - 24 - 12*Z )
      + LOGZ*LOG1PZ * ( 4 + 4*Z )
      + LOGZ*DI1MZ * ( 4 + 6*Z )
      + LOGZ*DIMZ * ( 32 )
      + LOGZ*ZETA2 * ( 4 + 10*Z )  - LOGZ * 4
      + LOG1PZ*DIMZ * (  - 48 - 24*Z )
      + LOG1PZ*ZETA2 * (  - 24 - 12*Z )
      + DI1MZ * ( 2 )
      + DIMZ * ( 4 + 4*Z )
      + ZETA2 * ( 2 + 2*Z )
      + TRI1MZ * ( 4 - 2*Z )
      + TRIMZ * (  - 24 + 20*Z )
      + S1MZ * ( 16 + 8*Z )
      + SMZ * (  - 48 - 24*Z )
      + ZETA3 * (  - 12 + 18*Z ) - 8 + 8*Z;
 return 2. * 1./16. * C_F * temp;
}
