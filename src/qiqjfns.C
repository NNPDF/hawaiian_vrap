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

#include "hardfns.h"

// "11" terms:
// main "hard" function:

double qq11(double y, double z){
  double u,t,r,Z2,z2,z3,z4, u2,u3,u4,u5,u6,u7, up12,
    ln2,lnz, lnr,lnu,lnopu,lnhfopu,lnumo, lnrpt,lnrmt, lnhfopz,lnoptr,lnomtr, 
    lnrmo, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;   z3 = z2*z;  z4 = z2*z2;  
 u2 = u*u;   u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;  u6 = u3*u3;  u7 = u4*u3;
 up12 = (1.+u)*(1.+u);
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnhfopu = log((1.+u)/2.);  lnumo = log(std::fabs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);   lnhfopz = log((1.+z)/2.); 
 lnoptr = log(1.+t*r);  lnomtr = log(1.-t*r);  lnrmo = log(std::fabs(r-1.));  
 temp = 
// q^2 terms (mostly):
// J's (only J21 here!):
- 8./3. * z2 * (1-z)/u/(z*up12-4*u2)
     * ( z2 - z*(2-z)*u - z*(6+z)*u2 + (6+2*z-z2)*u3 )
    * ( J21(y,z) 
	- (1+z) * (1-t*r)/(1+u)/r/(r+t)/(1-z)/(1-z)
        * lnz * lnhfopz
// from q^1 terms:
      + 1./z/(1+u) * ( log((r+t)/(1.+u)) + 2*u/(1-u) * lnz 
    - 1./2. * (1+u)/(1-u) * (1+t/r) * lnu ) )
// polylog's
+ 4./3. * z2 * (1-z) * (1+u)/(u-z) 
    * ( li2(z) - li2(-z) - lnz * log(1.+z) - 3./2. * Z2 )
- 4./3. * z * (1-z)/u3/(u-z) * ( 11*z3 - z2*(17-7*z)*u 
       + 2*z*(1-z)*(4-z)*u2 - (2-2*z+3*z2)*u3 )
    * ( li2(-u) + li2(t*r) + li2(-t/r)
      + li2((1.-t*r)/(1.+u)) )
- 8./3. * z3 * (1-z) * (1+u)/u/(u-z)
    * ( li2(t*r) - Z2 )
- 2./3. * z * (1-z)/u/(u-z) * ( 4*z2*(1+u) -(2+z)*u2 - u3 )
    * ( li2((1.-t*r)/(1+u)) - 1./2. * lnz * lnrmt )
// pure ln's
- 1./3. * z * (1-z)/u3/up12/(u-z) * ( 33*z3 - 3*z2*(17-29*z)*u
    + 3*z*(2-9*z)*(4-3*z)*u2 - 3*(2-15*z+35*z2-11*z3)*u3
    - (10-17*z+34*z2-6*z3)*u4 - (5+11*z+4*z2)*u5 - 6*z*u6 )
    * lnz*lnz
+ 1./3. * z * (1-z)/u3/(u-z) * ( 22*z3 - 2*z2*(17-7*z)*u
     + 4*z*(1-z)*(4-z)*u2 - 2*(2-2*z+3*z2)*u3 - (2+z)*u4 - u5 )
    * ( lnz * ( lnu + 4 * log((r+t)/(1.+u)) + log((r-t)/(r+t)) )
      - lnu * lnrpt )
+ 2./3. * z * (1-z)/u/(u-z) * ( 6*z2*(1+u) - 2*u*z + (2-z)*u2 + u3 )
    * lnz * log((r+t)/(1.+u))
+ 1./12. * z * (1-z)/u3/up12 * ( 22*z2 - 2*z*(6-29*z)*u
    + 2*(2-17*z+27*z2)*u2 + 2*(7-13*z+11*z2)*u3 + (3-6*z+4*z2)*u4
    - 2*(2+z)*u5 - u6 )
    * lnu*lnu
+ 2./3. * z * (1-z) * (2+u) * lnz * lnrmt
- 1./3. * z * (1-z) * u * (2+z + u)/(u-z)
    * ( lnz * lnrpt - lnu * log((r-t)/(r+t)) )
+ 1./3. * z * (1-z)/u3 * ( 22*z2 - 2*z*(6-7*z)*u + 2*(2-z+2*z2)*u2
    + 2*(1+z)*u3 + u4 )
    * ( (lnu + 2 * log((r+t)/(1.+u))) * lnomtr
      + log((u-z)/(1.+u))*log((u-z)/(1.+u)) - lnrmt*lnrmt )
// rest of q^1 terms:
- 2./9. * z * (1-z)*(1-z) * (4+7*z+4*z2) * (1+u)/(1+z)/(u-z)
    * log(1.-z)
- 4./9. * z * (1-z) * (10+12*z+3*z2) * (1+u)/(1+z)/(u-z)
    * lnhfopz
- 1./9. * z * (1-z)/(1+z)/u3/(1-u)/up12/(u-z) * ( 158*z3*(1+z)
   - 2*z2*(1+z)*(168+11*z)*u + z*(1+z)*(186-120*z-377*z2)*u2
   - (44-27*z-672*z2-524*z3+73*z4)*u3
   - (16+342*z+123*z2-388*z3-193*z4)*u4
   + (1+z)*(96-132*z-207*z2+101*z3)*u5
   + (106+168*z-51*z2-85*z3+20*z4)*u6 + (2-3*z-33*z2-24*z3)*u7 )
    * lnz
- 1./18. * z * (1-z)/u2/(1-u)/up12/(u-z) * ( 66*z2*(1-2*z)
    - 6*z*(12-48*z+25*z2)*u + (18-132*z+159*z2+88*z3)*u2
    - (6-30*z+249*z2-139*z3)*u3 - (69-180*z+216*z2-47*z3)*u4
    - (69-105*z+48*z2-8*z3)*u5 - 3*(6-11*z)*u6 )
    * lnu
//
+ 2./9. * z * (1-z)/(1+z)/u3/(1+u)/(u-z) * ( 158*z3*(1+z)
    - 2*z2*(1+z)*(135-121*z)*u + 3*z*(1+z)*(38-144*z+27*z2)*u2
    - (22-174*z-30*z2+171*z3+7*z4)*u3
    - (50-15*z-45*z2+20*z3+4*z4)*u4 - (10-9*z-21*z2-4*z3)*u5 )
    * log((r+t)/(1.+u))
- 2./9. * z * (1-z)*(1-z) * (1+u) * (4+7*z+4*z2)/(1+z)/(u-z)
    * lnrmt
- 1./9. * z * (1-z)/u2/up12/(u-z) * ( 66*z2*(1+2*z)
    - 6*z*(12+5*z-47*z2)*u + (3-2*z)*(6-24*z-97*z2)*u2
    + (24+24*z-168*z2+55*z3)*u3 - (3-30*z+30*z2-8*z3)*u4
    - 3*(2-z)*u5 )
    * (lnrmt + lnomtr)
//
+ 2./9. * t/r * z * (1-z)/u2/(u-z) * ( 66*z2 - 2*z*(40-21*z)*u
 + 12*(1-z)*(2-z)*u2 + 6*(1-z)*u3 - 3*u4 )
    * ( lnrmt + lnomtr - lnz + 1./2. * lnu )
+ 2./3. * t/r * z * (1-z)/u/(1-u)/up12/(u-z) * ( 2*z2 - 2*z*(5-z)*u
    + 2*(1-z)*(3+z)*u2 + (1-2*z-2*z2)*u3 + (7+4*z)*u4 - u5 - u6 )
    * lnu
// q^0:
+ 1./54. * z * (1-z)/(1+z)/u2/up12/(u-z) * ( 6*z2*(1+z)*(191+184*z)
    - 18*z*(1+z)*(64-65*z-109*z2)*u 
    - (74+1926*z+2937*z2-153*z3-598*z4)*u2
    - (546+1116*z+2778*z2+697*z3+409*z4)*u3
    - (519-315*z+594*z2-316*z3+176*z4)*u4
    - (236-117*z+111*z2-176*z3)*u5 )
- 1./27. * t/r * z * (1-z)/u2/(1+u)/(u-z) * ( 552*z2 - 2*z*(269-501*z)*u
   - 2*(27+476*z-210*z2)*u2 - 3*(24+189*z+7*z2)*u3 
   + 63*(1-z)*u4 - 18*u5 ) 
;
 return temp;
}

// u -> 1 limit of qq11_full:

double qq11_u1(double z){
 double Z2,t,z2,z3,z4, ln2,lnz,ln1mt,ln1pt,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;  z3 = z2*z;  z4 = z2*z2;
 ln2 = log(2.);  lnz = log(z);  ln1mt = log(1.-t);  
 ln1pt = log(1.+t); lnhf1pt = log((1.+t)/2.);
 temp =
// q^2:
 1./3. * z * ( 
  + 12. * z * (1-z) * J21_u1(z)
  + (1-z) * (40*z2-12*z+7)
    * ( 2. * li2((1.-t)/2.) + lnhf1pt * (lnhf1pt + 4*ln1mt) )
  + 8. * z * (li2(z)-li2(-z))
  - 8. * (10*z3-15*z2+5*z-1) * li2(-t)
  - 8. * (10*z3-13*z2+5*z-1) * li2(t)
  - 1./4. * (240*z3-320*z2+69*z-21) * lnz*lnz
  + 4. * (10*z3-13*z2+3*z-1) * lnz * ln1mt
  + 2. * (60*z3-78*z2+25*z-9) * lnz * ln1pt
  - 2. * z/(1+t)/(1+t) * (7+7*z+8*t) * lnz * log((1.+z)/2.)
  + 2. * (1-z) * (11-28*z+80*z2) * ln2 * lnz 
  - 4. * (1-z)*(10*z2-z+1) * Z2 )
// q^1:
+ 1./18. * z/(1+z) * ( 
+ (201*z4-480*z3*t+330*z3+8*z2*t-132*z2+380*t*z-174*z-108*t+71)
    * lnz
+ (1-t) * (-87*z2*t-71*z2-623*t*z+137*z+119*t-97+671*z3*t-289*z3)
    * ln1mt
+ 4. * (235*z4-196*z3-249*z2+150*z-40) * ln1pt
- 16. * (3*z2+12*z+10) * log(1.+z)
- 4. * (235*z4-204*z3-267*z2+108*z-72) * ln2 )
// q^0:
- 1./216. * z * (1-t)/(1+z) * ( 3079*z3*t-4733*z3+597*z2*t+1265*z2
			      -3991*t*z+4813*z+1051*t+1375 ) ;
  return temp;
}

// u -> z limit of qq11_full:

double qq11_uz(double y, double z){
 double u,r,t,Z2,z2,z3,z4, ln2,lnz,lnhf1pz, lnde,  temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;  z2 = z*z;   z3 = z2*z;  z4 = z2*z2;
 ln2 = log(2.);  lnz = log(z);   lnhf1pz = log((1.+z)/2.);
 lnde = log(r-t);
 temp = 1./3. * (1-z) * (
// q^2:
+ 14 * (2+2*z+z2) * ( li2((1.+z)/2.) - lnhf1pz * (lnde + ln2) )
+ 2 * z * (2+z) * lnz * lnde
+ 2 * (9*z2+22*z+28) * (li2(-z) + ln2 * lnz)
- 2 * z * (2+z) * ( li2(z) + lnz*lnz )
- 1./(1-z) * (29*z3+37*z2+12*z-70) * lnz * lnhf1pz
- z * (2+3*z) * Z2
// q^{1,0}:
+ 1./3./(1+z) * (8*z4+21*z3+45*z2-26*z-54) * (lnde + log(1.-z))
- 1./6./(1+z) * (40*z4+135*z3+411*z2+274*z-54) * lnz
- 8./3. * (z3+8*z2+z-12) * lnhf1pz
- 2./9. * (1-z)/(1+z) * (44*z3+89*z2+143*z+21) ) ;
  return temp;
}

// u -> 1 - z limit of qq11_full (note that it is y-independent):

double qq11_u1mz(double y, double z){
  double Z2,z2,z3, lnz, temp;
  Z2 = ZETA2;  z2 = z*z;  z3 = z2*z;    lnz = log(z);
  temp = 4./9. * ( 
+ 3 * z2 * ( li2(z) - li2(-z) + 3./2. * lnz*lnz
   - lnz * log(1.-z) - lnz * log(1.+z) - 3./2. * Z2 )
+ z/(1+z) * (
- (1-z) * (4+7*z+4*z2) * log(1.-z) + (2-9*z2-6*z3) * lnz
- (10+12*z+3*z2) * log((1.+z)/2.) - 1./3. * (1-z) * (25+7*z+22*z2) ) ) ;
  return temp;
}

// z -> 1 limit of qq11_full:

double qq11_z1(double y, double z){
  double a,a2,a3,a4,a5,y2,y3,y4,lna,lny,ln1my;
  a = 1.-z;  
  a2 = a*a;   a3 = a2*a;   a4 = a2*a2;   a5 = a3*a2;  
  y2 = y*y;   y3 = y2*y;   y4 = y2*y2; 
  lna = log(a);  lny = log(y);  ln1my = log(1.-y);
return a*(-2/3.*lna-1/3.*(-1+y)/y*ln1my+1/3.-1/3.*lny)
  +a2*(-1/3.*(-4+y)*lna-1/6.*(-1+y)*(y-3)/y*ln1my-1/6.*(-4+y)*lny+1/12.*y-1)
  +a3*(-1/6.*(y2+7-3*y)*lna-1/12.*(-1+y)*(-2*y+5+y2)/y*ln1my
    -1/12.*(y2+7-3*y)*lny+5/12.*y-19/72.*y2+1/2.)
  +a4*(1/36.*(8-24*y2+19*y3+6*y)*lna+1/72.*(-1+y)*(19*y3-5*y2+y+9)/y*ln1my
    +1/72.*(8-24*y2+19*y3+6*y)*lny-73/72.*y-17/216.+13/9.*y2-715/864.*y3)
  +a5*(-1/720.*(10*y2+520*y-345-570*y3+311*y4)*lna
    -1/1440.*(-1+y)*(-1+y)*(311*y3+52*y2-197*y+74)/y*ln1my
    -1/1440.*(10*y2+520*y-345-570*y3+311*y4)*lny
    -851/5760.-6913/8640.*y2+25057/86400.*y4+47/720.*y3+2719/4320.*y);
}

// Full, patched answer:

double qq11_full(double y, double z){
 // z = 1 patch:
 if (std::fabs(1.-z) < 0.05){ return qq11_z1(y,z); }
  // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001/sqrt(1.-z)){ return qq11_u1(z); }
 // u = z patch:
 if (y < 0.00000001/(1.-z)){ return qq11_uz(y,z); }
 // u = 1-z patch (different from u = z):
 if (1.-y < 0.00000001/(1.-z)){ return qq11_u1mz(1.-y,z); }
 // (z = (2*u/(1+u))^2 patch (average values on either side of strip):
  double eps = 1.0e-07;
  if (std::fabs(y-y_sing(z)) < eps){ 
    return 0.5 * ( qq11_full(y_sing(z)-2.*eps,z) 
		 + qq11_full(y_sing(z)+2.*eps,z) );
  }
  if (std::fabs(y-(1.-y_sing(z))) < eps){ 
    return 0.5 * ( qq11_full(1.-y_sing(z)-2.*eps,z) 
		 + qq11_full(1.-y_sing(z)+2.*eps,z) );
  }
 return qq11(y,z);
}

// The mu-dependent qq11 hard term:

double qq11_mu(double y, double z, double muFQ){
  double u,t,r,Z2,z2,z3,z4, u2,u3,u4,u5, lnmu, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2; z2 = z*z;  z3 = z2*z;  z4 = z2*z2;
 u2 = u*u;  u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;
 lnmu = 2.*log(muFQ);
//
 temp =  1./3. * z * (1-z) * lnmu * (
- 2/u3 * ( u4 + 2*(1+z)*u3 + 2*(2-z+2*z2)*u2 + 2*z*(7*z-6)*u + 22*z2 )
        * log((r+t)/t/(1+u))
+ 1/(u-z) * ( 4*z - (2-3*z)*u - u2 ) * log(z)
- u/(u-z) * ( 2+z + u ) * log(u)
//
- 1./3./(1+z)/u2/(1+u)/(1+u)/(u-z) * ( (3*z2+8*z3-3*z-2)*u5
   - (-18*z2+21+45*z-46*z3+8*z4)*u4 
   - (-162*z2+48+55*z4-137*z3+66*z)*u3
   - (194*z4-333*z2-60*z+26-57*z3)*u2 
   - 6*z*(1+z)*(47*z2-5*z-12)*u - 66*z2*(1+z)*(1+2*z) )
+ 2./3. * t/r/u2/(u-z) * ( 3*u4 + 6*(-1+z)*u3 - 12*(-1+z)*(-2+z)*u2
			 - 2*z*(-40+21*z)*u - 66*z2 ) );
 return temp;
}

//====================================================================
// "11" "Boost" soft limits:

// mu-dep. part:

double qq11_boost_mu(double z, double muFQ){
 double Z2,z2,z3, lnz, lnmu, temp;
 Z2 = ZETA2;  lnz = log(z);
 lnmu = 2.*log(muFQ);
 z2 = z*z;  z3 = z2*z;
 temp = 1./36./z * (
 lnmu*lnmu * ( (1-z) * (4+7*z+4*z2) + 6 * z * (1+z) * lnz )
+ lnmu * ( - 4 * (1-z) * (4+7*z+4*z2) * log(1.-z)
  + 12 * z * (1+z) * ( li2(z) - li2(-z) - lnz * log((1.+z)*(1.-z)) 
                     + 3./2. * lnz*lnz - 3./2. * Z2 ) 
  + 4 * (2-9*z2-6*z3) * lnz - 4 * (10+12*z+3*z2) * log((1.+z)/2.)
  - 4./3. * (1-z) * (25+7*z+22*z2) ) );
 return temp;
}

// now includes mu-dependent terms:

double qq11_boost_tot(double z, double muFQ){ 
 double Z2,Z3, z2,z3, ln2,lnz,lnz2,lnz3,ln1mz, lnhf1pz,lnhf1pz2, temp;
 Z2 = ZETA2;  Z3 = ZETA3;
 ln2 = log(2.);  lnz = log(z);  lnz2 = lnz*lnz;  lnz3 = lnz2*lnz;
 ln1mz = log(1.-z);
 lnhf1pz = log((1.+z)/2.);  lnhf1pz2 = lnhf1pz*lnhf1pz;
 z2 = z*z;  z3 = z2*z;
 temp = 1./9./z * (
// q^3
- 3 * z * (1+z) * ( 
   li3((1.-z)/2.) + li3(2.*z/(1.+z)) + li3((1.-z)/(1.+z))
 + ( 2 * log((1.-z)/z) + lnhf1pz ) * (li2(z) - li2(-z))
 - 13./12. * lnz3 + 3./2. * lnz2 * log((1.-z)*(1.+z))
 - 2. * lnz * ln1mz * log(1.+z) + 1./2. * log((1.-z)/2.) * lnhf1pz2
 - 1./3. * lnhf1pz2*lnhf1pz - ( 3. * log((1.-z)/z) + 1./2. * lnhf1pz) * Z2 
 - Z3 )
// q^2
- (2+3*z+6*z2+2*z3) * li2(z)
- (16+12*z+3*z2) * (li2(-z) + ln2 * lnz)
- 1./2. * (1+z) * (16+5*z+4*z2) * (li2((1.+z)/2.) - ln2 * lnhf1pz)
+ 1./8. * (8-9*z-81*z2-48*z3) * lnz2
- (6+3*z-12*z2-10*z3) * lnz * ln1mz
- (26+30*z+3*z2) * lnz * lnhf1pz
+ 1./2. * (24+27*z+3*z2-4*z3) * ln1mz * lnhf1pz
+ 1./4. * (36+45*z+15*z2+4*z3) * lnhf1pz2 
+ 1./2. * (4+7*z+4*z2) * ( 2*(1-z)*ln1mz*ln1mz + 3*z * Z2 )
// q^1
+ 2./3. * (1-z) * (25+7*z+22*z2) * ln1mz
- 1./6. * (50-129*z+30*z2-76*z3) * lnz
+ 1./3. * (20-9*z+51*z2) * lnhf1pz
// q^0
+ 1./36. * (1-z) * (268+97*z+136*z2) );
 return temp + qq11_boost_mu(z,muFQ);
}

//====================================================================
// "11" & "22" "Real" soft limits:

// qq11 DD1(n,y) terms, including conversion factor (for y = 0):
// auxiliary function (now includes mu-dep. terms):

double qq11_mu_aux(double y, double z, double muFQ){
 double z2,z3, lnmu;
 z2 = z*z;    z3 = z2*z;  
 lnmu = 2.*log(muFQ);
 return - 2./9. * z/(1+z) * ( (1-z) * (4+7*z+4*z2) + 6 * z * (1+z) * log(z) ) 
   * lnmu * DD1(0.,y) ;
}

double qq11_soft_aux(double y, double z, double muFQ){
 double Z2,z2, ln2,lnz,ln1mz;
 Z2 = ZETA2;
 ln2 = log(2.);  lnz = log(z);  ln1mz = log(1.-z);   z2 = z*z;
 return 4./9. * z/(1+z) * (
+ 1./2. * ( 6. * z * (1+z) * lnz + (1-z) * (4+7*z+4*z2) ) 
      * DD1(1.,y)
+ ( 3 * z * (1+z) * ( li2(-z) - li2(z) - 3./2. * lnz*lnz
                    + lnz * log((1.-z)*(1.+z)) + 3./2. * Z2 )
  - (2-9*z2-6*z2*z) * lnz + (1-z) * (4+7*z+4*z2) * ln1mz
  + (10+12*z+3*z2) * log((1.+z)/2.) + 1./3. * (1-z) * (25+7*z+22*z2) ) 
      * DD1(0.,y) )
  + qq11_mu_aux(y,z,muFQ);
}

double qq11_real_soft(double y, double z, double muFQ){
 return real_conv(0.,z) * qq11_soft_aux(y,z,muFQ);
}

// A remainder term, from "soft", treated as "hard" (no subtraction):

double qq11_soft_to_hard(double y, double z, double muFQ){
 return ( real_conv(y,z) - real_conv(0.,z) ) * qq11_soft_aux(y,z,muFQ);
}

//====================================================================
// The full qq11 real "hard" term, including conversion factor,
// and adding in all other terms multiplying "qq_11_lumi_ys":
// includes conversion factor:

double qq11_real_hard(double y, double z, double muFQ){
 return real_conv(y,z) * ( qq11_full(y,z) + qq11_mu(y,z,muFQ) )
     + qq11_soft_to_hard(y,z,muFQ) + qq11_real_soft(y,z,muFQ) ; 
}

//====================================================================
// "12" terms (there are no "soft" terms here).
// Main "hard" function, making use of symmetry:

double qq12_full_1(double y, double z){
  double u,t,r,Z2,z2, u2,u3,u4,u5, up12, 
    ln2,lnz, lnr,lnu,lnopu,lnhfopu,lnumo, lnrpt,lnrmt, lnhfopz,lnoptr,lnomtr, 
    lnrmo, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;  u2 = u*u;  u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;
 up12 = (1.+u)*(1.+u); 
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnhfopu = log((1.+u)/2.);  lnumo = log(std::fabs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);  lnhfopz = log((1.+z)/2.); 
 lnoptr = log(1.+t*r);  lnomtr = log(1.-t*r);  lnrmo = log(std::fabs(r-1.));  
 temp = 
// q^2 terms (mostly):
// J's:
- 4./3. * (1-z) * ( 4-2*z+z2 - 2*z*(3-z)*u + 3*z2*u2 )
    * ( J2(y,z) 
      - 1./2./(1+u) * ( lnomtr * lnhfopu + lnoptr*lnopu ) )
- 4./3. * (1-z) * (4 + 2*z + z2)
   * J3(y,z) 
//
- 8./3. * z2 * (1-z)/u/(z*up12-4*u2)
   * ( z2 - 2*z*(1-z)*u - z*(10-z)*u2 + 2*(5-z)*u3 )
    * ( J21(y,z) 
      - (1+z) * (1-t*r)/(1+u)/r/(r+t)/(1-z)/(1-z)
        * lnz * lnhfopz
// borrowed from q^1 terms:
      + 1./z/(1+u) * ( log((r+t)/(1.+u)) + 2.*u/(1-u) * lnz 
                    - 1./2. * (1+u)/(1-u) * (1+t/r) * lnu ) )
// polylog's
- 1./3. * z * (1-z)/u * ( 2*z - (2-z)*u )
   * ( li2(z) + lnz * log(1.-z) )
- 4./3. * (1-z) * (1+u)/(1-u)/u * ( z2 + 2*(1-z)*u )
   * li2(-u)
+ 2./3. * (1-z)/u/(1-u)/(1+u) * ( 8*z2 + (4-18*z+z2)*u
   + (28-4*z+z2)*u2 - z*(10-z)*u3 + 5*z2*u4 )
   * ( li2((1.-u)/(1.+u)) - ln2 * lnopu )
//
- 2./3. * (1-z)/(1-u)/(1+u) * ( 4-2*z+z2 - (4+8*z-z2)*u
    + z*(2+3*z)*u2 - z2*u3 )
   * ( li2(t*r) + 1./2. * lnz * lnomtr )
+ 2./3. * (1-z)/(1-u)/(1+u) * ( 3*z2-6*z+12 + (20-12*z+7*z2)*u
   - 3*z*(2-z)*u2 - z2*u3 )
   * li2(-t*r)
- 2./3. * (1-z)/(1-u)/(1+u) * ( 4-2*z+z2 + (2-z)*(6-z)*u - z*(14-3*z)*u2
   + 7*z2*u3 )
   * li2((1.-t*r)/2.)
+ 4./3. * z * (1-z) * (1+u)/u/(1-u) * ( z - (2+z)*u + z*u2 )
   * ( li2(-r*(t-1.)/(r+1.)) - li2((r-1.)/(r+t))
     + lnrpt * log(1.+t) )
- 2./3. * (1-z)*(1+u)/u/(1-u) * ( 4*z2 + (4-6*z-z2)*u )
   * li2(-(r-1)*(r+1)/(1+t*r))
// pure ln's
+ 1./6. * z * (1-z)/u/up12 * ( 17*z - 5*(6-7*z)*u - 2*(15-8*z)*u2 )
   * lnz*lnz
- 5./6. * z2 * (1-z) * (1-u) * (1+u)/u 
   * lnz * lnu
+ 1./3. * z * (1-z)/u * ( 10*z - (6-z)*u )
   * lnz * lnopu
- 1./3. * z * (1-z) * (1+u)/u/(1-u) * ( 2*z - (2+z)*u )
   * log(1.+t) * lnu
- 1./6. * (1-z)/u/up12 * ( 8*z2 - (18-22*z-9*z2)*u - (18-22*z+z2)*u2 )
   * lnu*lnu
+ 2./3. * (1-z)/u/(1-u)/(1+u) * ( 5*z2 + (2-14*z+z2)*u + 2*(7-z+z2)*u2 )
   * ln2 * lnu
- 2./3. * (1-z)/u * ( 3*z2 + 2*(1-2*z)*u ) 
   * lnumo * lnhfopu
+ 1./3. * (1-z)/u * ( 7*z2 - (2+6*z-z2)*u )
   * ln2*ln2
//
- 4./3. * (1-z)/u/(1-u)/(1+u) * ( 5*z2 + 2*(3-6*z+z2)*u 
    + (10-4*z+z2)*u2 )
   * lnr * lnopu
+ 4./3. * (1-z)/u/(1-u)/(1+u) * ( 5*z2 - 2*(3+2*z+2*z2)*u 
    + (6+8*z-3*z2)*u2 )
   * lnr * log(r+1.)
+ 2. * (1-z) * (1-u)/u/(1+u) * ( 4*z2 - (2+2*z-3*z2)*u )
   * lnr * lnrmo
//
- 1./3. * (1-z)/u/(1-u)/(1+u) * ( 2*z2 + 2*(4-2*z+3*z2)*u 
    + (2-3*z)*(2-z)*u2 + 2*(2+2*z-z2)*u3 - 5*z2*u4 )
    * lnoptr*lnoptr
+ 1./3. * (1-z)/(1-u)/(1+u) * ( 3*(4-2*z+z2) + (2-z)*(10-3*z)*u
    - z*(18-5*z)*u2 + 9*z2*u3 )
   * lnz * lnoptr
+ 1./6. * (1-z)/u/(1-u)/(1+u) * ( 4*z2 - (4+2*z-13*z2)*u 
    + 2*(28+5*z2)*u2 + (12-46*z-13*z2)*u3 + 10*z2*u4 )
   * lnu * lnoptr
+ 4./3. * (1-z)/(1-u) * ( 4-2*z+z2 - 4*z*u +  2*z2*u2 )
   * ln2 * lnoptr
+ 1./3. * (1-z)/u * ( 8*z2 - (4+2*z-5*z2)*u )
   * lnrpt * lnoptr
+ 1./3.*(1-z)/(1-u)/(1+u) * ( 4-2*z+z2 - (4-z2)*u + z*(10-z)*u2
   - 5*z2*u3 )
   * lnu * lnomtr
//
+ 1./3. * (1-z)/u * ( 3*z2 + (10-6*z-3*z2)*u )
   * Z2
// q^1
- 4./3. * z * (1-z) * (u + 4*z*(1+u)) * (2 + 4*u + 3*u2)/up12
    * J27(y,z)
+ 2./3. * z/u/up12 * ( 3*z - 2*z2 + u - (12-17*z+4*z2)*u*(1+u) )
    * lnz 
- 4./3. * t/r * z * (1-z)/u/(1-u)/up12
    * ( z - 2*(1-z)*u - (6-z)*u2 - 2*u3 )
    * lnu
- 2./3. * z * (1-z)/u/up12 * ( 2*z - (7+4*z)*u - 2*(4+3*z)*u2 )
    * log((1.+u)/r)
- 1./3. * z/u/(1-u)/(1+u) * ( 2*z*(2-3*z) + (1+3*z-2*z2)*u2 - 17*(1-z)*u )
    * lnu
- 2./3. * z/u/(1-u)/up12 * ( 2*z2 + (7-z-4*z2)*u - (3-3*z+4*z2)*u2
  - 4*(3-z-z2)*u3 + 2*u4*z2 + 2*u5*z )
    * log((r+t)/t)
// q^0
+ 8./3. * z * (1+u+u2)*(1-z)*(1-t*r)*(r-t)/r/up12
;
 return temp;
}

// u -> 1 limit of qq12_full:

double qq12_u1(double z){
 double Z2,t,z2, ln2,lnt,lnz,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;
 ln2 = log(2.);  lnt = log(t);  lnz = log(z);  lnhf1pt = log((1.+t)/2.);
 temp = 1./3. * ( 
// q^2:
- 2 * (1-z) * (3*z2+2-4*z) * ( 8 * J2_u1(z) + lnz * log(1.-z) )
- 8 * (1-z) * (z2+2*z+4) * J3_u1(z)
- 8 * (1-z) * z2 * (2*z-5) * J21_u1(z)
- 2 * z * (1-z) * (-2+3*z) * li2(z)
- 8 * (1-z)*(1-z) * li2(t) + 4 * (1-z) * (z2+6) * li2(-t)
- 8 * (1-z)*(1-z) * (1+4*z) * li2((1.-t)/2.)
//
+ z * (1-z) * (17*z-15) * lnz*lnz
+ 2 * (1-z) * (8+4*z-10*z2) * lnz * log(1.+t)
+ 4 * (1+z) * (2*z-5) * z2/(t+1)/(t+1) * lnz * log(1.+z)
- 2 * z/(t+1)/(t+1) * (15*z2*z+22*t*z2-12*z2-34*t*z-21*z+12*t+6) * ln2 * lnz
+ 4 * (1-z) * (z2+z-3) * lnhf1pt*lnhf1pt
+ 4 * (1-z) * (2*z2-3*z+5) * Z2
// q^{1,0}
- 18 * z * (1-z) * (1+8*z) * J27_u1(z)
+ 1./2. * t * (20*t*z2-24*z2+3*t*z+48*z+t-32) * lnz
+ t * (16*t*z2+24*z2-21*t*z-48*z+5*t+32) * lnhf1pt
+ 4 * t * (1-t) * (7*z2-3*t*z-10*z+3*t+8) );
  return temp;
}

// u -> z limit of qq12_full:

double qq12_uz(double y, double z){
 double u, Z2,z2,z3, ln2,lnz,ln1pz, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));
 Z2 = ZETA2;  z2 = z*z;   z3 = z2*z;
 ln2 = log(2.);  lnz = log(z); 
 ln1pz = log(1.+z);
 temp = 1./3. * z/(1+z) * ( 
- 4 * (1-z) * (1+z2) * ( J2rem(z) + li2(z) - 1./2. * lnz*lnz
      + (2. * lnz - ln1pz + ln2) * log(1.+2.*z) + 3./4. * ln1pz*ln1pz )
- 4 * (-2*z2+2*z3+z-3) * li2(-z)
- 2 * (3+z-3*z2+3*z3) * li2((1.+z)/2.)
+ 2 * (1-z) * (1+5*z2) * ln2 * lnz
- 2 * (-3*z2+5*z-9+3*z3) * lnz * ln1pz
+ 4 * (1+z) * (2-2*z+z2) * ( ln2 * log(1.-z) - ln1pz * log((1.-z)/2.) )
- (-7*z2+3*z+5+7*z3) * ln2*ln2 - 2 * (-2+z)*(z2+z+4) * Z2
- 4 * z * (1+z) * log((1.+z)/2./z)
 );
  return temp;
}

// z -> 1 limit of qq12_full (zero is good to O(delta^2)):

double qq12_z1(double y, double z){
  double a,a2,a3,a4,a5,a6,y2,y3,y4;
  a = 1.-z;  
  a2 = a*a;   a3 = a2*a;   a4 = a2*a2;   a5 = a3*a2;    a6 = a3*a3; 
  y2 = y*y;   y3 = y2*y;   y4 = y2*y2; 
return 7/12.*a2+a3*(47/108.+47/108.*y2-47/108.*y)
  +a4*(-371/1728.+731/864.*y-731/864.*y2)
  +a5*(5917/43200.*y4-5917/21600.*y3-613/800.*y2+39019/43200.*y-12073/43200.)
  +a6*(8249/43200.*y4-8249/21600.*y3-121/300.*y2+25673/43200.*y-383/1800.);
}

// Full, patched answer, using symmetry:

double qq12_full(double y, double z){
 // z = 1 patch:
 if (std::fabs(1.-z) < 0.02){ return qq12_z1(y,z); }
  // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001/sqrt(1.-z)){ return qq12_u1(z); }
 // u = z patch:
 if (y < 0.00000001/(1.-z)){ return qq12_uz(y,z); }
 // u = 1-z patch:
 if (1.-y < 0.00000001/(1.-z)){ return qq12_uz(1.-y,z); }
 // (z = (2*u/(1+u))^2 patch (average values on either side of strip):
  double eps = 1.0e-07;
  if (std::fabs(y-y_sing(z)) < eps){ 
    return 0.5 * ( qq12_full(y_sing(z)-2.*eps,z) 
		 + qq12_full(y_sing(z)+2.*eps,z) );
  }
  if (std::fabs(y-(1.-y_sing(z))) < eps){ 
    return 0.5 * ( qq12_full(1.-y_sing(z)-2.*eps,z) 
		 + qq12_full(1.-y_sing(z)+2.*eps,z) );
  }
 return qq12_full_1(y,z) + qq12_full_1(1.-y,z);
}

// The full qq12 real "hard" term, including conversion factor:

double qq12_real_hard(double y, double z){
 return real_conv(y,z) * qq12_full(y,z) ; 
}

//====================================================================
// quark-quark exchange "CE" terms. 
// There are "soft" terms here, and no y -> 1-y symmetry.
// Main "hard" function:

double qq_CE(double y, double z){
  double u,t,r,Z2,z2, u2,u3,u4, up12, 
    ln2,lnz, lnr,lnu,lnopu,lnhfopu,lnumo, lnrpt,lnrmt, lnhfopz, 
    lnrmo, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;  u2 = u*u;  u3 = u2*u;  u4 = u2*u2;
 up12 = (1.+u)*(1.+u); 
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnhfopu = log((1.+u)/2.);  lnumo = log(std::fabs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);  lnhfopz = log((1.+z)/2.); 
 lnrmo = log(std::fabs(r-1.));  
 temp = 1./9. * z * (1-z) * (
// q^2:
+ 2./u * ( z - (1-z)*u + u2 ) * ( 2. * (1+u)/u * J2(1.-y,z)
   - li2((1.-u)/2.) - 2 * li2((1.-u)/(1.+u)) 
   + li2(t/r) - li2((r-t)/2./r) 
   + 1./2. * log(u-z) * log(z/(1.+u)/(1.+u)) + 2. * log(2.*u) * lnopu
   - 1./2. * lnopu*lnopu - 1./2. * ln2 * lnu
   + lnumo * log((1.+u)/2./u/u)
   + 1./2. * lnrpt*lnrpt - lnrpt * lnu + lnrmt * (3./2. * lnu + ln2)
   - ln2*ln2 - 1./2. * Z2 )
+ 4. * (1+u2)/(u-z) * ( J3(1.-y,z) 
   + 1./4. * u/(1+u) * (lnz*lnz + 2. * lnu * lnopu) )
//
- 8. * z * (1+z2) * (1+u)/(1+z)/(1+z)/(u-z)
    * ( li2(-z) - 1./4. * lnz*lnz + lnz * log(1.+z) + 1./2. * Z2 )
+ 4./u/(1+u)/(u-z) * ( 2*z2 - 2*z*(1-z)*u + (1-2*z)*u2 + u4 )
    * ( li2(-u) - 5./8. * lnz*lnz - 1./4. * lnu*lnu 
      - 1./4. * lnz * lnu - 1./2. * lnz * lnopu )
+ 2./u/(1+u)/(u-z) * ( 5*z2 - 6*z*(1-z)*u + (3-6*z+z2)*u2 + u4 )
    * ( li2(-t/r) + 3./4. * lnz*lnz + lnz * lnopu 
      + 1./8. * lnu*lnu + Z2 )
// pure ln's:
- 2. * z/u * log((r+t)*(r+t)*u/z) * lnu
// q^1:
+ 4./u * ( 1 + 4*z*(1+u) ) * J27(1.-y,z)
+  2./u * ( 2*z -(1-2*z)*u ) * log((r+t)/(1.+u))
+ 1./u/(u-z) * ( 2*z2 - z*(1+6*z)*u + 3*(1+2*z)*u2 ) * lnz
- ( 2 + 4 * t/r/(u-z) * ( 2*z - (1-z)*u - u2 ) ) * lnu  
// q^0:
+ 4./(1+z)/(u-z) * ( z*(3+2*z+3*z2) - (1+3*z2)*u )
- 4. * t/r/(u-z) * ( 3*z - (1-z)*u - u2 )
 );
 return temp;
}

// u -> 1 limit of qq_CE_full:

double qq_CE_u1(double z){
 double Z2,t,z2, ln2,lnz,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;
 ln2 = log(2.);  lnz = log(z);  lnhf1pt = log((1.+t)/2.);
 temp = 1./9. * z * ( 
// q^2:
- 2 * z * (1-z) * ( - 8 * J2_u1(z) + 2 * li2((1.-t)/2.) - 2 * li2(t)
                  - lnhf1pt*lnhf1pt - lnz * log(1.-z) )
+ 8 * J3_u1(z)
- 16 * z * (1+z2)/(1+z)/(1+z) * ( li2(-z) + lnz * log(1.+z) )
+ 4 * (1-3*z+3*z2) * li2(-t)
+ 1./2. * (6*z+16*z2*z-5*z2+8*z2*z2+3)/(1+z)/(1+z) * lnz*lnz
+ 2 * (1-2*z)*(1-2*z) * log(2.) * lnz
+ 2 * (-7*z-4*z2+z2*z+5*z2*z2+1)/(1+z)/(1+z) * Z2
// q^{1,0}:
+ (1-z) * ( 4 * (1+8*z) * J27_u1(z) + 2 * (4*z-1) * lnhf1pt
          - (4*z2-5*z-3)/(1-z) * lnz
	    + 4 * (-3*t*z2+z2+2*t*z+4*z+t-1)/(1+t)/(1+z) ) );
 return temp;
}

// u -> z limit of qq_CE_full:

double qq_CE_uz(double y, double z){
 double u, Z2,z2,z3, ln2,lnz,ln1pz, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));
 Z2 = ZETA2;  z2 = z*z;   z3 = z2*z;
 ln2 = log(2.);  lnz = log(z);  ln1pz = log(1.+z);
 temp = 2./9. * z * (1-z) * (
+ 4 * z * ( J2rem(z) + li2(z) - 1./2. * li2((1.+z)/2.) 
          - log((1.+z)/2.) * log(1.+2.*z) + 3./4. * ln1pz*ln1pz 
          + 2 * lnz * log(1.+2.*z) - 3./4. * ln2*ln2 - 1./2. * Z2 )
- 4 * z * (5+2*z+z2)/(1+z)/(1+z) * ( li2(-z) + 1./2. * Z2 )
- (1+z+7*z2+3*z3)/(1+z)/(1+z) * lnz * ( lnz - 4 * ln1pz )
- 2 * (1+9*z) * lnz * ln1pz - 2 * (1+3*z) * ln2 * lnz
- (1-z) * (4+z+z2)/(1+z)/(1+z) * lnz
- 2 * (1+z2)/(1+2*z) * log((1.+z)/2/z/z)
+ 4 * z * (1-z)/(1+z) );
 return temp;
} 

// z -> 1 limit:

double qq_CE_z1(double y, double z){ 
  double a,a2,a3,a4,a5,a6,a7,y2,y3,y4,y5,y6;
  a = 1.-z;  a2 = a*a;
  a3 = a2*a;  a4 = a2*a2;  a5 = a3*a2;  a6 = a3*a3;  a7 = a4*a3;  
  y2 = y*y;  y3 = y2*y;  y4 = y2*y2;  y5 = y3*y2;   y6 = y3*y3;
  return -1/9.*y*(-1+y)*a3+1/18.*(-1+y)*(-1+y)*a4
   +a5*(1/72.+25/108.*y3-7/120.*y4-5/18.*y2+11/108.*y)
   +a6*(97/2160.*y4-2/81.*y5+2/27.*y3-125/648.*y2+13/108.*y-7/432.)
   +a7*(-149/4320.-793/22680.*y6-229/1296.*y2-227/1440.*y4+77/648.*y5
     +107/648.*y3+43/360.*y);
}

// Full, patched CE answer:

double qq_CE_full(double y, double z){
  // z = 1 patch:
  if (std::fabs(1.-z) < 0.05){ return qq_CE_z1(y,z); }
  // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001){ return qq_CE_u1(z); }
  // u = z patch:
 if (y < 0.0000001){ return qq_CE_uz(y,z); }
 else { return qq_CE(y,z); }
}

//====================================================================
// qq_CE "Boost" soft limits (total now includes mu-dependence):

double qq_CE_boost_mu(double z, double muFQ){
 double Z2,z2, lnz, lnmu, temp;
 Z2 = ZETA2;  lnz = log(z);
 lnmu = 2.*log(muFQ);  z2 = z*z;
 temp = 1./18. * lnmu * (
  - (1+z2)/(1+z) * ( 4 * (li2(-z) + lnz * log(1.+z)) - lnz*lnz + 2 * Z2 )
  + 2 * (1+z) * lnz + 4 * (1-z) );
 return temp;
}

double qq_CE_boost_tot(double z, double muFQ){ 
 double Z2,Z3, z2,z3, ln2,lnz,lnz2,lnz3,ln1mz,ln1pz,lnhf1pz,lnhf1pz2, temp;
 Z2 = ZETA2;  Z3 = ZETA3;
 ln2 = log(2.);  lnz = log(z);  lnz2 = lnz*lnz;  lnz3 = lnz2*lnz;
 ln1mz = log(1.-z);  ln1pz = log(1.+z); 
 lnhf1pz = log((1.+z)/2.);  lnhf1pz2 = lnhf1pz*lnhf1pz;
 z2 = z*z;  z3 = z2*z;
 temp = 1./9. * (
// q^3
 (1+z2)/(1+z) * ( li3(z) + li3(-1./2./z) - li3(1./2./(1.+z))
   + 3. * ( li3(1.-z) - li3(-z) + li3(-(1.+z)/z) )
   - 5./2. * ( li3((1.+z)/2.) + li3((1.-z)/2.) 
           + li3((1.-z)/(1.+z)) )
   - 1./4. * li3((1.-z)*(1.+z))
   + log((1.+z)/2./z/z) * ( li2(-z) - li2(-(1.+z)/z) )
   + 4 * ln1mz * ( li2(-z) + lnz * ln1pz )
   + 11./12. * lnz*lnz*lnz - 1./2. * lnz*lnz * ( ln1mz + 8 * ln1pz )
   + lnz * ( lnhf1pz * ln1pz - 1./2. * ln2*ln2 )
   + 7./12. * lnhf1pz2*lnhf1pz - 5./4. * ln1mz * lnhf1pz2
   + 9./4. * ln2 * lnhf1pz2 + 2 * ln2*ln2 * ln1pz 
   - Z2 * ( 9./2. * lnz - 2 * ln1mz - 4 * lnhf1pz - ln2 ) 
   - 5./6. * ln2*ln2*ln2 + 11./4. * Z3 )
// q^2
+ (1+z) * ( li2(z) + li2(-z) 
          + lnz * log(z*(1.+z)/(1.-z)) - 1./2. * Z2 )
+ 2 * li2(-z) + 2 * ln2 * lnz + Z2 
// q^{1,0}
- 4 * ( lnhf1pz + (1-z) * ln1mz ) + 1./4. * (5-19*z) * lnz
- 15./4. * (1-z) );
 return temp + qq_CE_boost_mu(z,muFQ);
}

//====================================================================
// qq_CE "Real" soft limits:

// qq_CE DD1(n,y) terms, including conversion factor (for y = 0):
// auxiliary function:

double qq_CE_soft_aux(double y, double z){
 double Z2,z2, lnz;
 Z2 = ZETA2;   z2 = z*z;  lnz = log(z);
 return 8./9. * z*z * (
   (1+z2)/(1+z)/(1+z) * ( li2(-z) + lnz * log(1.+z) 
                     - 1./4. * lnz*lnz + 1./2. * Z2 )
   - 1./2. * lnz -  (1-z)/(1+z) ) * DD1(0.,y) ;
}

double qq_CE_real_soft(double y, double z){
 return real_conv(0.,z) * qq_CE_soft_aux(y,z);
}

// A remainder term, from "soft", treated as "hard" (no subtraction):

double qq_CE_soft_to_hard(double y, double z){
 return ( real_conv(y,z) - real_conv(0.,z) ) * qq_CE_soft_aux(y,z);
}

//====================================================================
// The full qq_CE real "hard" term, including conversion factor,
// and adding in all other "CE" terms multiplying "qq_ident_lumi_ys":

double qq_CE_real_hard(double y, double z){
 return real_conv(y,z) * qq_CE_full(y,z) 
      + qq_CE_soft_to_hard(y,z) + qq_CE_real_soft(y,z) ; 
}

//====================================================================
// quark-quark exchange "CF" terms. 
// There are no "soft" terms here, and there is a y -> 1-y symmetry.
// Main "hard" function, using symmetry:

double qq_CF_full_1(double y, double z){
  double u,t,r,Z2,z2, u2,u3,u4,u5, up12, 
    ln2,lnz, lnr,lnu,lnopu,lnhfopu,lnumo, lnrpt,lnrmt, lnhfopz,lnoptr, 
    lnrmo, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;  u2 = u*u;  u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;
 up12 = (1.+u)*(1.+u); 
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnhfopu = log((1.+u)/2.);  lnumo = log(std::fabs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);  lnhfopz = log((1.+z)/2.); 
 lnoptr = log(1.+t*r);  lnrmo = log(std::fabs(r-1.));  
 temp = 1./9. * z * (1-z) * (
// q^2:
+ 2 * (1+u)/(1-u)/u2 * ( z2 - z*(2+z)*u + (1+z2)*u2 )
 * ( 4 * li2(-u) - 2 * li2((1.-u)/(1.+u)) - 2 * li2((1.-u)/2.)
   + 2 * li2((1.-u)/(1.+t*r))
   - lnu * lnoptr + lnoptr*lnoptr 
   + 2 * ln2 * lnopu - ln2*ln2 + 2 * Z2 )
+ 4/(1-u)/u2 * ( z2 - 2*z*u + (2-2*z+z2)*u2 )
    * ( li2(t/r) + li2((r-t)*r/(1.+u))
      + 1./2. * log(u-z) * lnz - lnrmt * lnr
      + lnrpt * lnopu - 1./4. * lnz*lnz - lnz * lnopu
      + 3./2. * lnu * lnopu 
      - 3./2. * Z2 )
- 4/(1-u)/u2 * ( 2*z2 - 4*z*u + 2*(2-2*z+z2)*u2
   + (2-2*z+z2)*u3 - 2*z*u4 + z2*u5 )
    * ( li2(-t*r) + 1./2. * lnr * lnrpt )
// pure ln's
- 2 * (1+z2) * lnu*lnu
- 4 * (1-z)* (1-z) * (1+u)/(1-u)
    * lnu * lnopu
- 2/u2 * ( 2*z2 - 2*z*(2-z)*u + (1-4*z+z2)*u2 )
    * ( lnopu*lnopu + log((r+t)/r) * lnoptr )
- 1/(1-u) * ( 2*(1+z2) - z*(2+z)*u - 2*z*u2 + z2*u3 )
    * lnu * lnrpt
// q^1:
+ 2/u2 * ( 13*z2 - 2*z*(5-2*z)*u + (3-2*z+2*z2)*u2
   - 2*z*u3 )
    * log((r+t)/t/(1.+u))
// q^0:
+ 1./2./u * ( 30*z + 52*z2 - (7+16*z+3*z2)*u ) 
- t/r * (1+u)/u  * ( 26*z - (5+7*z)*u ) );
 return temp;
}

// u -> 1 limit of qq_CF_full:

double qq_CF_u1(double z){
 double Z2,t,z2, lnz,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;
 lnz = log(z);  lnhf1pt = log((1.+t)/2.);
 temp = 2./9. * z * (1-z) * (
- 4 * (1-z) * (1-3*z) * ( li2(-t) - li2(t) - li2((1.-t)/2.)
  + 1./4. * lnz*lnz - 1./2. * lnz * log(1.-z) + 1./2. * lnhf1pt*lnhf1pt
  + log(2.) * lnz + 3./2. * Z2 )
+ (3-10*z+15*z2) * ( 2 * lnhf1pt - lnz )
- 1./2. * (1-t) * (7+3*t-27*z+33*t*z) );
 return temp;
}

// z -> 1 limit:

double qq_CF_z1(double y, double z){ 
  double a,a2,a3,a4,a5,a6,a7,y2,y3,y4,y5,y6;
  a = 1.-z;  a2 = a*a;
  a3 = a2*a;  a4 = a2*a2;  a5 = a3*a2;  a6 = a3*a3;  a7 = a4*a3;  
  y2 = y*y;  y3 = y2*y;  y4 = y2*y2;  y5 = y3*y2;   y6 = y3*y3;
  return a3*(-2/9.-4/9.*y2+4/9.*y)+a4*(17/108.-4/9.*y+4/9.*y2)
    +a5*(13/324.+31/81.*y3-31/324.*y-31/324.*y2-31/162.*y4)
    +a6*(373/25920.+35/216.*y3-19/864.*y-35/432.*y4-17/288.*y2)
    +a7*(-4147/32400.*y6+251/43200.+13459/64800.*y3-1019/129600.*y
         +4147/10800.*y5-54929/129600.*y4-691/21600.*y2);
}

// Full, patched CF answer, using the symmetry:

double qq_CF_full(double y, double z){
  // z = 1 patch:
  if (std::fabs(1.-z) < 0.05){ return qq_CF_z1(y,z); }
  // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001){ return qq_CF_u1(z); }
 else { return qq_CF_full_1(y,z) + qq_CF_full_1(1.-y,z); }
}

// The full qq_CF real "hard" term, including conversion factor:

double qq_CF_real_hard(double y, double z){
 return real_conv(y,z) * qq_CF_full(y,z) ; 
}

//====================================================================
// quark-(anti)-quark axial vector "AB" terms. 
// There are no "soft" terms here, and there is a y -> 1-y symmetry.
// Main "hard" function, using symmetry:

double qq_AB_ax_full_1(double y, double z){
  double u,t,r,Z2,z2,z3,z4,z5, u2,u3, up12, lnz,lnu, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;   z2 = z*z;  z3 = z2*z;  z4 = z2*z2;  z5 = z3*z2;  
 u2 = u*u;  u3 = u2*u;
 up12 = (1.+u)*(1.+u); 
 lnz = log(z);  lnu = log(u);  
 temp = 4./3. * z * (1-z) * ( 
// q^2
- 8 * z * u /(z*up12-4*u2)/(z*up12-4*u2) * ( z2*(1+z) - z*(1+3*z-2*z2)*u
    + z*(1-z)*(3-z)*u2 - 2*(1-z)*(1-z)*u3 )
    * ( J21(y,z) 
      - (1+z) * (1-t*r)/(1+u)/r/(r+t)/(1-z)/(1-z)
        * lnz * log((1.+z)/2.)
// from q^1 terms:
      + 1./z/(1+u) * ( log((r+t)/(1.+u)) + 2*u/(1-u) * lnz 
    - 1./2. * (1+u)/(1-u) * (1+t/r) * lnu ) )
+ 1./2. * z*u/up12 * log(z/u) * log(z*u)
// q^1
+ 2./(1+u)/(z*up12-4*u2) * ( 2*z2 - z*(3-z)*u + (2-3*z+z2)*u2 )
   * log((r+t)/t/(1.+u))
- 2. * u/(1-z)/(1-z)/up12/(z*up12-4*u2)/(z*up12-4) * (
      z*(9-10*z+7*z2+8*z3-2*z4) + 4*(1+z2+4*z4-2*z5)*u
    - (4-7*z+26*z2+7*z3-8*z4+6*z5)*u2 )
   * log(z/u)
- 4./(1+u)/(1-u)/(z*up12-4*u2)/(z*up12-4) * ( 
   + u/(1-z)/(1-z) * ( 2*z*(3-6*z+7*z2-z3) - z * (3+2*z-15*z2+16*z3-2*z4)*u
                 - (4-17*z+30*z2-23*z3+14*z4-6*z5)*u2 )
   + t * r/(1+u) * ( - z*(1+z)*(4-z) + z*(4+4*z2-7*z)*u 
                   + (8+20*z-13*z2+7*z3)*u2 - (16-4*z+9*z2-4*z3)*u3 ) )
   * lnu
// q^0
- 2. * r * (r-t) * (1-t*r)/(1-z)/up12/(z*up12-4*u2)/(z*up12-4)
   * ( z*(1-z)*(4-z) - 4*z2*(7-z)*u - (24-4*z+23*z2-3*z3)*u2 
     - 4. * t * r * (1+u) * ((2+z2)*(1+u) + 2*z*u) ) );
 return temp;
}

// u -> 1 limit of qq_AB_ax_full:

double qq_AB_ax_u1(double z){
 double Z2,t,z2, lnz, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;   lnz = log(z);
 temp = 1./3. * z * (
+ 8 * z * (1-2*z+2*z2) * ( J21_u1(z) 
               - 1./2. * (1+z)/(1-z)/(1+t)/(1+t) * lnz * log((1.+z)/2.) )
+ z * (1-z) * lnz*lnz + 2 * (1+3*z2)/(1-z) * lnz + 4 * z * log((1.+t)/2.)
- 2 * (1+t) * (t*z-5*z+5*t-3) ) ;
 return temp;
}

// z -> 1 limit of qq_AB_ax_full:

double qq_AB_ax_z1(double y, double z){
  double a,a2,a3,a4,a5,a6,y2,y3,y4;
  a = 1.-z;  
  a2 = a*a;   a3 = a2*a;   a4 = a2*a2;   a5 = a3*a2;   a6 = a3*a3;  
  y2 = y*y;   y3 = y2*y;   y4 = y2*y2; 
return 4/9.*a2 + a3*(-2/3.*y-5/9.+2/3.*y2) + a4*(5/9.*y-5/9.*y2+1/270.)
  +a5*(4/135.*y-3/20.*y2+13/360.+13/54.*y3-13/108.*y4)
  +a6*(-1/180.*y+23/108.*y3-109/1080.*y2+131/5040.-23/216.*y4);
}

// z -> 0 limit of qq_AB_ax_full:

double qq_AB_ax_z0(double y, double z){
  double u, t, r;
  t = sqrt(z);  u = (z+y*(1.-z))/(1.-y*(1.-z));   r = sqrt(u);
 return 2./3. * z/(1+u)/(1+u) * ( 4 * u * (log(z)+3)
    + t/r * ( - (1-u)*(1+8*u+u*u) * log(u) + 2 * (1-3*u-3*u*u+u*u*u) ) ); 
}

// Full, patched AB_ax answer, using the symmetry:

double qq_AB_ax_full(double y, double z){
 // z = 1 patch:
 if (std::fabs(1.-z) < 0.05){ return qq_AB_ax_z1(y,z); }
 // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001/sqrt(1.-z)){ return qq_AB_ax_u1(z); }
 // z = 0 patch:
 if (std::fabs(z) < 0.000000001){ return qq_AB_ax_z0(y,z); }
 // (z = (2*u/(1+u))^2 patch (average values on either side of strip):
  double eps = 1.0e-05;
  if ((std::fabs(y-y_sing(z)) < eps) && (y_sing(z) > 2.*eps)){ 
    return 0.5 * ( qq_AB_ax_full(y_sing(z)-2.*eps,z) 
		 + qq_AB_ax_full(y_sing(z)+2.*eps,z) );
  }
  if ((std::fabs(y-(1.-y_sing(z))) < eps) && (y_sing(z) > 2.*eps)){ 
    return 0.5 * ( qq_AB_ax_full(1.-y_sing(z)-2.*eps,z) 
		 + qq_AB_ax_full(1.-y_sing(z)+2.*eps,z) );
  }
 return qq_AB_ax_full_1(y,z) + qq_AB_ax_full_1(1.-y,z);
}

// The full qq_AB_ax real "hard" term, including conversion factor:

double qq_AB_ax_real_hard(double y, double z){
 return real_conv(y,z) * qq_AB_ax_full(y,z) ; 
}

//====================================================================
// quark-(anti)-quark axial vector "CD" terms. 
// There are no "soft" terms here, and there is a y -> 1-y symmetry.
// Main "hard" function, using symmetry:

double qq_CD_ax_full_1(double y, double z){
  double u,t,r,Z2,z2, u2,u3,u4,u5, up12, 
    ln2,lnz, lnr,lnu,lnopu,lnumo, lnrpt,lnrmt, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;  u2 = u*u;  u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;
 up12 = (1.+u)*(1.+u); 
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnumo = log(std::fabs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);
 temp =  4./3. * z * ( 
// q^2 
+ (1-z) * (
- 1./2./(1+u) * ( 2-z - 2*(1+z)*u + z*u2 ) * (
 + 2 * (1+u) * J2(y,z) - li2((1.-u)/(1.+u)) - li2((1.-t*r)/2.)
 + li2(t*r) - log((1.+u)/2.) * log(1.-t*r) - log(1.+t*r) * lnopu )
//
- (2+z) * ( J3(y,z) 
          + 1./4./(1+u) * ( lnz*lnz - 2 * lnu * lnopu ) )
- 2 * z/u/(z*up12-4*u2)
   * ( z2 - 2*z*(1-z)*u - z*(10-z)*u2 + 2*(5-z)*u3 )
    * ( J21(y,z) 
      - (1+z) * (1-t*r)/(1+u)/r/(r+t)/(1-z)/(1-z)
        * lnz * log((1.+z)/2.)
// from q^1 terms:
    + 1./z/(1+u) * ( log((r+t)/(1.+u)) + 2*u/(1-u) * lnz 
    - 1./2. * (1+u)/(1-u) * (1+t/r) * lnu ) )
// polylogs:
- 1./4./u * ( 2*z - (2-z)*u ) 
   * ( li2(z) + lnz * log(1.-z) )
+ (1+u)/u/(1-u) * (z-2*u) * li2(-u)
- 4/(1+u)/(1-u) * (z-2*u) 
   * ( li2((1.-u)/2.) + li2((1.-u)/(1.+u)) )
+ u * (2-z*u)/(1-u) * li2(t*r)
- 1./2./u/(1-u) * ( 4*z - (14-3*z)*u - 6*u2 + 3*z*u3 )
   * li2(-t*r)
+ (2-z) * u/(1-u) 
   * ( li2((r-t)/2./r) - li2((1.-u)/(1.+t*r)) )
- (1+u)/u/(1-u) * ( z - (2+z)*u + z*u2 )
   * ( li2((r-1.)/(r+t)) + li2((1.-t)/(r+1.))  
     - log(r+1.) * log(1.+t) )
// pure lns:
+ 1./8. * 1/u/up12 * ( 3*z - (2+3*z)*u - 2*(1+4*z)*u2 )
   * lnz*lnz
+ 1./2. * 1/(1+u)/(1-u) * ( 2-3*z + (6-z)*u )
    * lnz * lnu
- 1./4./u * ( 2*z - (2-3*z)*u ) 
    * lnz * lnopu 
+ 1./2. * (1-4*u+u2)/u/(1+u)/(1-u) * ( z - 2*u ) 
    * lnu * lnumo
- 1./4. * z * (1-4*u+u2)/u * lnumo * log((1.+u)/2.)
- 1./2. * z/u * (1-6*u) * ln2 * log((1.+u)/r)
- 1./4. * z * (1-u)/u * lnu * lnumo
- z * log((1.+u)/r)*log((1.+u)/r) 
+ 1./4. * 1/u/(1+u)/(1-u) * ( z - 2*(6-z)*u - (4-5*z)*u2 )
    * lnu * lnopu
- 1./16. * 1/u/up12 * ( z - 2*(4+13*z)*u - (8+31*z)*u2 )
    * lnu*lnu
+ 1./2. * 1/(1+u)/(1-u)  * ( 2+3*z - (6+z)*u ) 
    * ln2 * lnu
//
+ 1./2. * (1+u)/u/(1-u) * ( 2*z - (2+z)*u )
   * ( log(1.+t) * log((r+t)/r/(1.+t*r)) + lnr * log(std::fabs(r-1.)/r) )
- 1./4. * 1/u/(1+u)/(1-u) * ( 3*z - (10-z)*u - (8-3*z)*u2 - 3*(2-z)*u3
   + 2*z*u4 )
   * lnrpt*lnrpt
+ 1./2. * 1/u/(1+u)/(1-u) * ( z - 8*z*u + (12+z)*u2 - 2*(2-z)*u3 )
   * lnu * lnrpt
- (2-z) * u/(1-u) * log(2.*z) * lnrpt
- 1./2. * (2-z) * log((r+t)/r) * log(1.+t*r)
//
+ 1./4. * 1/u/(1+u)/(1-u) * ( z - (2-5*z)*u - (8+z)*u2 + (2-z)*u3 )
   * lnz * log(u-z)
- 1./4. * 1/u/(1+u)/(1-u) * ( 3*z - (6+z)*u + z*u2 - (2-z)*u3 )
   * lnu * lnrmt
+ 1./4./u * ( z - (2+5*z)*u ) * ln2*ln2 + 1./4./u * ( z + (6-5*z)*u ) * Z2 
)
// q^1
- (1-z) * u2/up12 * (u + 4*z*(1+u)) * J27(y,z) 
+ 1./4./u/up12 * ( 2*z2 - (3-19*z+12*z2)*u - 2*(2-12*z+9*z2)*u2 )
   * lnz
- 1./2. * (1-z)/u/up12 * ( 6*z - (11-12*z)*u - 6*(2-z)*u2 ) 
   * log((1.+u)/r)
+ 1./2./u/up12/(1-u) * ( 2*z*(2-3*z) - (9-11*z+4*z2)*u 
 + (1+4*z2-z)*u2 + 2*(7-7*z+2*z2)*u3 + 2*(1-3*z+z2)*u4 - 2*z*u5 )
   * lnrpt
+ ( 1./4./u/(1+u)/(1-u) * ( 10*z2 -8*z + 17*(1-z)*u - (1-z+2*z2)*u2 )
  - t/r * (1-z)/u/up12/(1-u) * ( z - 2*(1-z)*u - (6-z)*u2 - 2*u3 ) )
   * lnu
// q^0
- 2. * r * (1-z) * (1-t*r) * (r-t)/up12 ) ;
 return temp;
}

// u -> 1 limit of qq_CD_ax_full:

double qq_CD_ax_u1(double z){
 double Z2,t,z2, lnz,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;
 lnz = log(z);  lnhf1pt = log((1.+t)/2.);
 temp = 1./3. * z * (1-z) * (
// q^2
+ 16 * z * J2_u1(z) - 8 * (2+z) * J3_u1(z) 
+ 8 * z * (5-2*z) * J21_u1(z)
+ 2 * (2-3*z) * li2(z) - 8 * (1-2*z) * li2(t)
+ 4 * (4-3*z) * li2(-t) - 8 * li2((1.-t)/2.)
- 3 * (1+z) * lnz*lnz + 2 * z * lnz * log(1.-z)
- 4 * z * (1+z)*(5-2*z)/(1-z)/(1+t)/(1+t) * lnz * log(1.+z)
+ 4 * (2-z) * lnz * log(1.+t) - 4 * (1-z) * lnhf1pt*lnhf1pt
+ 2./(1-z)/(1+t)/(1+t) * (z2*z+10*t*z2+4*z2-14*t*z+5*z+4*t+2) * log(2.) * lnz
+ 4 * (3-4*z) * Z2
// q^1
- 2 * (1+8*z) * J27_u1(z)
- (3-16*t-19*z+8*t*z+16*z2)/(1-z) * lnhf1pt
+ 1./2. * (-12*z2+8*t*z+43*z-16*t-7)/(1-z) * lnz
// q^0
- 4 * (5*t*z-10*t-z+1)/(1+t) );
 return temp;
}

// z -> 1 limit of qq_CD_ax_full:

double qq_CD_ax_z1(double y, double z){
  double a,a2,a3,a4,a5,a6,y2,y3,y4;
  a = 1.-z;  
  a2 = a*a;   a3 = a2*a;   a4 = a2*a2;   a5 = a3*a2;   a6 = a3*a3;  
  y2 = y*y;   y3 = y2*y;   y4 = y2*y2; 
return 7/12.*a2+a3*(-97/108.-47/108.*y+47/108.*y2)
  +a4*(155/864.*y-83/1728.-155/864.*y2)
  +a5*(-14581/43200.*y-151/21600.*y2+5527/43200.-4961/14400.*y4+4961/7200.*y3)
  +a6*(-3409/14400.*y+2639/21600.+17551/21600.*y3-1831/10800.*y2
       -17551/43200.*y4);
}

// Full, patched CD_ax answer, using the symmetry:

double qq_CD_ax_full(double y, double z){
 // z = 1 patch:
 if (std::fabs(1.-z) < 0.05){ return qq_CD_ax_z1(y,z); }
  // u = 1 patch:
 if (std::fabs(y-0.5) < 0.0001){ return qq_CD_ax_u1(z); }
 // (z = (2*u/(1+u))^2 patch (average values on either side of strip):
  double eps = 1.0e-07;
  if (std::fabs(y-y_sing(z)) < eps){ 
    return 0.5 * ( qq_CD_ax_full(y_sing(z)-2.*eps,z) 
		 + qq_CD_ax_full(y_sing(z)+2.*eps,z) );
  }
  if (std::fabs(y-(1.-y_sing(z))) < eps){ 
    return 0.5 * ( qq_CD_ax_full(1.-y_sing(z)-2.*eps,z) 
		 + qq_CD_ax_full(1.-y_sing(z)+2.*eps,z) );
  }
  return qq_CD_ax_full_1(y,z) + qq_CD_ax_full_1(1.-y,z);
}

// The full qq_CD_ax real "hard" term, including conversion factor:

double qq_CD_ax_real_hard(double y, double z){
 return real_conv(y,z) * qq_CD_ax_full(y,z) ; 
}
