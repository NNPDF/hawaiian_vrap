/* =============================================================
   The functions required for the gg initial-state 
   NNLO contributions to the DY rapidity distributions.
   ============================================================= */

#include "hardfns.h"

double gg_full_1(double y, double z){
  double u,t,r,Z2,z2, u2,u3,u4,u5,u6,u7, up12,up13,up14,up15, 
    ln2,lnz, lnr,lnu,lnopu,lnhfopu,lnumo, lnrpt,lnrmt, lnoptr,lnomtr, 
    lnrmo, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;
 z2 = z*z;  u2 = u*u;  u3 = u2*u;  u4 = u2*u2;  u5 = u3*u2;
 u6 = u3*u3;  u7 = u4*u3;
 up12 = (1.+u)*(1.+u);   up13 = up12*(1.+u);
 up14 = up12*up12;    up15 = up13*up12;  
 ln2 = log(2.);  lnz = log(z);  lnr = log(r);  lnu = log(u); 
 lnopu = log(1.+u);  lnhfopu = log((1.+u)/2.);  lnumo = log(abs(u-1.));
 lnrpt = log(r+t);  lnrmt = log(r-t);
 lnoptr = log(1.+t*r);  lnomtr = log(1.-t*r);  lnrmo = log(abs(r-1.));  
 temp = 
// q^2:
// J's:
+ 1./8. * z * (1-z)/u2/(1+u) * ( 1 - 2*u*z + u2*(1+2*z+2*z2) )
    * ( (1+u) * J2(y,z) - 1./4. * li2(2.*u/(1.+u)) 
      + 3./4. * li2((1.-u)/2.) + 1./2. * ln2 * lnomtr )
+ 1./8. * z * (1-z)/u2 * ( 1 - 2*u*z - 7*u2*(1+2*z+2*z2) )
    * J3(y,z)
// polylog's
+ 1./4. * z2 * (1-z) * (1-u) * (1+u)/u2 * ( 2*z + u*(1+2*z) )
    * li2(-u)
- 1./16. * z * (1-z) * (1-u)/u2/(1+u) * ( 1 + u*(1-2*z) + u2*(1+z2) )
    * li2(-r)
- 1./64. * z * (1-z)*(1-u)/u2/(1+u) * ( 1+(u+u3)*(1-2*z)+2*u2*(1+z2)+u4 )
    * ( 3./2. * li2(-u) + li2((r-1.)/r) 
      - 3./2. * li2(-1./u) - li2(1.-r) )
//
- 1./4. * z * (1-z)/u2/(1-u) * ( 1-2*z2 - 2*u*z*(2+z) + u2*(1+6*z2) )
    * li2(t*r)
+ 1./2. * z * (1-z) * u/(1+u) * ( 2+5*z+6*z2 + u*z*(1+3*z) + u2*z2 )
    * li2(-t*r)
- 1./16. * z * (1-z)/u2/(1-u)/(1+u) * ( 3 + u*(1-6*z) + u2*(3-4*z+6*z2)
  + u3*(1-6*z+2*z2) )
    * ( li2((1.-t*r)/2.) - li2(t*r) - li2(-t*r) )
- 1./8. * z * (1-z)/u2/(1-u) * ( 1-4*z2 - 2*u*z*(7+2*z) - 3*u2*(1+2*z+2*z2)
   - u3*(12-24*z-16*z2) - 8*u4*(2+z) + 32*u5 )
    * li2((1.-t*r)/(1.+u))
//
+ 2. * t/r * z * (1-z) * (1+2*z) * (1+u)
    * ( li2(t*r) - li2(-t*r) + li2((1.-t*r)/2.)
      + 1./2. * li2((1.-t*r)/(1.+u)) 
      + 1./4. * lnoptr*lnoptr - 1./4. * lnomtr*lnomtr 
      - 1./2. * lnoptr * ( lnz + lnopu + lnr - lnrmt - lnomtr )
      + 1./2. * lnomtr * ( lnz - 1./2. * lnrmt + 5./2. * lnr
                         - 2. * ln2 ) 
      + 1./4. * lnopu*lnopu - 1./8. * lnu*lnu - 5./4. * Z2 )
// pure ln's
+ 1./16. * z * (1-z)/u2/up14 * ( 6*z2 + 2*z*(7+18*z)*u
  + 2*(3+52*z2+38*z)*u2 + (7+152*z2+123*z)*u3 + (87*z2+106*z+10)*u4 )
    * lnz*lnz
- 1./8. * z2 * (1-z) * (1-u) * (1+u)/u2 * ( 6*z + u*(7+6*z) )
   * lnz * lnu
+ 1./4. * z * (1-z)/u2 * ( 4*(4+z2) + (10*z+3+10*z2)*u2 + 4*(2*z2+z+2)*u )
    * lnz * lnopu
+ 1./256. * z * (1-z)/u2/up14 * ( 5+96*z2 + (15+512*z2+182*z)*u
    + 2*(406*z+13+536*z2)*u2 + (1792*z2+2490*z+321)*u3
    + (161+1140*z+992*z2)*u4 )
    * lnu*lnu
+ 1./64. * z * (1-z)/u2 * ( 133+16*z2 + (67+32*z2+6*z)*u 
   + (58*z+21+69*z2)*u2 )
    * lnopu*lnopu
+ 3./32. * z * (1-z)/u2 * ( 1 - (1+2*z)*u + (1+z)*(1+z)*u2 )
    * lnumo * lnhfopu
- 1./8. * z * (1-z)/u2/(1-u)/(1+u) * ( 1+4*z2 - (1-10*z-8*z2)*u 
    + (4*z2+12*z+3)*u2 - (5+18*z2+18*z)*u3 )
    * lnu * lnopu
- 1./64. * z*(1-z)/u2 * (3 + (5-6*z)*u + (1-3*z)*(3-z)*u2 )
    * ln2*ln2
+ 1./64. * z * (1-z) * (1-u)/(1+u)/u2 * ( 2 + 2*(1-2*z)*u + (5+5*z2+6*z)*u2 )
    * lnrmo * lnu
+ 3./64. * z * (1-z) * (1+z)*(1+z) * (1-u)/(1+u) 
    * log((r+1.)/2.) * lnu
+ 1./32. * z * (1-z)/u2/(1+u) * ( 8*u3+41*u2+96*u+48*u3*z2
  +56*u3*z+34*z2*u2-18*u*z+24*u4*z+18*u2*z+24*u4*z2+8*u5*z2+65 )
    * lnoptr*lnoptr
+ 1./32. * z * (1-z)/u2/(1-u)/(1+u) * ( 65*u3+115*u2+16*u4-63*u
 -46*u3*z2-94*u3*z-58*z2*u2+26*u*z+64*u4*z-68*u2*z+40*u4*z2
 +56*u5*z+48*u5*z2+24*u6*z2-125 )
    * lnz * lnoptr
+ 1./16. * z * (1-z)/u2/(1-u)/(1+u) * ( 31*u3+53*u2+8*u4-33*u
 -18*u3*z2-34*u3*z-38*z2*u2+22*u*z+32*u4*z-28*u2*z+24*u4*z2
 +24*u5*z+16*u5*z2+8*u6*z2-67 )
    * (lnopu + lnr) * lnoptr
+ 1./4. * z * (1-z)/u2/(1-u)/(1+u) * ( u6*z2+2*u5*z2-u4*z2
 -2*u3*z2+2*z2*u2-2*u*z+u5*z-3*u3*z+1+u2 )
    * lnu * lnoptr
+ 1./2. * z * (1-z)/u2 * ( 8+4*u-2*u*z+4*z2*u2+u2+4*u2*z
 +2*u3*z2+3*u3*z+u4*z2 )
    * ( lnrmt + lnomtr ) * lnoptr
+ 1./8. * z * (1-z)/u2/(1-u) * (1-2*u*z+u2*(1-2*z+2*z2))
    * ln2 * lnoptr
+ 1./4. * z * (1-z)/u * (4*z2*u+2*u2*z+2*z+u+4*u*z)
    * (lnrmt + lnomtr) * lnomtr 
+ 1./32. * z * (1-z)/u2/(1-u)/(1+u) * ( -3*u3-17*u2+16*u4-3*u
 -22*u3*z2-6*u3*z-66*z2*u2-30*u*z+64*u4*z-52*u2*z+56*u4*z2
 +40*u5*z+16*u5*z2+8*u6*z2-1 )
    * lnz * lnomtr
- 1./16. * z * (1-z)/u2/(1+u) * ( 8*z2+96*u+96*u4+48*u3+49*u2
 +64*u5+65+80*u3*z2+72*u3*z+74*u2*z+82*z2*u2+6*u*z+24*z2*u
 +8*u5*z2+24*u4*z2+8*u4*z )
    * lnopu * lnomtr
- 1./32. * z * (1-z)/u2/(1-u)/(1+u) * ( -16*z2+64*u5-63*u3-9*u2
 -116*u4+u+128*u6+18*u3*z2+66*u3*z-26*z2*u2-46*u*z-32*z2*u
 +48*u4*z-52*u2*z+40*u4*z2-32*u5*z+16*u5*z2+8*u6*z2+3 )
    * lnu * lnomtr
+ 1./32. * z * (1-z)/u2 * ( 5 + (7-50*z-16*z2)*u + 6*(1-5*z-3*z2)*u2 )
    * Z2
// q^1:
- 1./16. * z * (1-z) * (u+4*z*(1+u)) /u/up13
  * ( - 2 + 21*u + 46*u2 + 14*u3 + 36*u*(1+u)*z )
     * J27(y,z)
- 1./32. * z * (1-z)/u2/up13 * ( 64*z2-4*u+588*u4+20*u3-30*u2
 +1437*u5+1278*u6+392*u7+280*u3*z2+376*u3*z+126*u2*z
 +310*z2*u2-2*u*z+224*z2*u-36*u6*z+76*u5*z2+206*u4*z2
 -8*u5*z+276*u4*z )
    * lnrpt
+ 1./4. * z * (1-z)/u/up12 * ( 16+31*u3+28*u2+16*u4+31*u
 -42*u2*z-16*u*z+3*z+2*u4*z2+3*u4*z-13*u3*z2-16*u3*z
 -32*z2*u2-13*z2*u+2*z2 )
    * lnrmt 
- 1./128. * z * (1-z)*(1-u)/u2/up13 * ( 280*u3*z+3577*u3
 +232*u3*z2+6184*u2+304*z2*u2+624*u2*z+3328*u-16*u*z+784 )
    * lnu
+ 1./64. * z * (1-z)/u2/up12 * ( 260*u3*z2+452*z2*u2+320*z2*u
 +128*z2+496*u3*z+312*u2*z-76*u*z+784+1050*u2+83*u3+1764*u )
    * lnopu
+ 1./64. * z * (1-z)/u2/up14 * ( 64*z2-250*u-932*u4-2023*u3-1017*u2
 +1280*u3*z2+666*u3*z-104*u2*z+426*z2*u2-94*u*z+224*z2*u
 +690*u4*z2+676*u4*z )
    * lnz 
+ 1./8. * t/r * z * (1-z) * (1+u)/u * ( 
    (1-4*u-34*u*z+u2) * log((r+t)/2./t)
  - 4 * (8*u2+u2*z-11*u*z-7*u+8+z) * log((r-t)/t)
   - (1-u)/up15 * ( 1+4*z - 2*(2-5*z)*u - (59+32*z)*u2 + 16*z*u3 )
    * lnu )
// q^0:
- 1./32. * z * (1-z)/u/up14 * ( -523*u3*z2-235*z2*u2-42*z2*u
 +48*z2-399*u3*z+286*z+802*u2*z+731*u*z+528+3195*u2+1407*u3+1938*u )
+ 1./16. * t/r * z * (1-z)/u/up13 * ( -62*u*z-202*u2*z-332*u3*z
	+24*z+824*u+1097*u2+321*u3+264 ) 
;
 return temp;
}

// u -> 1 limit of gg_full:

double gg_u1(double z){
 double Z2,t,z2,ln2,lnt,ln1mt,lnhf1pt, temp;
 Z2 = ZETA2;   t = sqrt(z);   z2 = z*z;
 ln2 = log(2.);  lnt = log(t);  ln1mt = log(1.-t);  lnhf1pt = log((1.+t)/2.);
 temp = 1./128. * z * (1-z) * (
// q^2
+ 64 * (z2+1) * J2_u1(z) - 64 * (3+8*z+7*z2) * J3_u1(z) 
+ 16 * (43*z2-128*t*z+16*z-64*t+15) * li2(-t)
+ 16 * (51*z2+192*t*z+56*z+96*t+91) * li2((1.-t)/2.) 
+ 16 * (23*z2+128*t*z+16*z+64*t-5) * li2(t)
+ 128 * (1-2*t+2*z)*(1-2*t+2*z) * ln1mt*ln1mt
+ 4 * (385*z2+319*z+23) * lnt*lnt
+ 256 * (7*z2+13+5*z+4*t+8*t*z) * ln1mt * lnhf1pt
+ 128 * (27+22*z2+14*z) * ln2 * lnt
- 16 * (89*z2-128*t*z+128*z-64*t+21) * ln1mt * lnt
- 16 * (133*z2+128*t*z+96*z+64*t+201) * log(1.+t) * lnt 
+ 8 * (105+40*z+57*z2+32*t+64*t*z) * lnhf1pt* lnhf1pt
- 8 * (35*z2+320*t*z+80*z+160*t-17) * Z2
// q^1
- 2 * (79+72*z)*(1+8*z) * J27_u1(z)
+ 32 *(27*z-72*t+61)*(1-z) * ln1mt
+ (1342*z2-64*t*z+572*z+2368*t-2111) * lnt
- (64*t+2240*t*z+3617+796*z+1160*z2) * lnhf1pt
// q^0
+ 2 * ( -355*z+188*z2-1767-572*t*z+2506*t ) 
);
  return temp;
}

// u -> z limit of gg_full:

double gg_uz(double y, double z){
 double u,r,t,Z2,z2,ln2,lnz,ln1mz,ln1pz,lnhf1pz, lnde,  temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;  z2 = z*z;
 ln2 = log(2.);  lnz = log(z);  ln1mz = log(1.-z);  
 ln1pz = log(1.+z);   lnhf1pz = log((1.+z)/2.);   lnde = log(r-t);
 temp = 1./32. * z * (1-z) * (
+ 8 * (1-2*z+2*z2) * ( lnde*lnde + 2 * ln1mz * lnde 
                      + ln1mz*ln1mz - 3 * lnz * ln1mz ) 
+ 8 * (-1+10*z+6*z2) * lnz * lnde 
- 16 * (3+4*z+6*z2) * lnhf1pz * lnde 
- 16 * (1-z) * (1-5*z) * (ln1mz + lnde) 
+ 4 * z * (1+z) * ( J2rem(z) - log(1.+2.*z) * log((1.+z)/2./z2) )
- 4 * (4+7*z+23*z2) * li2(z) 
+ 16 * (5+6*z+6*z2) * li2(-z)
+ 2 * (24+31*z+47*z2) * li2((1.+z)/2.)
+ 2 * (1-37*z-31*z2) * lnz*lnz
+ 2 * (52+65*z+69*z2) * lnz * lnhf1pz 
+ z * (1+5*z) * lnhf1pz*lnhf1pz + 8 * (10+11*z+11*z2) * ln2 * lnz
- 2 * (24+29*z+45*z2) * ln2 * lnhf1pz + 8 * (1-27*z-14*z2) * lnz
- 2 * (4-23*z-7*z2) * Z2 + 8 * lnhf1pz * (29 * z2 + 11 * z + 4)
- 208 * z * (1-z) );
  return temp;
}

// z -> 1 limit of gg_full:

double gg_z1(double y, double z){
  double a,a2,a3,a4,a5,y2,y3,y4,lna,lna2,lny,lny2,ln1my,ln1my2;
  a = 1.-z;  
  a2 = a*a;   a3 = a2*a;   a4 = a2*a2;   a5 = a3*a2;  
  y2 = y*y;   y3 = y2*y;   y4 = y2*y2; 
  lna = log(a);  lny = log(y);  ln1my = log(1.-y);
  lna2 = lna*lna;   lny2 = lny*lny;  ln1my2 = ln1my*ln1my; 
return 
   a*(-1/2.*ZETA2+1/4.*ln1my2+1/2.*lny*ln1my+lna*ln1my+lna2+lna*lny+1/4.*lny2)
 +a2*(3/2.*ZETA2-3/4.*ln1my2+3/2.*ln1my-3/2.*lny*ln1my-3*lna*ln1my-3/4.*lny2
   +3/2.*lny-3*lna*lny-3*lna2+3*lna-3/4.)
 +a3*(-(-y+2+y2)*ZETA2+1/2.*(-y+2+y2)*ln1my2-1/2.*(2*y2+7-2*y)*ln1my
   +(-y+2+y2)*lny*ln1my+2*(-y+2+y2)*lna*ln1my+1/2.*(-y+2+y2)*lny2
   -1/2.*(2*y2+7-2*y)*lny+2*(-y+2+y2)*lna*lny+2*(-y+2+y2)*lna2
   -(2*y2+7-2*y)*lna-7/4.*y+7/4.*y2+105/32.)
 +a4*(1/2.*(-y+2+y2)*ZETA2-1/4.*(-y+2+y2)*ln1my2
   +1/24.*(21*y2+65-21*y)*ln1my-1/2.*(-y+2+y2)*lny*ln1my-(-y+2+y2)*lna*ln1my
   -1/4.*(-y+2+y2)*lny2+1/24.*(21*y2+65-21*y)*lny-(-y+2+y2)*lna*lny
   -(-y+2+y2)*lna2+1/12.*(21*y2+65-21*y)*lna+119/96.*y-119/96.*y2-1933/576.)
 +a5*(-7/4.*y2*(-1+y)*(-1+y)*ZETA2+7/8.*y2*(-1+y)*(-1+y)*ln1my2
   -1/96.*(30*y-564*y3+282*y4+252*y2+35)*ln1my+7/4.*y2*(-1+y)*(-1+y)*lny*ln1my
   +7/2.*y2*(-1+y)*(-1+y)*lna*ln1my+7/8.*y2*(-1+y)*(-1+y)*lny2
   -1/96.*(30*y-564*y3+282*y4+252*y2+35)*lny+7/2.*y2*(-1+y)*(-1+y)*lna*lny
   +7/2.*y2*(-1+y)*(-1+y)*lna2-1/48.*(30*y-564*y3+282*y4+252*y2+35)*lna
   -23/384.*y-1191/128.*y3+3965/4608.+3619/768.*y2+1191/256.*y4);
}

// Full, patched answer, using symmetry:

double gg_full(double y, double z){
 // z = 1 patch:
 if (abs(1.-z) < 0.01){ return gg_z1(y,z); }
  // u = 1 patch:
 if (abs(y-0.5) < 0.0001/sqrt(1.-z)){ return gg_u1(z); }
 // u = z patch:
 if (y < 0.00000001/(1.-z)){ return gg_uz(y,z); }
 // u = 1-z patch:
 if (1.-y < 0.00000001/(1.-z)){ return gg_uz(1.-y,z); }
 // (NO z = (2*u/(1+u))^2 patch REQUIRED HERE!)
 return gg_full_1(y,z) + gg_full_1(1.-y,z);
}

//================================================================
// The mu-dependent gg terms

double gg_mu_1(double y, double z, double muFQ){   
  double u,t,r,Z2,z2, u2,u3, lnmu;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 Z2 = ZETA2;     z2 = z*z;  u2 = u*u;  u3 = u2*u;
 lnmu = 2.*log(muFQ);
 return 1./4. * z * (1.-z)/u * ( 
1./2. * ( 2*z + (1+4*z+4*z2)*u + 2*z*u2 )
      * ( lnmu*lnmu - 4 * lnmu * log((u-z)/t) )
- t/r * (1+2*z) * u * (1+u)
    * ( lnmu*lnmu - 4 * lnmu * log(2.*(r-t)/(r+t)) )
- 2 * (1+2*u)/u * ( z2 + z*u - 2*z*u2 + 4*u3 )
    * lnmu * log(r+t)
+ lnmu * (
+ z/u * ( z + (1+2*z)*u ) * log(z*u)
+ 2/u * ( 8+z2 + (4+z+2*z2)*u + (1+2*z)*(1+2*z)*u2 )
    * log((1.+u)/r)
- 1/(1+u)/(1+u) * ( 16+3*z+2*z2 + (31-16*z-13*z2)*u + (14-21*z-16*z2)*u2 )
+ t/r * (1+u) * ( 16+2*z - (7+11*z)*u ) ) );
}

double gg_mu(double y, double z, double muFQ){
  return gg_mu_1(y,z,muFQ) + gg_mu_1(1.-y,z,muFQ);
}

// The full gg real "hard" term, including conversion factor and
// mu-dependence:

double gg_real_hard(double y, double z, double muFQ){
 return real_conv(y,z) * ( gg_full(y,z) + gg_mu(y,z,muFQ) ) ; 
}
