/* =============================================================
   Special functions required for the "hard" NNLO contributions 
   to the DY rapidity distributions.
   ============================================================= */

#include "hardfns.h"

// For the "real" phase space integration, factor needed
// to go from Kirill's conventions to my program's conventions:

double real_conv(double y, double z){
  return (1.+z)/4./z/z/(1. + (1.-z)*(1.-z)/z * y*(1.-y)) ;
}

//====================================================================
// "plus" distributions (without the subtraction):

double DD1(double n, double y){ return pow(log(y),n)/y ; }

//====================================================================
// The J functions:

// J3:

double J3(double y, double z){
 double u, t, r, d1, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));    t = sqrt(z);  r = sqrt(u);
 d1 = sqrt(u+4.*z*(1.+u));
 temp = 1./(1.+u) * ( - 0.25 * log(u/z)*log(u/z)
     - 0.25 * ( - 2. * log(std::fabs(d1-r-2*t*(1+u))/(d1-r))
               - 2. * log((d1+r+2*t*(1+u))/(d1+r))
               + 2. * log(1+t*r) + log(z/u) ) * log(1+u)
// Im part comes from left 2 polylogs; just set Im part to zero:
     + li2(-2.*t*(1.+u)/(r-d1)) + li2(-2.*t*(1.+u)/(d1+r))
     - li2(-2.*z*(1.+u)/r/(r-d1)) - li2(-2.*z*(1.+u)/r/(d1+r)) );
 return temp;
}

// Functions needed for J2:

double M1(double x, double xp, double xm){
 double l1,l2,l3,l4,temp;
 l1 = log(std::fabs((x+xp)/(x-xp)));  l2 = log(std::fabs((x+xm)/(x-xm)));
 l3 = log(std::fabs((x+xp)/2/xp));    l4 = log(std::fabs((x+xm)/2/xm));
 temp =  - li2((x+xp)/2./xp) + li2((x+xm)/2./xm)
         - li2(-(x-xp)/(xp-xm)) + li2((x+xp)/(xp-xm))
         + 0.25 * (l1*l1-l2*l2) - 0.5 * (l3*l3-l4*l4)
         + 0.5 * l1 * log(std::fabs((x+xm)*(x-xm)/(xp-xm)/(xp-xm))) ;
 return temp;
}

double I_analytic(double y, double z){
 double u,r,r1,r2,d1,xp,xm,xmax,xmin,term1,p1,p2,l1,l2,term2;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   r = sqrt(u);
 r1 = sqrt(1.+u); r2 = sqrt(5.+4.*u);   d1 = sqrt(u+4.*z*(1.+u));
 xp = -(2.+u-2.*r1)/u;   xm = 1./xp;
 xmax = r2 + 2.*r1;   xmin = (d1 + 2.*sqrt(z*(1.+u)))/r;
 term1 = 0.5/u * (M1(xmax,xp,xm)-M1(xmin,xp,xm));
//
 p1 = - (2.+u) + u*r2;    p2 = r*d1;
 l1 = log(std::fabs(-(2.+u+u*r2)/p1));
 l2 = log((2.+u+p2)/(2.+u-p2));
 term2 = 0.75/u * ( 
     l1 * log(std::fabs(u*(u+2.)*(r2-1.)/p1/(1.+u)))
   - log(std::fabs(p1*(r+d1)/r/(1.+r2)/(2.+u-p2))) * log(1.+u)
   - l2 * log(r*(u+2.)*(d1-r)/(2.+u-p2)/(1.+u))
   - 1./6. * ( l1*l1 - l2*l2 )
   + li2(u*(u+2.)*(1.+r2)/p1) - li2(-r*(u+2.)*(r+d1)/(2.+u-p2))
   - li2((2.+u+p2)/(2.+u-p2)/(1.+u)) + li2(-(2.+u+u*r2)/p1/(1.+u))
 );
//
 return term1 + term2;
}

// J27 (also an auxiliary function for J2):

double J27(double y, double z){
 double u, t, r, r1, r2, d1;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 r1 = sqrt(1.+u);   r2 = sqrt(5.+4.*u);  d1 = sqrt(u+4.*z*(1.+u));
 return - 0.5/r/d1 * ( 
	   log(std::fabs((r2-r1 + r1*r2-3.-2.*u)/(r2-r1 - r1*r2+3.+2.*u)))
         + log(std::fabs((d1-2.*t*r1+r1*r - r1*d1+2.*t+2.*t*u+r)
	          /(d1-2.*t*r1+r1*r + r1*d1-2.*t-2.*t*u-r)))        
         + 3. * log(std::fabs((r - d1+2.*t*r1)*(1. + r2-2.*r1)
                       /(r + d1-2.*t*r1)/(1. - r2+2.*r1))) );
}

// finally, J2:

double J2(double y, double z){
 double u,r,d1,term0,term1, term2;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   r = sqrt(u);    d1 = sqrt(u+4.*z*(1.+u));
 term0 = - J3(y,z);
 term1 = - d1/r * log(std::fabs((2.+u+r*d1)/(2.+u-r*d1))) * J27(y,z);
 term2 = I_analytic(y,z);
 return term0 + u/(1.+u) * (term1 + term2) ;
}

// J1 is not independent, but here it is anyway:

double J1(double y, double z){
 double u;
 u = (z+y*(1.-z))/(1.-y*(1.-z));
 return - (1.+u)/u * J3(1.-y,z)  
         + li2(-u) + log(u) * (log(1.+u) - 0.5 * log(u)) + 0.5*PISQ6 ;
}

// J21:

double J21(double y, double z){
 double u, t, r, x1, tempsoft, tempdiff; 
 complex rtx1, a1p, a1m, temp;
 u = (z+y*(1.-z))/(1.-y*(1.-z));   t = sqrt(z);   r = sqrt(u);
 x1 = 4.*u*u/z/(1.+u)/(1.+u) - 1.;
 if (x1 < 0) { rtx1 = I*sqrt(-x1); }
 else { rtx1 = sqrt(x1); }
 a1p = 2./z/(1.+1./u)/(1.+rtx1*I);   
 a1m = 2./z/(1.+1./u)/(1.-rtx1*I);
 temp = I/z/(1.+u)/rtx1 * (
      - log(t*a1p) * log((a1p-1.)/(1.-t*r*a1p))
      + log(t*a1m) * log((a1m-1.)/(1.-t*r*a1m))
      - CLi2(1.-a1p) + CLi2(1.-t*r*a1p) 
      + CLi2(1.-a1m) - CLi2(1.-t*r*a1m) );
 tempsoft = - (1.+z)*(1.-u*z)*log((1.+z)/2.)*log(z)
              /(1.+u)/r/(r+t)/(1.+t*r)/(1.-z)/(1.-z);
 tempdiff = real(temp) - tempsoft;
 return tempdiff;
}

// u -> limits of selected functions:

double J3_u1(double z){ return J3(0.5,z); }

double J27_u1(double z){
 double t,x1;
 t = sqrt(z);
 x1 = sqrt(1.+8.*z);
 return - 0.25/x1 * ( - 3 * log((x1+1.)/(x1-1.))
                      + log((3.*x1-1.+8.*t)/(3.*x1+1.-8.*t)) );
}

double J2_u1(double z){
 double t, term0, term1, x1, cp, cm, dp, dm, rt2, rt8, tmrt8, ln2, term2;
 t = sqrt(z);
 x1 = sqrt(1.+8*z);   rt2 = sqrt(2.);   rt8 = 2.*rt2;  tmrt8 = 3.-rt8;
 ln2 = log(2.);
 cp = x1+rt8-3.+2.*t*rt2;   
 cm = x1-rt8+3.+2.*t*rt2;
 dp = 3.*x1+6.*t*rt2+1.-8.*t-rt8*x1;
 dm = 3.*x1+6.*t*rt2-1.-8.*t-rt8*x1;
 term0 = - J3_u1(z);
 term1 = - x1 * log((3.+x1)/(3.-x1)) * J27_u1(z);
 term2 = 
   - 0.5 * li2(-dm/2.)
   + 0.5 * li2(-(dp-4*(3*rt2-4))/4/(3*rt2-4)) 
   - 0.5 * li2((dm+4*(3*rt2-4))/4/(3*rt2-4)) 
   + 0.5 * li2(-cp/2./tmrt8)
   - 0.75 * li2(3.*(x1+1)/(x1-3.)) 
   - 0.75 * li2(-(3.+x1)/2./(x1-3.))
   + 0.5 * ln2 * log(std::fabs(dm/cp)) + 0.5 * log(cp/cm) * log(tmrt8*2.*rt8)
   - 0.5 * log(tmrt8) * log(cp) - 0.25 * log(cp/cm) * log(std::fabs(dp*dm))
   + 0.25 * log(std::fabs(cp/dm)) * log(std::fabs(cp*dm)) 
   + 0.125 * log(std::fabs(dm/dp*cm/cp))*log(std::fabs(dm/dp*cp/cm))
   + 0.125 * log((3.+x1)/(3.-x1)) * log((3.+x1)/(3.-x1)) 
   - 0.75 * ln2 * log(4./3.*(x1+1.)/(3.-x1))
   - 0.75 * log((3.+x1)/(3.-x1)) * log(3./2.*(x1-1.)/(3.-x1))
   - 0.5 * li2(-3/rt8) - 0.5 * li2(-rt8/tmrt8) 
   - 3./16. * ln2*ln2 - 0.75 * ln2*log(3.) + 0.5 * log(6.) * log(tmrt8)   
   + 1.25 * PISQ6 ;
 return term0 + 0.5 * (term1 + term2);
}

double J21_u1(double z){ return J21(0.5,z); }

// u -> z limits:

double J2rem(double z){
 double u, r, r1, r2, bp, bm, cp, cm, dp, dm, ee, gg, temp;
 u = 1./z;  r = sqrt(u);  r1 = sqrt(1.+u);  r2 = sqrt(5.+4.*u);
 bp = u*r2+2.*u*r1+2.+u-2.*r1;   bm = u*r2+2.*u*r1-2.-u+2.*r1;
 cp = 2.+u+u*r2;   cm = -2.-u+u*r2;
 dp = 2.+2.*u+u*r1+2.*r1;   dm = -2.-2.*u+u*r1+2.*r1;
 ee = 2.*r2+u*r2-2.*r2*r1+4.*r1+2.*u*r1-4.-5.*u;
 gg = 2.+u-2.*r1;
 temp = + 0.5 * li2(bm*gg/4./dm) - 0.5 * li2(-bp*gg/4./dm)
        + 0.5 * li2(-ee/2./u) - 0.5 * li2(-bm/2./gg)
        - 0.5 * li2(gg*r1/dm) + 0.5 * li2(-(u+2)*gg/2./dm)
        + 0.5 * li2(-2.*r1/gg)
        + 0.75 * li2(-cp/cm/(1.+u)) + 0.75 * li2(u*(u+2.)*(1.+r2)/cm)
//
  - 0.75 * log((u+2.)*(u+2.)/2./r/(1.+u)) 
         * log((u+2.)*(u+2.)/(1.+u)/(1.+u)/2./r)
  - 0.25 * log(bm/2./gg)*log(bm/2./gg) 
  - 0.125 * log(ee/(ee+2.*u))*log(ee/(ee+2.*u))
  + 0.125 * log(bm/bp)*log(bm/bp) + 0.25 * log(ee/2./u)*log(ee/2./u)
  + 0.25 * log(u*u*(ee+2.*u)*ee/16./dm/dm) * log(bm/bp)
   //
  + 0.125 * log((u+2.)*(u+2.)/2./r/(1.+u))
          * log((u+2.)*(u+2.)/2./r/(1.+u)) 
  + 0.375 * log((u+2.)*(u+2.)/(1.+u)/(1.+u)/2./r)
          * log((u+2.)*(u+2.)/(1.+u)/(1.+u)/2./r)
  + 0.375 * log((u+2.)*(u+2.)/2./r) * log((u+2.)*(u+2.)/2./r)
  - 0.75 * log(1.+u) * log(cm*(u+2.)/(1.+r2)/2./r/u)
   //
  - 0.125 * log(cp/cm)*log(cp/cm) 
  - 0.125 * log(2.*r1/(u+2.))*log(2.*r1/(u+2.)) 
  + 0.25 * log(2.*r1/gg)*log(2.*r1/gg)
  - 0.125 * log(dp/r/u*gg/(u+2.))*log(dp/r/u*gg/(u+2.))
  - 0.25 * log(r*u*u*dp*gg/4./dm/dm/(u+2.)) * log(2.*r1/(u+2.))
  + 0.75 * log(cp/cm) * log(u*(u+2.)*(r2-1.)/cm/(1.+u)) - PISQ/8. ;
 return temp;
}

// There is also a mild (?) singularity for z = (2*u/(1+u))^2,
// or y = y_sing(z), where

double y_sing(double z){ return (sqrt(z)*(1.+z)-2.*z)/2./(1-z); }
