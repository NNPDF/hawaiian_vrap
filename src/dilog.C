#include "dilog.h"

// routines only valid for real values of the argument; 
// imaginary parts not presently given.

double abs(double x){
  if (x < 0.0){ return -x; }
  else return x;
}

// this version uses 't Hooft and Veltman's change of variable
// good to ~ 10^(-16)

double li2(double x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1.) return PISQ6;
 if (x <= x_0){ 
   double temp = log(abs(1.0-x));
   return -li2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double z = - log(1.0-x);
   double temp = z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
                  *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
                  *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
                  *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
                   ))))))));
   return temp; }
   else if (x < x_2) return - li2(-x) + li2(x*x)/2.0 ;
   else { return PISQ6 - li2(1.0-x) 
                  - log(abs(x))*log(abs(1.0-x)) ; }
}

// COMPLEX version using 't Hooft and Veltman's change of variable

complex CLi2(complex x){
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1) return PISQ6;
 if (real(x) >= x_2) { return PISQ6 - CLi2(1.-x) - log(x)*log(1.-x) ; }
 if ((abs(imag(x)) > 1.) || (real(x)*real(x) + imag(x)*imag(x) > 1.2))
   return - CLi2(1/x) - 0.5 * log(-x) * log(-x) - PISQ6 ;
 if (real(x) <= x_0){ 
   complex zz = log(1.-x);
   return -CLi2(-x/(1.-x)) - zz*zz/2. ; }
 else if (real(x) < x_1){
   complex z = - log(1.-x);
   complex temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
                  *(1.-5.*z*z/294.*(1.-7.*z*z/360.
                  *(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
                  *(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)
                   ))))))));
   return temp; }
   else return - CLi2(-x) + CLi2(x*x)/2. ;
}

double ReLi2(complex x){ return real(CLi2(x)) ; }

double ImLi2(complex x){ return imag(CLi2(x)) ; }

// maple's definition

double dilog(double x){
 return li2(1.-x);
}

// the trilog, li3.
// good to ~ 10^(-14)

double li3(double x){
 double x_0 = -1.0;
 double x_1 = -0.85;
 double x_2 = 0.25;
 double x_3 = 0.63;
 double x_4 =  1.0;
 if (x == 1.) return ZETA3;
 if (x == -1.) return - 0.75 * ZETA3;
 if (x <= x_0){ 
   double lnx = log(-x);
   return li3(1.0/x) - PISQ6*lnx - lnx*lnx*lnx/6.0; }
 else if (x < x_1){
   return li3(x*x)/4.0 - li3(-x); }
   else if (x < x_2){
     double z = - log(1.0-x);
     double temp = z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
                    *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
                    *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
                    *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
                    *(1.0-53598548.0*z/524808375.0
                    *(1.0+22232925.0*z/107197096.0
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return li3(x*x)/4.0 - li3(-x); }
       else if (x < x_4){
         double ln1x = log(1.0-x); 
         return -li3(1.0-x) - li3(-x/(1.0-x)) + ZETA3 + PISQ6*ln1x
	   - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0; }
       else { 
         double lnx = log(x);
         return li3(1./x) + 2.0*PISQ6*lnx - lnx*lnx*lnx/6.0; }
}

// S12(x) -- only gets the real part correct for x > 1:

double S12(double x){
  double c;
  if (x > 1){ 
    return - li3(1.-x) + ZETA3 + log(x-1.) * li2(1.-x) 
           + 0.5 * log(x) * ( log(x-1.)*log(x-1.) - PISQ ); }
  else if (x==1){ return ZETA3; }
  else if ((0 < x) && (x < 1)){ 
    return - li3(1.-x) + ZETA3 + log(1.-x) * li2(1.-x) 
           + 0.5 * log(x) * log(1.-x)*log(1.-x); }
  else if (x==0){ return 0.; }
  else if (x<0){ 
    c = 1./(1.-x);
    return - li3(c) + ZETA3 + log(c) * li2(c) 
           + 0.5 * log(c)*log(c) * log(1.-c) 
           - 1./6. * log(c)*log(c)*log(c); }
  else { return 0.; }
}

// the tetralog, li4.
// good to ~ 1.7 * 10^(-12) (worst areas:  x = 0.9501, also x = - 0.9701)

double li4(double x){
 double x_0 = -1.0;
 double x_1 = -0.97;
 double x_2 = 0.25;
 double x_3 = 0.95;
 double x_4 =  1.0;
 if (x == -1) return -0.35 * PISQ6*PISQ6;
 if (x == 1) return 0.4 * PISQ6*PISQ6;
 if (x <= x_0){ 
   double lnx = log(-x);
   return - li4(1./x) - 0.5 * PISQ6*lnx*lnx
          - 1./24. * lnx*lnx*lnx*lnx - 0.7 * PISQ6*PISQ6; }
 else if (x < x_1){
   return li4(x*x)/8. - li4(-x); }
   else if (x < x_2){
     double z = - log(1.-x);
     double temp = z*(1.-7.*z/16.*(1.-151.*z/567.*(1.0-411.*z/2416.
                    *(1.-24986.*z/256875.*(1.-805.*z/49972.
                    *(1.+583406.*z/1159683.*(1.-7455.*z/137272.
                    *(1.+659921.*z/2444175.*(1.-251559.*z/2639684.
                    *(1.+24259894.*z/136410197.*(1.-30625595.*z/218339046.
                    *(1.+2134239258113.*z/16698772722450.
                    *(1.-1640443805715.*z/8536957032452.
                     )))))))))))));
     return temp; }
     else if (x < x_3){
       return li4(x*x)/8. - li4(-x); }
       else if (x < x_4){
        double y = 1.-x; 
        double lny = log(y);
        return 0.4*PISQ6*PISQ6 + ZETA3 * log(1.-y) + 0.5*PISQ6 * y*y 
         + (-11./36.+0.5*PISQ6+1./6.*lny) * y*y*y
         + (11./24.*PISQ6+1./4.*lny-19./48.) * y*y*y*y
         + (7./24.*lny+5./12.*PISQ6-599./1440.) * y*y*y*y*y
         + (137./360.*PISQ6+5./16.*lny-79./192.) * y*y*y*y*y*y
	 + (7./20.*PISQ6-3343./8400.+29./90.*lny) * y*y*y*y*y*y*y
	  + (363./1120*PISQ6-21977./57600.+469./1440.*lny) * y*y*y*y*y*y*y*y;
                        }
       else { 
         double lnx = log(x);
         return - li4(1/x) + PISQ6*lnx*lnx 
                - 1./24. * lnx*lnx*lnx*lnx + 0.8 * PISQ6*PISQ6; }
}

// the li2 below does not use 't Hooft and Veltman's change of variable
// -> not quite as good

double myli2(double x){
 int n = 25;
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == 1) return PISQ6;
 if (x <= x_0){ 
   double temp = log(abs(1.0-x));
   return -myli2(-x/(1.0-x)) - temp*temp/2 ; }
 else if (x < x_1){
   double temp = 0.0;
   double xk = 1.0;
   for (int k=1; k<=n; k++){ 
     xk = xk*x;
     temp = temp + xk/k/k; }
   return temp; }
   else if (x < x_2) return - myli2(-x) + myli2(x*x)/2.0 ;
   else { return PISQ6 - myli2(1.0-x) 
                  - log(abs(x))*log(abs(1.0-x)) ; }
}

// maple's definition

double mydilog(double x){
 return myli2(1.0-x);
}

// the 3-mass triangle, i3_3m(x,y,z), presently only for the pure 
// Minkowski region, x,y,z > 0:

double i3_3m(double s1, double s2, double s3){
  double x, y, lambda, rho, mylamsq, mylam, eta1, eta2, eta3, rho2;
  x = s1/s3;
  y = s2/s3;
  lambda = sqrt((1.0-x-y)*(1.0-x-y) - 4.0*x*y);
  rho = 2.0/(1.0-x-y+lambda);
  mylamsq =  s1*s1+s2*s2+s3*s3-2*(s1*s2+s2*s3+s3*s1);
  mylam = sqrt(mylamsq);
// Euclidean region:
  if (s1<0 && s2<0 && s3<0) {
    if (mylamsq > 0) { 
      double temp = 2.0 * ( li2(-rho*x) + li2(-rho*y) ) 
	            + log(abs(rho*x))*log(abs(rho*y)) 
	            + log(y/x)*log(abs((1.0+rho*y)/(1.0+rho*x)))
	            + PI_DEF*PI_DEF/3.0 ; 
      if (rho<0) temp = temp - PI_DEF*PI_DEF; 
      return  -temp/s3/lambda;
    }
    else {
       eta1 = (s1-s2-s3)/2.0;
       eta2 = (s2-s3-s1)/2.0;
       eta3 = (s3-s1-s2)/2.0;
       rho2 = sqrt(eta1*eta2 + eta2*eta3 + eta3*eta1);
       return  1.0/rho2 * ( fastCl(2.0*atan(rho2/eta1))
                          + fastCl(2.0*atan(rho2/eta2)) 
		          + fastCl(2.0*atan(rho2/eta3)) ); }
  }
  else if (s1>0 && s2>0 && s3>0) return -i3_3m(-s1,-s2,-s3);
  else return 0;
}

// Auxiliary function for the 3-mass triangle:

double fastCl(double x){
  double xb, x_lower, x_upper, x_switch;
  x_switch = 0.75*PI_DEF;
  x_lower = -x_switch;
  x_upper = 2.0*PI_DEF+x_lower;
// put x in the range [x_lower,x_upper]:
  if (x > x_upper) return fastCl(x-2.0*PI_DEF); 
  if (x < x_lower) return fastCl(x+2.0*PI_DEF);
// dividing line for the two expansions is x_switch: 
  if (x < x_switch){          // expansion around 0:
     return  -x*log(abs(x)) + x * ( 1.0 + x*x/72.0 * ( 1.0 + x*x/200.0 
           * ( 1.0 + x*x*5.0/441.0 * ( 1.0 + x*x*7.0/480.0
	   * ( 1.0 + x*x*2.0/121.0 ) ) ) ) ); }
  else {
    xb = x-PI_DEF;
    return -xb*log(2.0) + xb*xb*xb/24.0 * ( 1.0 + xb*xb/40.0  
           * ( 1.0 + xb*xb/21.0 * ( 1.0 + xb*xb/288.0
           * ( 1.0 + xb*xb*62.0/55.0 * ( 1.0 + xb*xb*691.0/9672.0
           * ( 1.0 + xb*xb*5461.0/72555.0 ) ) ) ) ) );
  }
}      
