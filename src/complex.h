
#ifndef COMPLEX_H
#define COMPLEX_H


#include <math.h>
#include <iostream>
#include "thbasics.h"

class complex{

   private:
      double r, i;

   public:
      inline complex(double re=0.0, double im= 0.0): r(re), i(im) {}
      inline complex(const complex & a): r(a.r), i(a.i){}

      friend inline double real(const complex & a);
      friend inline double imag(const complex & a);
      friend inline double norm(const complex & a);
      friend inline double arg(const complex & a);
              /*   should be in the range  -PI to PI  */
      friend inline double ABS(const complex & a);

      inline complex & operator = (const complex & a){
         if (r != a.r || i != a.i) {
            r = a.r; i = a.i;
         }
         return *this;}
      inline complex  operator - () const {complex A = *this; A.r = -r; 
                       A.i = -i; return A;}
      inline complex& operator += (double a){ r += a; return *this;}
      inline complex& operator -= (double a){ r -= a; return *this;}
      inline complex& operator *= (double a){ r *= a; i *= a;
         return *this;}
      inline complex& operator /= (double a){ r /= a; i /= a;
         return *this;}
      inline complex& operator += (const complex &a) {
         r +=a.r ; i += a.i; return *this;}
      inline complex& operator -= (const complex &a) {
         r -=a.r; i -= a.i; return *this;}
      inline complex& operator *= (const complex &a) { double rr = r;
         r = rr*a.r - i*a.i; i = i*a.r + rr*a.i; return *this;}
      inline complex& operator /= (const complex &a){
	      double m = norm(a); double rr = r;
	      r = (rr * a.r+ i*a.i)/m; i = (i*a.r - rr*a.i)/m; return *this;}
      inline complex& operator ++ () {r += 1.0; return *this;}
      inline complex& operator -- () {r -= 1.0; return *this;}

      friend std::istream& operator >> (std::istream& is, complex& a);
      friend std::ostream& operator << (std::ostream& os,  const complex& a);

      friend inline complex operator + (const complex &a, const complex &b);
      friend inline complex operator + (const complex &a, double b);
      friend inline complex operator + (double b, const complex &a);
      friend inline complex operator - (const complex &a, const complex &b);
      friend inline complex operator - (const complex &a, double b);
      friend inline complex operator - (double b, const complex &a);
      friend inline complex operator * (const complex &a, const complex &b);
      friend inline complex operator * (const complex &a, double b);
      friend inline complex operator * (double b, const complex &a);
      friend inline complex operator / (const complex &a, const complex &b);
      friend inline complex operator / (const complex &a, double b);
      friend inline complex operator / (double b, const complex &a);
      friend inline int operator == (const complex &a, const complex &b);
      friend inline int operator != (const complex &a, const complex &b);
      friend inline complex conj(const complex & a);
      friend inline complex exp(const complex & a);
      friend inline complex log(const complex &a);
      friend inline complex cos(const complex & a);
      friend inline complex sin(const complex & a);

 /*  constructor from polar coordinates  */
      friend inline complex polar(double m,double ar);
};

const complex I(0.0,1.0);

inline double real(const complex & a){ return a.r;}
inline double imag(const complex & a){ return a.i;}
inline double norm(const complex & a){ return (a.r*a.r + a.i* a.i);}
inline double arg(const complex & a){ return atan2(a.i,a.r); }

inline complex operator + (const complex &a, const complex &b){
   return complex(a.r+b.r,a.i+b.i);
}
inline complex operator + (const complex &a, double b){
   return complex(a.r + b, a.i);
}
inline complex operator + (double b, const complex &a){
   return complex(a.r + b, a.i);
}

inline complex operator - (const complex &a, const complex &b){
   return complex(a.r-b.r,a.i-b.i);
}
inline complex operator - (const complex &a, double b){
   return complex(a.r - b, a.i);
}
inline complex operator - (double b, const complex &a){
   return complex(b -a.r, -a.i);
}

inline complex operator * (const complex &a, const complex &b){
   return complex(a.r*b.r-a.i*b.i,a.i*b.r + a.r*b.i);
}
inline complex operator * (const complex &a, double b){
   return complex(a.r* b, a.i*b);
}
inline complex operator * (double b, const complex &a){
   return complex(b *a.r, b*a.i);
}

inline complex operator / (const complex &a, const complex &b){
   double m = norm(b);
   return complex((a.r*b.r+a.i*b.i)/m,(a.i*b.r- a.r*b.i)/m);
}
inline complex operator / (const complex &a, double b){
   return complex(a.r/b, a.i/b);
}
inline complex operator / (double b, const complex &a){
   double m = norm(a);
   return complex(b *a.r/m, -b*a.i/m);
}

inline int operator == (const complex &a, const complex &b) {
   return((a.r ==b.r) && (a.i == b.i) ? 1 : 0 );}

inline int operator != (const complex &a, const complex &b) {
   return((a.r == b.r) && (a.i == b.i) ? 0 : 1);}

inline complex polar(double m, double argu){
   return  complex(m*cos(argu),m*sin(argu));}

inline complex conj(const complex & a) {
   return complex(a.r,-a.i);}

inline complex sqrt(const complex & a) {
   double m = sqrt(sqrt(norm(a)));
   return polar(m,0.5*arg(a));}

inline complex exp(const complex & a){
   double m = exp(a.r);
   return  polar(m,a.i);}

inline complex log(const complex &a){
   return complex( 0.5 * log(norm(a)), arg(a));}

inline complex pow(double a,const complex & b){
   return exp(log(a)*b);}

inline complex pow(const complex & a,double b){
   return exp(log(a)*b);}

inline complex pow(const complex & a,const complex & b){
   return exp(log(a)*b);}

inline complex cos(const complex & a){
   return complex( cos(a.r)*cosh(a.i), sin(a.r)*sinh(a.i)) ;}

inline complex sin(const complex & a){
   return complex( sin(a.r)*cosh(a.i), cos(a.r)*sinh(a.i)) ;}

inline double ABS(const complex & a){
   return sqrt(a.r*a.r + a.i* a.i);}


#endif    /*  COMPLEX_H  */

