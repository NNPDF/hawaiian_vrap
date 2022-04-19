
#include "random.h"

/*  ran1 returns a uniform random deviate between 0.0 and 1.0 using the 
         algorithm ran1 in Numerical Recipes, p. 280. 

    To initialize the random number sequence, call ran1(n), where n 
        is a NEGATIVE integer     */

double ran1(int n){
  static short init=1;
  const int IM = 2147483647;
  const double AM = 1.0/IM;
  const int IA = 16807;
  const int IQ = 127773;
  const int IR = 2836;
  const int NTAB = 32;
  const int NDIV = (1 + (IM-1)/NTAB);
  const double EPS = 1.2e-10;
  const double RNMX = 1.0 - EPS;

  static int idum = 1;
  static int iy;
  static long iv[NTAB];
  int j,k;

  if (n < 0) {
    init = 1;
    idum = -n;
  }
  if (init) {
    for (j = NTAB+7; j>=0; j--) {
      k = idum/IQ;
      idum = IA*(idum - k*IQ) - k* IR;
      if (idum<0) idum += IM;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
    init = 0;
  }
  k = idum/IQ;
  idum = IA*(idum - k*IQ) - k* IR;
  if (idum < 0) idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = idum;
  double r = AM*iy;
  if (r > RNMX) return RNMX;
  return r;
}
      
/*  ran2 returns a uniform random deviate between 0.0 and 1.0 using the 
         algorithm ran2 in Numerical Recipes, p. 282. 

    To initialize the random number sequence, call ran2(n), where n 
        is a NEGATIVE integer     */

double ran2(int n){
  static short init=1;
  const int IM1 = 2147483563;
  const int IM2 = 2147483399;
  const double AM = 1.0/IM1;
  const int IMM1 = IM1-1;
  const int IA1 = 40014;
  const int IA2 = 40692;
  const int IQ1 = 53668;
  const int IQ2 = 52774;
  const int IR1 = 12211;
  const int IR2 = 3791;
  const int NTAB = 32;
  const int NDIV = (1 + IMM1/NTAB);
  const double EPS = 1.2e-10;
  const double RNMX = 1.0 - EPS;

  static int idum = 1;
  static int idum2;
  static int iy;
  static long iv[NTAB];
  int j,k;

  if (n < 0) {
    init = 1;
    idum = -n;
  }
  if (init) {
    idum2 = idum;
    for (j = NTAB+7; j>=0; j--) {
      k = idum/IQ1;
      idum = IA1*(idum - k*IQ1) - k* IR1;
      if (idum<0) idum += IM1;
      if (j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
    init = 0;
  }
  k = idum/IQ1;
  idum = IA1*(idum - k*IQ1) - k* IR1;
  if (idum < 0) idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2*(idum - k*IQ2) - k* IR2;
  if (idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if (iy < 1) iy += IMM1;
  double r = AM*iy;
  if (r > RNMX) return RNMX;
  return r;
}
      




