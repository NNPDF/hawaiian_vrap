
#ifndef RANDOM_H
#define RANDOM_H

/*   random number generators    */


double ran1(int i = 0);

/*  ran() returns a uniform random deviate between 0.0 and 1.0 using the 
         algorithm ran1 in Numerical Recipes, p. 280.
   
      Initialize the random number sequence by calling ran1(n), where
         n is a NEGATIVE integer.  If you forget, the sequence will 
         be initialized at n=(-1).  */


double ran2(int i = 0);

/*  ran2() returns a uniform random deviate between 0.0 and 1.0 using the 
         algorithm ran2 in Numerical Recipes, p. 282.
   
      Initialize the random number sequence by calling ran2(n), where
         n is a NEGATIVE integer.  If you forget, the sequence will 
         be initialized at n=(-1).  */

/*  according to Numerical Recipes, ran1 is 1.5 times faster than ran2
       but has a period of 2X 10^9.  The authors recommend ran2, with a 
        period of 2 X 10^18, whenever a single computation requires more
         than 10^8 random numbers. */

#endif    /*  RANDOM_H  */

