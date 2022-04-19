
#ifndef VEGAS_H
#define VEGAS_H

/*   a special class for the integration program vegas   */

#include "DClasses.h"


const double VEGASALPH = 1.5;
const int VNDMX = 50;

     
class VegasCounter;



class VegasGrid {

  public:

   int nl,nh,ndim,nd; double xnd;
       /*  vectors on the integration domain have j = nl .. nh ,
               ndim = nh-nl+1,
            the grid points per dimension are i= 1.. nd  (nd < VNDMX);
                   xnd = (double) nd                   */
    
   DVector xl, xh, dx;
        /*  lower and upper limits of the region of integrations, and the 
                width of this region   */
   DMatrix xi;
         /*  the grid itself:   
             xi[j][i] is the ith grid point in the jth dimension;
               xi[j][0] = 0 ,  ...   xi[j][nd] = 1    */

   int npg, ng; 
   double calls, xjac, dxg, dv2g;
        /* vegas makes calls calls per iteration x it iterations;
            the other three items control the precise number of calls:
     without stratified sampling, ng =calls, dxg = nd, xjac = volume/calls,
                                       dv2g = 1/(calls-1)
     with stratified sampling, ng is smaller, dxg = nd/ng , 
           npg = ncall/(ng)^ndim, giving  and no. of calls = npg* ng^ndim
                 xjac = volume/calls,    dv2g = 1/(npg-1)          */

   double si, swgt, schi;
   int it;
      /*   running estimates of the integral and its chi2:
             integral = si/swgt ; 
             sd = sqrt(1.0/swgt) ; 
             chi2 = (schi - si* si/swgt)/(it-1) 
       where it = no of iterations run                  */
        
    int mds,nprt;
        /*  controls:   mds = 0 disables stratified sampling
                        nprt = 0 outputs a complete record of parameters
                               1 outputs also the various intermediate grids
			       */

    VegasGrid(DVector xl, DVector xh, int ncalls, int ranseed);
    VegasGrid(const VegasGrid & A);

    VegasGrid & operator = (const VegasGrid & A);

  /*  this function reinitializes the grid, keeping the same integration
         region and using the number of calls ncall */

    void reset(int ncall); 

  /*   this function clears the accumulated answers from a given VegasGrid */

    void clear();



/* this function changes the number of calls associated with a VegasGrid */

    void morecalls(int nc);

     /*   this function returns a VegasGrid modified for random sampling */
   
    VegasGrid sampling(int ncalls)   const;

  /*   here are some helper functions that organize the operation of vegas */
    
    void reset();
    void refine(DMatrix d);
    void rebin(double rc, DVector & r, int j);
    DVector choose(double & xvol,VegasCounter & VC);
    DVector choose(double & xvol,VegasCounter & VC, IVector & ia);
    void record(double fb, double f2b, double & ti, double & tsi);
    void finish(double ti, double & tsi, double & integral, double & sd, 
                          double & chi2);
    
    void integralinfo(int it,double ti,double tsi,double integral,
                   double sd,double chi2);        

}; 

class VegasCounter{ 

     friend class VegasGrid;

     public:  
       int step, done;
       int nl,nh,npg,ng;
       int kct;
       IVector kg;
      /* set of counters determining the distribution of points in 
            stratified sampling: 
                kct  =  1, 2,  ... , npg  
                kg[j] = 1, 2, ..., ng  
  for each set of kg[j], when kct = npg, we set srecord = 1
  when all of the  kg[j] reach ng, we set done = 1 and the iteration is done */

   
       VegasCounter(const VegasGrid & grid);
       void reset();
       void resetstep();
       void move();


};
   

# endif    /*  VEGAS_H   */
