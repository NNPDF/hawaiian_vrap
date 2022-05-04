
#ifndef NACLASSES_H
#define NACLASSES_H

/*   general purpose theory classes    */

#include "DClasses.h"
#include "vegas.h"
#include "pineappl_interface.h"

/*   numerical analysis classes    :    
         derive from these classes and overload the virtual functions     
                     to use the associated methods     */

class Curve{
   public:
      virtual double curve(double x)=0;

      double zero(double a, double b, double acc = 1.0e-6);
         /*   find a zero between a and b with accuracy acc  */
      double integral(double a, double b, double acc = 1.0e-5);
      double integral0(double a, double b, double acc = 1.0e-5);
         /*   computes the integral under the curve, from a to b, with 
            accuracy acc  ;  integral0 assumes that the function to be 
          integrated is zero on the boundary   */

      double trapint(double a, double b, int j);
      double trapint0(double a, double b, int j);
         /*  helper functions for the above */ 
};

double interpolate(DVector & Xa, DVector & Ya, double x, 
                double & dy, int nu=0);

/*  interpolates a data set (Xa[],YA[]) to a point x, returning the 
          corresponiding y, and returning an error estimate in dy.
       (If nu is defined, the program considers elements of Xa, Ya only
          up to  element nu.) */


class Surface{
   public:
      virtual double surface(DVector & V)=0;
      pinerap::CheffPanopoulos *vegas_piner;

      double vegas(VegasGrid & grid,int its, double & sd, double & chi2);
      double vegasint(VegasGrid & grid,int its, 
		              double & sd, double & chi2);
      double sampleint(VegasGrid & grid,int its, 
                              double & sd,  double & chi2);

          /*  Monte Carlo integration of the volume under 
                surface(V)  over a VegasGrid.
            The functions return the values of the integral, an estimate
          of the standard deviation (in sd) and a chi^2 per degree of 
           freedom (in chi2). 
         
         vegas accepts an initialized VegasGrid and adapts it
         vegasint integrates on an existing VegasGrid using stratified 
                       sampling
         sampleint integrates on an existing VegasGrid using random
                      (importance) sampling only
           (Both vegasint and sampleint change some parameters of the
              grid, so it might be advisable to create a copy of an 
               existing grid to feed to these functions.)  
	                                                      */

     DVector  sample(VegasGrid & grid, double  & weight);

     /*  sample chooses a random point from an existing VegasGrid, 
            returns this point, and returns the weight to be 
             assigned to the point as its value.     
                                                        */

};
     
class Flow{
   public:
      virtual DVector flow (double t, DVector & y)= 0;

      DVector stepintegrate(const DVector & y0, double t1, double t2,
                     DVector & tt, DMatrix & yt);
      DVector integrate(const DVector & y0, double t1, double t2,
               DVector & tt, DMatrix & yt, int & steps, double eps,
                 double dtsave, double hmin = 1.0e-8); 

        /* integrates the differential equation 
                  dy/dt = flow(t,y)
            from t1 to t2, returning the values of t in tt and the values
            of y in yt.  These objects should be initialized to 
               tt(0,nsteps)   and yt(nl,nh,0,nsteps)
           where nsteps is the number of steps for stepintegrate and the 
          maximum number of steps for integrate
   
       stepintegrate takes steps 4th-order Runge-Kutta steps
       integrate takes the number of 5th-order Runge-Kutta steps 
            necessary to achieve the accuracy eps and returns the 
              number of steps take in steps.  It saves results at
            intervals dtsave.  */

      /*   helper functions for the above:   */

      DVector rk4(const DVector & y, const DVector & dydt, double t, double h);
      DVector rkck(const DVector & y, const DVector & dydt, double t, 
                           double h,DVector & yerr);
      DVector rkqs(const DVector & y,const DVector & dydt,double & t, 
                     double hh, double eps, DVector & yscal, double & hdid, 
                      double & hnext);
   
};






#endif    /*  NACLASSES_H  */

