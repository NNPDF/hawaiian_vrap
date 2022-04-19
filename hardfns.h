/* =============================================================
   Special functions required for the "hard" NNLO contributions 
   to the DY rapidity distributions.
   ============================================================= */

#ifndef HARDFNS_H
#define HARDFNS_H

#include <iomanip>
#include "dilog.h"
#include "complex.h"

double real_conv(double y, double z); // Kirill -> my conventions
double DD1(double n, double y);  // "plus" distributions (w.o. subtraction)
double J3(double y, double z);
// Functions needed for J2:
double M1(double x, double xp, double xm);
double I_analytic(double y, double z);
double J27(double y, double z);
double J2(double y, double z);
// J1 is not independent, but here it is anyway:
double J1(double y, double z);
double J21(double y, double z);
// u -> 1 limits:
double J3_u1(double z);
double J27_u1(double z);
double J2_u1(double z);
double J21_u1(double z);
double J2rem(double z);
// For singularity at z = (2*u/(1+u))^2:
double y_sing(double z);
# endif    /*  HARDFNS_H  */
