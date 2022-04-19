#ifndef DILOG_H
#define DILOG_H

//#include "NAClasses.h"
//#include "DClasses.h"
#include "complex.h"

const double PI_DEF = 3.1415926535897932385;
const double PISQ   =  9.86960440108935861883;
const double ZETA2  =  1.64493406684822643647;
const double PISQ6  =  1.64493406684822643647;
const double PI4    = 97.409091034002437236440;
const double ZETA3  =  1.2020569031595942855;
const double ZETA4  =  1.082323233711138191516;
const double ZETA5  =  1.036927755143369926331;

double li2(double x);
double dilog(double x);
double ReLi2(complex x);
double ImLi2(complex x);
complex CLi2(complex x);
double li3(double x);
double S12(double x);  // only gets real part correct for x > 1.
double li4(double x);
double myli2(double x);
double mydilog(double x);
double i3_3m(double x, double y, double z);
double fastCl(double x);

# endif    /*  DILOG_H   */
