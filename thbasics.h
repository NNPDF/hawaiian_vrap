
#ifndef THBASICS_H
#define THBASICS_H

/*  header file for general purpose routines for theory classes,

       incorporating the utility routines from

     Numerical Recipes, by Press, Teukolsky, Vetterling, and Flannery

             C version,  2nd ed  1992

                     modified for C++ by MEP               */


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdarg.h>

const double PI = M_PI;
const double E = M_E;

const double TINY =  1.0e-30;


template<class T>
inline
T SQR(T a){return a*a;}

double ABS(double a); 
int  ABS(int a);

double MIN(double a,double b);
double  MAX(double a, double b);
double SIGN(double a, double b);

int MIN(int a, int b);
int  MAX(int a, int b);
int SIGN(int a, int b);

int INT(double a);

void therror(char error_text[]);
/*        standard error handler  */


#endif    /*  THBASICS_H  */

