
#include "thbasics.h"

using namespace std;

double ABS(double a){return (a >0 ? a: -a);}
int  ABS(int a) {return (a >0 ? a: -a);}

double MIN(double a,double b) {return (a < b ? a : b);}
double  MAX(double a, double b) {return  (a > b ? a: b);}
double SIGN(double a, double b){return (b < 0 ? -ABS(a): ABS(a));}

int MIN(int a, int b) {return (a < b ? a : b);}
int  MAX(int a, int b) {return  (a > b ? a: b);}
int SIGN(int a, int b){return (b < 0 ? -ABS(a): ABS(a));}

int INT(double a){return (int) a;}

void therror(const char error_text[]) {
/*        standard error handler  */
   cerr << "TH classes run-time error..." << endl ;
   cerr << "     " << error_text << endl;
   cerr << "   ...  now exiting to system ... " << endl;
   exit(1);
}


