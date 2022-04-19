
#include "complex.h"

std::istream & operator >> (std::istream& is, complex & a){
   return is >> a.r  >> a.i;   }
std::ostream & operator << (std::ostream & os, const complex & a){
   return os << "( "<< a.r <<" + " << a.i << " I )" ; }

