
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <cmath>

using namespace std;

int main(int argc, char *argv[]) {

  bool doNLO=false, doNNLO=false;
  if     (argc > 2) doNNLO = true;
  else if(argc > 1) doNLO = true;

  ifstream infile("input/E866P.dat");
  vector<double> Y, M, sigma;
  double a,b,c,d;
  
  while(infile.good()) {
    infile >> d >> a >> b >> c >> d;
    Y.push_back(a);
    M.push_back(b);
    sigma.push_back(c);
  }
  Y.pop_back();
  M.pop_back();
  sigma.pop_back();
  infile.close();
  
  double sqrtS = 38.8;

  if(doNLO) {
    ostringstream osn;
    osn << "output/E866P_NLO_comparison.dat";
    ofstream ofile(osn.str());
    //  
    for(int i=0; i<Y.size(); i++) {
      ostringstream os;

      cout << M[i] << "   " << Y[i] << endl;

      os << "Vrap inputE866nlo.dat " << M[i] << " " << Y[i] << " > output/outA.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outA.dat > output/outA2.dat");
      ifstream infile2("output/outA2.dat");
      double vrap;
      infile2 >> vrap;
      infile2.close();
      //vrap /= pow(M[i],3);
      //vrap *= pow(sqrtS,3);
      ofile << setw(4) << i
	    << setw(10) << Y[i]
	    << setw(10) << M[i]
	    << setw(10) << vrap
	    << endl;
    }
    ofile.close();
  }


  if(doNNLO) {
    ostringstream osn;
    osn << "output/E866P_NNLO_comparison.dat";
    ofstream ofile(osn.str());
    //  
    for(int i=0; i<Y.size(); i++) {
      ostringstream os;

      cout << M[i] << "   " << Y[i] << endl;

      os << "Vrap inputE866nnlo.dat " << M[i] << " " << Y[i] << " > output/outA.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outA.dat > output/outA2.dat");
      ifstream infile2("output/outA2.dat");
      double vrap;
      infile2 >> vrap;
      infile2.close();
      //vrap /= pow(M[i],3);
      //vrap *= pow(sqrtS,3);
      ofile << setw(4) << i
	    << setw(10) << Y[i]
	    << setw(10) << M[i]
	    << setw(10) << vrap
	    << endl;
    }
    ofile.close();
  }

  return 0;
}
