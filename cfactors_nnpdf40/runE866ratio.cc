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

  ifstream infile("input/E866R.dat");
  vector<double> Y, M;
  double a,b;
  string dummy;
  
  while(infile.good()) {
    infile >> dummy >> dummy >> a >> b >> dummy;
    Y.push_back(a);
    M.push_back(b);
  }
  Y.pop_back();
  M.pop_back();
  infile.close();
  
  double sqrtS = 38.8;


  if(doNLO) {
    ostringstream osn;
    osn << "output/E866R_NLO_comparison.dat";
    ofstream ofile(osn.str());
    //  
    for(int i=0; i<Y.size(); i++) {
      double vrapP, vrapD;
      ostringstream os;

      cout << M[i] << "   " << Y[i] << endl;

      os << "Vrap inputE866nlo.dat " << M[i] << " " << Y[i] << " > output/outB.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outB.dat > output/outB2.dat");
      ifstream infile2("output/outB2.dat");
      infile2 >> vrapP;
      infile2.close();
      os.str("");
      os << "Vrap inputE866deutnlo.dat " << M[i] << " " << Y[i] << " > output/outB.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outB.dat > output/outB2.dat");
      infile2.open("output/outB2.dat");
      infile2 >> vrapD;
      infile2.close();
      ofile << setw(4) << i
	    << setw(10) << Y[i]
	    << setw(10) << M[i]
	    << setw(15) << vrapD/vrapP
	    << setw(15) << vrapP
	    << setw(15) << vrapD
	    << endl;
    }
    ofile.close();
  }


  if(doNNLO) {
    ostringstream osn;
    osn << "output/E866R_NNLO_comparison.dat";
    ofstream ofile(osn.str());
    //  
    for(int i=0; i<Y.size(); i++) {
      double vrapP, vrapD;
      ostringstream os;

      cout << M[i] << "   " << Y[i] << endl;

      os << "Vrap inputE866nnlo.dat " << M[i] << " " << Y[i] << " > output/outB.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outB.dat > output/outB2.dat");
      ifstream infile2("output/outB2.dat");
      infile2 >> vrapP;
      infile2.close();
      os.str("");
      os << "Vrap inputE866deutnnlo.dat " << M[i] << " " << Y[i] << " > output/outB.dat";
      system(os.str().c_str());
      system("awk 'END {print $NF}' output/outB.dat > output/outB2.dat");
      infile2.open("output/outB2.dat");
      infile2 >> vrapD;
      infile2.close();
      ofile << setw(4) << i
	    << setw(10) << Y[i]
	    << setw(10) << M[i]
	    << setw(10) << vrapD/vrapP
	    << setw(10) << vrapP
	    << setw(10) << vrapD
	    << endl;
    }
    ofile.close();
  }

  return 0;
}
