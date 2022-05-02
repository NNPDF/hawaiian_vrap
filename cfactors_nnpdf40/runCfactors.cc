#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
#include <vector>
#include <cmath>

using namespace std;



void run(string expname) {

    ostringstream os;

    os << "output/" << expname << "_NLO_comparison.dat";
    //cout << os.str() << endl;
    ifstream infile1(os.str());
    os.str("");
    os << "output/" << expname << "_NNLO_comparison.dat";
    //cout << os.str() << endl;
    ifstream infile2(os.str());
    os.str("");

    vector<double> sigma1, sigma2;
    double c,d;

    while(infile1.good()) {
        if(expname=="E866R") {
            infile1 >> d >> d >> d >> c >> d >> d;
        }
        else {
            infile1 >> d >> d >> d >> c;
        }
        sigma1.push_back(c);
    }
    sigma1.pop_back();
    infile1.close();

    while(infile2.good()) {
        if(expname=="E866R") {
            infile2 >> d >> d >> d >> c >> d >> d;
        }
        else {
            infile2 >> d >> d >> d >> c;
        }
        sigma2.push_back(c);
    }
    sigma2.pop_back();
    infile2.close();

    os << "output/" << expname << "_Cfactors.dat";
    //cout << os.str() << endl;
    ofstream ofile(os.str());
    //
    for(int i=0; i<sigma2.size(); i++) {
        ofile << setw(10) << sigma2[i]/sigma1[i]
            << endl;
    }
    ofile.close();

}

int main(int argc, char *argv[]) {

    run("E866P");
    run("E866R");
    run("E605");

    return 0;
}
