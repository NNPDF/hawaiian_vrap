#ifndef _H_LHAPDF_H
#define _H_LHAPDF_H

struct pdfArray {
		double tbar;
		double bbar;
		double cbar;
		double sbar;
		double ubar;
		double dbar;
		double gluon;
		double d;
		double u;
		double s;
		double c;
		double b;
		double t;
		double x;
};


void LHApdfInit(int pdft,int iset);

void LHAComputePdf(double x,double Q,pdfArray& pa);

#endif /* _H_LHAPDF_H */




