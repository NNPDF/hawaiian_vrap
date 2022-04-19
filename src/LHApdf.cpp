#include "LHApdf.h"
#include "LHAPDF/LHAPDF.h"
#include "pdf.h"
#include <iostream>

using namespace std;

// in pdf.C
extern int i_mrst, pdf_index, A_order;

void LHAComputePdf(double x,double Q,pdfArray& pa){
	LHAPDF::xfx(x,Q,&pa.tbar);
	pa.x=x;
	pa.u/=x;
	pa.ubar/=x;
	pa.d/=x;
	pa.dbar/=x;
	pa.s/=x;
	pa.sbar/=x;
	pa.c/=x;
	pa.cbar/=x;
	pa.b/=x;
	pa.bbar/=x;
	pa.t/=x;
	pa.tbar/=x;
	pa.gluon/=x;
}


void LHApdfInit(int pdft,int i_set){

	LHAPDF::initLHAPDF();
	
	pdf_index = pdft;
	if (pdft == cteq) { 
		LHAPDF::initPDFSetByName("cteq6m.LHpdf");	
	}
	
	else if (pdft == mrst) {
		switch (i_set){
			case 2: case 3: case 4: case 5: {
				LHAPDF::initPDFSetByName("MRST2001nlo.LHpdf");	
				LHAPDF::initPDF(i_set-1);
			} break;
			case 6: case 7: case 8: case 9: {
				LHAPDF::initPDFSetByName("MRST2002nnlo.LHpdf");	
				LHAPDF::initPDF(i_set-5);
			} break;
			case 15: {
				LHAPDF::initPDFSetByName("MRST2004nlo.LHpdf");	
				LHAPDF::initPDF(1);
			} break;
			case 16: {
				LHAPDF::initPDFSetByName("MRST2004nnlo.LHgrid");	
				LHAPDF::initPDF(1);
			} break;
			default: cerr << "Case not handeled in LHApdfInit"<< endl;
		}
		
	} // for MRST we only init. i here.


// for Alekhin, the second arument must be "order_flag" and we use it
// to pick LO,NLO,NNLO:
	else if (pdft == alekhin) { 
		A_order = i_set + 1; 
	} 


	}






