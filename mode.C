#include <string>
#include <iostream>
#include "LHAPDF/LHAPDF.h"


using namespace std;

extern int order_flag;
extern double alpha_s_Z;

// mode names have to be one word

void set_mode(const string& mode){
//
	if ( mode == "MRST_2001_LO"){
		order_flag = 0;  alpha_s_Z = 0.130;  
		LHAPDF::initPDFSetByName("MRST2001lo.LHgrid");	
		LHAPDF::initPDF(1); 
		return;
	}
	if ( mode == "MRST_2001_NLO"){
		order_flag = 1;  alpha_s_Z = 0.119015;  
		LHAPDF::initPDFSetByName("MRST2001nlo.LHgrid");	
		LHAPDF::initPDF(1); 
		return;
	}
	if ( mode == "MRST_2001_NNLO"){
		order_flag = 2;  alpha_s_Z = 0.1155;  
		LHAPDF::initPDFSetByName("MRST2001nnlo.LHgrid");	
		LHAPDF::initPDF(1); 
		return;
	}
//
	if ( mode =="MRST_2003c_NLO"){
		order_flag = 1;  alpha_s_Z = 0.1165;  
		LHAPDF::initPDFSetByName("MRST2003cnlo.LHgrid");	
		LHAPDF::initPDF(1);		
		return;
	}  
	if ( mode =="MRST_2003c_NNLO"){
		order_flag = 2;  alpha_s_Z = 0.1153;		
		LHAPDF::initPDFSetByName("MRST2003cnnlo.LHgrid");	
		LHAPDF::initPDF(1);		
		return;
	} 
// 
	if ( mode =="MRST_2004_NLO"){
		order_flag = 1;  alpha_s_Z = 0.120;   
		LHAPDF::initPDFSetByName("MRST2004nlo.LHgrid");	
		LHAPDF::initPDF(1);
		return;
	}
	if ( mode == "MRST_2004_NNLO"){
		order_flag = 2;  alpha_s_Z = 0.1167;  
		LHAPDF::initPDFSetByName("MRST2004nnlo.LHgrid");	
		LHAPDF::initPDF(1); 
		return;
	}
//
	if ( mode =="MRST_2006_NNLO"){
		order_flag = 2;  alpha_s_Z = 0.1191;		
		LHAPDF::initPDFSetByName("MRST2006nnlo.LHgrid");	
		LHAPDF::initPDF(0);		
		return;
	}
// 
	if ( mode =="MSTW_2008_LO"){
		order_flag = 0;  alpha_s_Z = 0.13939;		
		LHAPDF::initPDFSetByName("MSTW2008lo68cl.LHgrid");	
		LHAPDF::initPDF(0);		
		return;
	}  
	if ( mode =="MSTW_2008_NLO"){
		order_flag = 1;  alpha_s_Z = 0.12018;		
		LHAPDF::initPDFSetByName("MSTW2008nlo68cl.LHgrid");	
		LHAPDF::initPDF(0);		
		return;
	}  
	if ( mode =="MSTW_2008_NNLO"){
		order_flag = 2;  alpha_s_Z = 0.11707;		
		LHAPDF::initPDFSetByName("MSTW2008nnlo68cl.LHgrid");	
		LHAPDF::initPDF(0);		
		return;
	}  
	if ( mode =="ABKM_2009_NLO"){
	  order_flag = 1;  alpha_s_Z = 0.1179;   // from Sven Moch	
		LHAPDF::initPDFSetByName("abkm09_5_nlo.LHgrid");
		LHAPDF::initPDF(0);
		//      double mZ = 91.1876;
		//     alpha_s_Z = LHAPDF::alphasPDF(mZ); Fix later	
		return;
	}  
	if ( mode =="ABKM_2009_NNLO"){
	  order_flag = 2;  alpha_s_Z = 0.1135; // from Sven Moch
		LHAPDF::initPDFSetByName("abkm09_5_nnlo.LHgrid");	
		LHAPDF::initPDF(0);
		// double mZ = 91.1876;
		// alpha_s_Z = LHAPDF::alphasPDF(mZ); Fix later	
		return;
	}  
	cerr << "Mode \'" << mode << "\' not known." << endl;

}

void set_mode(const string& filename,int iset){
		LHAPDF::initPDFSetByName(filename);	
		LHAPDF::initPDF(iset);		
		return;

	
}
