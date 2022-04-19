#include "pdf.h"
#include "a02.h"
#include "thbasics.h"
int i_mrst, pdf_index, A_order;
/* Interface to CTEQ5 and MRST2001 and Alekhin02 PDFs and alpha_s.
  New MRST interpolation formula, courtesy of Jeppe Andersen.
  We use mrst.cc and mrst.h to access the 2001 data files.
  MRST2001 indexing:
  LO:    (hep-ph/0201127)
   i=1:                       a_s(M_Z) = 0.130,  lo2002.dat [ -> cor01.dat ]
  NLO:   (hep-ph/0110215)
   i=2:  central gluon, a_s,  a_s(M_Z) = 0.119,  alf119.dat [ -> cor02.dat ]
   i=3:  lower a_s,           a_s(M_Z) = 0.117,  alf117.dat [ -> cor03.dat ]
   i=4:  higher a_s,          a_s(M_Z) = 0.121,  alf121.dat [ -> cor04.dat ]
   i=5:  fit jet data,        a_s(M_Z) = 0.121,  j121.dat   [ -> cor05.dat ]
  NNLO:  (hep-ph/0201127)
   i=6:  `average' evolution, a_s(M_Z) = 0.1155, vnvalf1155.dat [cor06.dat]
   i=7:  `fast' evolution,   a_s(M_Z) = 0.1155, vnvalf1155a.dat [cor07.dat]
   i=8:  `slow' evolution,   a_s(M_Z) = 0.1155, vnvalf1155b.dat [cor08.dat]
   i=9:  "mode 4", fit jet data,  a_s(M_Z) = 0.1180, vnvalf1180j.dat [cor09.dat]
  OLD NLO MRST99 FILES (copied from pdf99):  
   i=10:  MRST99 set 1, central gluon, a_s(M_Z) = 0.1175, COR01 [ -> cor10.dat]
   i=11:  MRST99 set 2, higher gluon, a_s(M_Z) = 0.1175, COR02 [ -> cor11.dat]
   i=12:  MRST99 set 3, lower gluon, a_s(M_Z) = 0.1175, COR03 [ -> cor12.dat]
  We checked our implementation of these pdfs using hep-ph/0201127,
  figure 3 for g(x): LO, NLO (average) and NNLO (average, fast and slow), and
  figure 4 for g(x), u(x) and d(x) [NNLO(average)/NLO(average) ratio].
// Added 1/24/2005:
  NLO AND NNLO MRST2003c ("conservative") FILES (see hep-ph/0308087):
    i=13:  NLO MRST2003c, a_s(M_Z) = 0.1165, [ -> mrst2003cnlo.dat]
    i=14:  NNLO MRST2003c, a_s(M_Z) = 0.1153, [ -> mrst2003cnnlo.dat]
  NLO AND NNLO MRST2004 ("conservative") FILES (see hep-ph/0410230):
    i=15:  NLO MRST2004, a_s(M_Z) = 0.120, [ -> mrst2004nlo.dat]
    i=16:  NNLO MRST2004, a_s(M_Z) = 0.1167, [ -> mrst2004nnlo.dat]
//
 Alekhin02 indexing:
  i=1:  NNLO (kord=3), variable-flavor (kschem=1), nominal (kset=0) PDFs.
*/

// arrays, etc. for Alekhin to write to:
double pdfs[10], dpdfs[15][10], delpdf[10];
int npdf, npar, n;

#if USE_OLD

// New initialization of MRST distributions 
// (may have to change if we want more):
class c_mrst parton1("cor01.dat");   // MRST2001 LO
class c_mrst parton2("cor02.dat");   // MRST2001 NLO
class c_mrst parton5("cor05.dat");   // MRST2001 NNO mode 4, fit jet data
class c_mrst parton6("cor06.dat");   // MRST2001 NNLO (default mode 1)
class c_mrst parton7("cor07.dat");   // MRST2001 NNLO `fast' evolution
class c_mrst parton8("cor08.dat");   // MRST2001 NNLO `slow' evolution
class c_mrst parton9("cor09.dat");   // MRST2001 NNLO mode 4, fit jet data
class c_mrst parton10("cor10.dat");   // MRST99 NLO set 1
class c_mrst parton11("cor11.dat");   // MRST99 NLO set 2
class c_mrst parton12("cor12.dat");   // MRST99 NLO set 3
class c_mrst parton13("mrst2003cnlo.dat");  // MRST2003c NLO
class c_mrst parton14("mrst2003cnnlo.dat"); // MRST2003c NNLO
class c_mrst parton15("mrst2004nlo.dat");  // MRST2004 NLO
class c_mrst parton16("mrst2004nnlo.dat"); // MRST2004 NNLO

// pointers to MRST distributions (couldn't get this to work???):
//class c_mrst * pdfchoice[2];
//pdfchoice[0] = &parton1;
//pdfchoice[1] = &parton2;


// Generic pdf initialization:

void pdf_init(pdf_type pdft, int i_set) {
  pdf_index = pdft;
  if (pdft == cteq) { pdfinit(i_set); }
  else if (pdft == mrst) { i_mrst = i_set; } // for MRST we only init. i here.
// for Alekhin, the second arument must be "order_flag" and we use it
// to pick LO,NLO,NNLO:
  else if (pdft == alekhin) { A_order = i_set + 1; } 
}

// call pdfs in a standard notation, after the above initialization:

double u(double x, double Q) {
  if (pdf_index==cteq) { return pdf(1,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 1);
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 4);
    return (pdfs[1] + pdfs[4])/x;
  }
  else { if (i_mrst==1) { return ( parton1.mrstparton(1,x,Q)
                                 + parton1.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==2) { return ( parton2.mrstparton(1,x,Q)
                                      + parton2.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==5) { return ( parton5.mrstparton(1,x,Q)
                                      + parton5.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==6) { return ( parton6.mrstparton(1,x,Q)
                                      + parton6.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==7) { return ( parton7.mrstparton(1,x,Q)
                                      + parton7.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==8) { return ( parton8.mrstparton(1,x,Q)
                                      + parton8.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==9) { return ( parton9.mrstparton(1,x,Q)
                                      + parton9.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==10) { return ( parton10.mrstparton(1,x,Q)
                                       + parton10.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==11) { return ( parton11.mrstparton(1,x,Q)
                                       + parton11.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==12) { return ( parton12.mrstparton(1,x,Q)
                                       + parton12.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==13) { return ( parton13.mrstparton(1,x,Q)
                                       + parton13.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==14) { return ( parton14.mrstparton(1,x,Q)
                                       + parton14.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==15) { return ( parton15.mrstparton(1,x,Q)
                                       + parton15.mrstparton(4,x,Q) )/x; }
         else if (i_mrst==16) { return ( parton16.mrstparton(1,x,Q)
                                       + parton16.mrstparton(4,x,Q) )/x; }
         else return 0.;
  }
}

double ubar(double x, double Q) {
  if (pdf_index==cteq) { return pdf(-1,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 4);
    return pdfs[4]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(4,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(4,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(4,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(4,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(4,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(4,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(4,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(4,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(4,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(4,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(4,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(4,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(4,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(4,x,Q)/x; 
         else return 0.;
  }
}

double d(double x, double Q) {
  if (pdf_index==cteq) { return pdf(2,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 2);
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 6);
    return (pdfs[2] + pdfs[6])/x; 
  }
  else { if (i_mrst==1) { return ( parton1.mrstparton(2,x,Q)
                                 + parton1.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==2) { return ( parton2.mrstparton(2,x,Q)
                                      + parton2.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==5) { return ( parton5.mrstparton(2,x,Q)
                                      + parton5.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==6) { return ( parton6.mrstparton(2,x,Q)
                                      + parton6.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==7) { return ( parton7.mrstparton(2,x,Q)
                                      + parton7.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==8) { return ( parton8.mrstparton(2,x,Q)
                                      + parton8.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==9) { return ( parton9.mrstparton(2,x,Q)
                                      + parton9.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==10) { return ( parton10.mrstparton(2,x,Q)
                                       + parton10.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==11) { return ( parton11.mrstparton(2,x,Q)
                                       + parton11.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==12) { return ( parton12.mrstparton(2,x,Q)
                                       + parton12.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==13) { return ( parton13.mrstparton(2,x,Q)
                                       + parton13.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==14) { return ( parton14.mrstparton(2,x,Q)
                                       + parton14.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==15) { return ( parton15.mrstparton(2,x,Q)
                                       + parton15.mrstparton(8,x,Q) )/x; }
         else if (i_mrst==16) { return ( parton16.mrstparton(2,x,Q)
                                       + parton16.mrstparton(8,x,Q) )/x; }
         else return 0.;
  }
}

double dbar(double x, double Q) {
  if (pdf_index==cteq) { return pdf(-2,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 6);
    return pdfs[6]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(8,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(8,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(8,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(8,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(8,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(8,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(8,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(8,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(8,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(8,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(8,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(8,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(8,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(8,x,Q)/x; 
         else return 0.;
  }
}

double str(double x, double Q) {
  if (pdf_index==cteq) { return pdf(3,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 5);
    return pdfs[5]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(6,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(6,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(6,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(6,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(6,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(6,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(6,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(6,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(6,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(6,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(6,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(6,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(6,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(6,x,Q)/x; 
         else return 0.;
  }
}

double chm(double x, double Q) {
  if (pdf_index==cteq) { return pdf(4,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 7);
    return pdfs[7]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(5,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(5,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(5,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(5,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(5,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(5,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(5,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(5,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(5,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(5,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(5,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(5,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(5,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(5,x,Q)/x; 
         else return 0.;
  }
}

double bot(double x, double Q) {
  // ONLY ZERO OUT b(x) IF YOU WANT TO e.g. COMPARE WITH MRST01 sigma_W,Z:
  // return 0.;
  if (pdf_index==cteq) { return pdf(5,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 8);
    return pdfs[8]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(7,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(7,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(7,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(7,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(7,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(7,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(7,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(7,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(7,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(7,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(7,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(7,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(7,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(7,x,Q)/x; 
         else return 0.;
  }
}

double g(double x, double Q) {
  if (pdf_index==cteq) { return pdf(0,x,Q); }
  if (pdf_index==alekhin) {
    a02(x, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 3);
    return pdfs[3]/x;
  }
  else { if (i_mrst==1) return parton1.mrstparton(3,x,Q)/x; 
         else if (i_mrst==2) return parton2.mrstparton(3,x,Q)/x; 
         else if (i_mrst==5) return parton5.mrstparton(3,x,Q)/x; 
         else if (i_mrst==6) return parton6.mrstparton(3,x,Q)/x; 
         else if (i_mrst==7) return parton7.mrstparton(3,x,Q)/x; 
         else if (i_mrst==8) return parton8.mrstparton(3,x,Q)/x; 
         else if (i_mrst==9) return parton9.mrstparton(3,x,Q)/x; 
         else if (i_mrst==10) return parton10.mrstparton(3,x,Q)/x; 
         else if (i_mrst==11) return parton11.mrstparton(3,x,Q)/x; 
         else if (i_mrst==12) return parton12.mrstparton(3,x,Q)/x; 
         else if (i_mrst==13) return parton13.mrstparton(3,x,Q)/x; 
         else if (i_mrst==14) return parton14.mrstparton(3,x,Q)/x; 
         else if (i_mrst==15) return parton15.mrstparton(3,x,Q)/x; 
         else if (i_mrst==16) return parton16.mrstparton(3,x,Q)/x; 
         else return 0.;
  }
}


#endif 
// Alekhin writes an alpha_s output too:

double Alekhin_alpha_s(double Q){
  a02(0.1, Q*Q, pdfs, dpdfs, &npdf, &npar, A_order, 1, 0, 0);
  return pdfs[0];
}
