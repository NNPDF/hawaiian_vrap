/* Interface to CTEQ5 and MRST2001 PDFs and alpha_s.
  We use mrst99.C and mrst99.h as a kludge to access the 2001 data files.
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
   i=9:  fit jet data,       a_s(M_Z) = 0.1180, vnvalf1180j.dat [cor09.dat]
  DUMMY FILES (copied from pdf99):   cor10.dat, cor11.dat, cor12.dat
*/
#ifndef PDF_H
#define PDF_H
enum pdf_type { cteq = 1, mrst = 2, alekhin = 3 } ;
enum collider { pp = 1, ppbar = 2 } ;
// Controls number of factors of e_q^2 weighting the q-\bar{q} luminosity
// function:
enum process { gluon = 0, DY = 1, gamgam = 2 } ;
// initialize a pdf.
void pdf_init(pdf_type pdft, int i_set);
// call pdfs in a standard notation, after the above initialization:
double u(double x, double Q);
double ubar(double x, double Q);
double d(double x, double Q);
double dbar(double x, double Q);
double str(double x, double Q);
double chm(double x, double Q);
double bot(double x, double Q);
double g(double x, double Q);
// Alekhin writes an alpha_s output too:
double Alekhin_alpha_s(double Q);
# endif    /*  PDF_H  */
