#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int       kords   = -1;
static int       kschems = -1;
static int       ksets   = -1;

#define nxb   99
#define nq    20
#define np     9
#define nvar  15

static double    f [nq + 1][nxb][np + 1]; 
static double    df[nq + 1][nxb][np + 1][nvar];

static double    dels = 0.0, delx = 0.0, x1 = 0.0, delx1 = 0.0, xlog1 = 0.0;
static int       nxbb = 0;

void a02(double xb, double q2, double pdfs[10], double dpdfs[15][10],
         int *npdf, int *npar, int kord, int kschem, int kset, int i0)
{
//
//     This is a code for the parton distributions with account 
//     of their experimental (stat+syst) and theoretical uncertainties. 
//     The Q**2 range is 1.4 < Q**2 < 2e8, the x range is 1e-7 < x < 1. 
//
//  Input parameters:
//        kord=1 -- the LO PDFs
//        kord=2 -- the NLO PDFs
//        kord=3 -- the NNLO PDFs
//      
//        kschem=0 -- the fixed-flavor-number (FFN) scheme 
//        kschem=1 -- the variable-flavor-number (VFN) scheme
//
//        kset=0 -- nominal PDFs
//        kset=1 -- PDFs with mass of c-quark increased from 1.5 to 1.75 GeV
//        kset=2 -- PDFs with the strange sea suppression factor increased from 
//                  0.42 to 0.52
//        kset=3 -- PDFs with the choice B (slow evolution) for the NNLO kernel 
//                  (used with kord=3 only)
//
//  Output parameters:
//     The array pdfs contains fitted values of the strong coupling constant 
//     and the parton distributions at given x and Q:
//        pdfs(0) -- \alpha_s
//        pdfs(1) -- valence u-quarks 
//        pdfs(2) -- valence d-quarks
//        pdfs(3) -- gluons 
//        pdfs(4) -- sea u-quarks 
//        pdfs(5) -- s-quarks 
//        pdfs(6) -- sea d-quarks 
//        pdfs(7) -- c-quarks
//        pdfs(8) -- b-quarks
//        pdfs(9) -- t-quarks
//     npdf is the number of PDFs returned (npdf=6 for the FFN PDFs and 9 for 
//     the VFN ones).
//     Output array dpdfs[ipar][ipdf] contains derivatives of \alpha_s and
//     the PDFs on the fitted parameters with the number of the parameters 
//     returned in npar. With the derivatives of \alpha_s included one can take 
//     into account the correlations of the fitted PDFs with \alpha_s as well.
//     All derivatives are transformed to the orthonormal 
//     basis of eigenvectors of the parameters error matrix. For this reason 
//     the variation of the PDFs in the derivatives directions can be performed 
//     independently. For example the dispersion of the i-th PDF can be stored 
//     in delpdf using the code 
//
//-----------------
//          double delpdf = 0.0;
//          for (int k = 0; k < npar; k++) {
//            delpdf += dpdfs[k][i] * dpdfs[k][i];
//          }
//-----------------
//     and its random value can be stored in rpdf using the code 
//-----------------
//          double rpdf = pdfs[i]          
//          for (int k = 0; k < npar; k++) {
//             double s = 0.;
//             for (int m = 0; m < 96; m++) {
//                s += (2.0 * drand48() - 1.0) / sqrt(32.0);
//             }
//             rpdf += s * dpdfs[k][i];
//          }
//-----------------
//          
//         Reference: Phys. Rev. D67, 14002 (2003) [hep-ph/0211096]
//      
//         Comments: alekhin@sirius.ihep.su                      
//                                                               
//     Initial version: Nov 2002      
//     Revision of May 2003: interpolation scheme was simplified without 
//                           loosing of the quality; value of \alpha_s(Q)
//                           is now stored in pdfs(0), the error in 
//                           \alpha_s(Q) and their correlation with the   
//                           PDFs are included in the derivatives stored in 
//                           dpdfs; the grid is expanded to the lower and
//                           higher values of Q.
//    Revision of June 2003: minor change in the calling sequence (xb and q2
//                           are now left unchanged even if they are out of the
//                           range allowed by the parameterization).   
//    Revision of Sep 2003:  the grid spacing is changed; parabolic 
//                           interpolation has been implemented; the grid was 
//                           expanded down to Q^2=0.8 Gev^2.
//    Revision of Nov 2003:  Translated into C (Willy Langeveld)
//                           Modified for speed in Vrap application;
//                           x=0.3 bug fixed (L.Dixon)
//
   const double xmin = 1e-7;
   const double xmax = 1.0;
   const double qsqmin = 0.8;
   const double qsqmax = 2e8;
   int n, m, i, k;
   double x, qsq, a, ss, b;
//
// put in your local address of the PDFs files in LOCDIR
// or set environment variable ALEKHIN_DIR [WL]
//
   const char *locdir = "/afs/slac.stanford.edu/u/th/lance/c/pdf/alekhin/";

   const char *pdford[3]   = { "1", "2", "3" };
   const char *pdfschem[2] = { "ffn", "vfn" };
   const char *pdfset[4]   = { "", "_mc", "_ss", "_kr" };
//
// Check for environment variable [WL]:
//
   char *p = getenv("ALEKHIN_DIR");
   if (p) locdir = p;

   if (kschem == 0) *npdf = 6;
   else             *npdf = 9;

   *npar = nvar;

   if ((kords != kord) || (kschems != kschem) || (ksets != kset)) {
      char file[512];
      FILE *fp;

      kords   = kord;
      kschems = kschem;
      ksets   = kset;    

      dels = (log(log(qsqmax / 0.04)) - log(log(qsqmin / 0.04))) / (double) (nq - 1);

      nxbb = nxb / 2;
      x1 = 0.3;
      xlog1 = log(x1);
      delx = (log(x1) - log(xmin)) / (double) (nxbb - 1);
      delx1 = (1 - x1) * (1 - x1) / (double) (nxbb + 1);

      sprintf(file, "%sa02.pdfs_%s_%s%s", locdir, pdford[kord - 1],
                                          pdfschem[kschem], pdfset[kset]);
      fp = fopen(file, "r");
      if (fp) {
         for (n = 0; n < nxb - 1; n++) {
            for (m = 0; m < nq; m++) {
               for (i = 0; i <= *npdf; i++) {
                  fscanf(fp, "%lg", &f[m][n][i]);
               }
            }
         }
         fclose(fp);
      }
      else {
         printf("The PDF set is unavailable (FILE: %s)\n", file);
         exit(1);
      }

      sprintf(file, "%sa02.dpdfs_%s_%s", locdir, pdford[kord - 1],
                                          pdfschem[kschem]);
      fp = fopen(file, "r");
      if (fp) {
         for (n = 0; n < nxb - 1; n++) {
            for (m = 0; m < nq; m++) {
               for (i = 0; i <= *npdf; i++) {
                  for (k = 0; k < *npar; k++) {
                     fscanf(fp, "%lg", &df[m][n][i][k]);
                  }
               }
            }
         }
         fclose(fp);
      }
      else {
         printf("The PDF set is unavailable (FILE: %s)\n", file);
         exit(1);
      }

      for (i = 1; i <= *npdf; i++) {
         for (m = 0; m < nq; m++) {
            f[m][nxb - 1][i] = 0.0;
            for (k = 0; k < *npar; k++) {
               df[m][nxb - 1][i][k] = 0.0;
            }
         }
      }

      for (m = 0; m < nq; m++) {
         f[m][nxb - 1][0] = f[m][nxb - 2][0];
         for (k = 0; k < *npar; k++) {
            df[m][nxb - 1][0][k] = df[m][nxb - 2][0][k];
         }
      }
   }

   if ((q2 < qsqmin) || (q2 > qsqmax)) {
      printf("  WARNING:  Q^2 VALUE IS OUT OF RANGE %g\n", q2);
   }
   if ((xb < xmin) || (xb > xmax)) {
      printf("  WARNING:   X  VALUE IS OUT OF RANGE %g\n", xb);
   }

   x = (xb > xmin) ? xb : xmin;
   x = (xb < xmax) ? xb : xmax;

   qsq = (q2 > qsqmin) ? q2 : qsqmin;
   qsq = (q2 < qsqmax) ? q2 : qsqmax;

   if (x > x1) {
      double xd = (1 - x1) * (1 - x1) - (1 - x) * (1 - x);
      n  = (int) (xd / delx1) + nxbb - 1;
      a  = xd / delx1 - n + nxbb - 1;
   }
   else {
      double xd = log(x) - xlog1;
      n  = (int) (xd / delx) + nxbb - 2;
      a  = xd / delx - n + nxbb - 1;
   }

   ss = log(log(qsq / 0.04)) - log(log(qsqmin / 0.04));
   m = (int) (ss/dels);
   b = ss / dels - m;

   for (i = 0; i <= *npdf; i++) {
     if (i==i0) {
      if ((n > 0) && (m > 0) && (n != 48)) { 
         pdfs[i] =  f[m  ][n  ][i] * (1 + a * b - a * a - b * b)
                  + f[m+1][n+1][i] * a * b
                  + f[m  ][n+1][i] * a * (a - 2 * b + 1) / 2.
                  + f[m+1][n  ][i] * b * (b - 2 * a + 1) / 2.
                  + f[m  ][n-1][i] * a * (a - 1) / 2.
                  + f[m-1][n  ][i] * b * (b - 1) / 2.;
      }
      else {
         pdfs[i] =  (1 - a) * (1 - b) * f[m][n  ][i] + (1 - a) * b * f[m+1][n  ][i]
                  +    a    * (1 - b) * f[m][n+1][i] +    a    * b * f[m+1][n+1][i];
      }
      /*  remove dpdfs functionality, trade for pdfs speed!
      for (k = 0; k < *npar; k++) {
         if ((n > 0) && (m > 0) && (n != 48)) {
            dpdfs[k][i] =  df[m  ][n  ][i][k] * (1 + a * b - a * a -b * b)
                         + df[m+1][n+1][i][k] * a * b
                         + df[m  ][n+1][i][k] * a * (a - 2 * b + 1) / 2.
                         + df[m+1][n  ][i][k] * b * (b - 2 * a + 1) / 2.
                         + df[m  ][n-1][i][k] * a * (a - 1) / 2.
                         + df[m-1][n  ][i][k] * b * (b - 1) / 2.;
         }
         else {
            dpdfs[k][i] =  (1 - a) * (1 - b) * df[m][n  ][i][k] + (1 - a) * b * df[m+1][n  ][i][k]
                         +    a    * (1 - b) * df[m][n+1][i][k] +    a    * b * df[m+1][n+1][i][k];
         }
      }
      */
     }
   }
   return;
}
