

using namespace std;

#include "DClasses.h"

const int OFFSET = 1;
  
  /*  the indexing and storage allocation for vectors and matrices follow the 
        prescriptions of Numerical Recipes in C, by 
            Press, Teukolsky, Vetterling, and Flannery  */

DVector::DVector() { 
    therror(" to construct a DVector, use the constructor with parameters");
}

DVector::DVector(int xnl, int xnh){
       init(xnl,xnh);
       zero();
}

DVector::DVector(const DVector & A){
       init(A.nl, A.nh);
       load(A);
}

DVector::DVector(const IVector & A){
       init(A.nl, A.nh);
       loadI(A);
}

DVector::~DVector(){ 
	destroy();
}

DVector & DVector::operator = (const DVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible DVector");
        }
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = A.ptr[i];
	return * this;
}


int DVector::l() const {return nl;}
int DVector::h() const {return nh;} 

void DVector::fill(int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++) ptr[i] = va_arg(a,double);
        va_end(a);
}

void DVector::zero() { 
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = 0;
}


DVector DVector::operator - () const {
        DVector A(nl,nh);
        int i;
        for (i = nl; i<= nh; i++) A.ptr[i] = -ptr[i]; 
        return A;
}

DVector & DVector::operator *= (double  A) {
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] * A;
	return * this;
}

DVector & DVector::operator /= (double  A) {
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] / A;
	return *this;
}

DVector & DVector::operator += (const DVector & A) {
        if (A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible DVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] +  A.ptr[i];
	return * this;
}
 
DVector & DVector::operator -= (const DVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible DVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] - A.ptr[i];
	return * this;
}

 
void DVector::init(int xnl,int xnh){
      nl = xnl; nh = xnh;
      if ((nh-nl)<0) 
         therror(" the number of components of a DVector should be positive ");
      ptr = new double[nh-nl+1 + OFFSET];
      if (!ptr) therror("allocation failure in DVector()");
      ptr -= (nl - OFFSET); 
}

void DVector::load(const DVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void DVector::loadI(const IVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void DVector::destroy() {
	ptr += (nl-OFFSET);
	delete [] ptr;
}	



std::ostream & operator << (std::ostream & os, const DVector & A){
    int i;
    os << "<";
    for (i=A.nl;i<= A.nh;i++) os << "\t" << A.ptr[i];
    os << ">" << endl;
    return os;
}


int operator ==  (const DVector &A,const DVector &B) {
  if ( A.nl != B.nl  || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i;
  for (i=A.nl; i<= A.nh; i++) {  if (A.ptr[i] != B.ptr[i]) return 0;};
  return 1;
}


int operator !=  (const DVector &A,const DVector &B) {
  return (A == B ? 0 : 1);
}

int eqacc(const DVector & A, const DVector &B, double acc) {
  if ( A.nl != B.nl  || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i;
  for (i=A.nl; i<= A.nh; i++) {
    if (ABS(A.ptr[i] - B.ptr[i]) > acc) return 0; }; 
  return 1;
}

DVector operator + (const DVector & A, const DVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to add incompatible DVectors ");}
  DVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] + B.ptr[i];
  return C; 
}
    
DVector operator - (const DVector & A, const DVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to subtract incompatible DVectors ");}
  DVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] - B.ptr[i];
  return C; 
}

DVector operator * (double A, const DVector & B){
  DVector C(B.nl,B.nh);
  int i;
  for (i=B.nl; i<= B.nh; i++) C.ptr[i] = A * B.ptr[i];
  return C; 
}

DVector operator * (const DVector & A, double  B){
  DVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] * B;
  return C; 
}

DVector  operator / (const DVector & A, double  B){
  DVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] / B;
  return C; 
}

double  operator * (const DVector & A, const DVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to dot incompatible DVectors ");}
  double s = 0;
  int i;
  for (i=A.nl; i<= A.nh;i++) s += A.ptr[i] * B.ptr[i];
  return s;
}



DMatrix::DMatrix() { 
      therror("to construct a DMatrix, use the constructor with parameters");
}

DMatrix::DMatrix(int xml, int xmh, int xnl, int xnh){
       init(xml,xmh,xnl,xnh);
       zero();
}

DMatrix::DMatrix(const DMatrix & A){
       init(A.ml,A.mh,A.nl, A.nh);
       load(A);
}

DMatrix::DMatrix(const IMatrix & A){
       init(A.ml,A.mh,A.nl, A.nh);
       loadI(A);
}

DMatrix::DMatrix(int xnl, int xnh, double a) {
       init(xnl,xnh,xnl,xnh);
       int i,j;
       for (i=ml; i<= mh; i++) {
	 for (j=nl;j<=nh; j++) ptr[i][j] = 0;
	 ptr[i][i] =  a;
       }
}

DMatrix::~DMatrix(){ 
	destroy();
}

DMatrix &  DMatrix::operator = (const DMatrix & A) {
	if (A.ml != ml || A.mh != mh ||A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible DMatrix");
        }
	int i,j;
	for (i = ml; i<=mh; i++) 
           for (j = nl; j<=nh; j++)  ptr[i][j] = A.ptr[i][j];
	return * this;
}

int DMatrix::rl() const {return ml;}
int DMatrix::rh() const {return mh;}
int DMatrix::cl() const {return nl;}
int DMatrix::ch() const {return nh;} 

void DMatrix::fillRow(int m, int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++) ptr[m][i] = va_arg(a,double);
        va_end(a);
}

void DMatrix::readRow(int n,const DVector & A){
        if (A.nl != nl || A.nh != nh)
         therror(" attempt to fill a DMatrix with an incompatible DVector");
        int i;
        for (i= nl; i<= nh; i++)  ptr[n][i] = A.ptr[i];
}

void DMatrix::readColumn(int n, const DVector & A){
        if (A.nl != ml || A.nh != mh)
         therror(" attempt to fill a DMatrix with an incompatible DVector");
        int i;
        for (i= ml; i<= mh; i++)  ptr[i][n] = A.ptr[i];
}

void DMatrix::zero() {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = 0;
}

DVector DMatrix::row(int i) const{
	DVector A(nl,nh);
	int j;
        for (j=nl;j<=nh;j++) A[j] = ptr[i][j];
	return A;
}

DVector DMatrix::column(int j) const{
	DVector A(ml,mh);
	int i;
        for (i=ml;i<=mh;i++) A[i] = ptr[i][j];
	return A;
}

DMatrix DMatrix::operator - () const {
        DMatrix A(ml,mh,nl,nh);
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) A.ptr[i][j] = - ptr[i][j];
	return A;

}

DMatrix & DMatrix::operator *= (double A) {
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) ptr[i][j] = ptr[i][j] * A;
	return * this; 
}

DMatrix & DMatrix::operator /= (double A) {
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) ptr[i][j] = ptr[i][j] / A;
	return *this;
}

DMatrix & DMatrix::operator += (const DMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible DMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] +  A.ptr[i][j];
	return * this;
}

DMatrix & DMatrix::operator -= (const DMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible DMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] -  A.ptr[i][j];
	return * this;
}

void DMatrix::init(int xml, int xmh, int xnl,int xnh){
	ml = xml; mh = xmh; nl = xnl; nh = xnh;
        if ((mh-ml)<0) 
           therror(" the number of rows of a DMatrix should be positive ");
        if ((nh-nl)<0) 
           therror(" the number of columns of a DMatrix should be positive ");
	int nrow = mh-ml+1;
        int ncol = nh-nl+1;
	ptr = new double * [nrow + OFFSET];
	if (!ptr)  therror ("allocation failure in DMatrix()");
	ptr += OFFSET;
        ptr -= ml;    /*  now ptr has the value s.t. ptr + ml is the correct
                              starting address   */
	ptr[ml] = new double[nrow*ncol + OFFSET];
	if (!ptr[ml]) therror("allocation failure, stage 2, in DMatrix()");
	ptr[ml] += OFFSET;
	ptr[ml] -= nl;
             /* now evaluate the pointers to rows  */
	for (int i = ml+1; i<= mh; i++) ptr[i] = ptr[i-1] + ncol;
} 

void DMatrix::load(const DMatrix & A) {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = A.ptr[i][j];
}

void DMatrix::loadI(const IMatrix & A) {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = A.ptr[i][j];
}

void DMatrix::destroy() {
	double * w = ptr[ml] + nl - OFFSET;
	delete [] w;
	ptr += (ml - OFFSET);
	delete [] ptr;
}

ostream & operator << (ostream & os, const DMatrix & A){
    int i,j;
    for (i=A.ml;i<=A.mh;i++) {
       os << endl <<  "<";
       for (j=A.nl;j<= A.nh;j++) os << "\t" << A.ptr[i][j];
    os << "\t>" << endl;
    }
    os << endl;
    return os;
}

int operator ==  (const DMatrix &A,const DMatrix &B) {
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) if (A.ptr[i][j] != B.ptr[i][j]) return 0;
  return 1;
}

int operator !=  (const DMatrix &A,const DMatrix &B) {
  return (A == B ? 0 : 1);
}


int eqacc(const DMatrix &A,const DMatrix &B,double acc) {
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) 
       if (ABS(A.ptr[i][j] -B.ptr[i][j])> acc) return 0;
  return 1;
}

DMatrix operator + (const DMatrix & A, const DMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to add incompatible DMatrices ");
  DMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j]+ B.ptr[i][j];
  return C; 
}
    
DMatrix operator - (const DMatrix & A, const DMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to subtract incompatible DMatrices ");
  DMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] - B.ptr[i][j];
  return C;
} 

DMatrix operator * (double A, const DMatrix & B){
  DMatrix C(B.ml,B.mh,B.nl,B.nh);
  int i,j;
  for (i=B.ml; i<= B.mh; i++)  
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A * B.ptr[i][j];
  return C; 
}

DMatrix operator * (const DMatrix & A, double B){
  DMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] * B;
  return C; 
}

DMatrix operator / (const DMatrix & A, double B){
  DMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] / B;
  return C; 
}

DVector operator * (const DVector & A, const DMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible DVector and DMatrix ");
  DVector C(B.nl,B.nh);
  int i,j;
  for (j=B.nl; j<= B.nh; j++) {
    C.ptr[j] = 0;
    for (i=B.ml;i<=B.mh; i++) C.ptr[j] += A.ptr[i] * B.ptr[i][j];
  }
  return C; 
}

DVector operator * (const DMatrix & A, const DVector & B){
  if (A.nl != B.nl || A.nh != B.nh)
    therror(" attempt to multiply incompatible DMatrix and DVector ");
  DVector C(A.ml,A.mh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++) {
    C.ptr[i] = 0;
    for (j=B.nl;j<=B.nh; j++) C.ptr[i] += A.ptr[i][j] * B.ptr[j];
  }
  return C; 
}

DMatrix operator * (const DMatrix & A, const DMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible DMatrices ");
  DMatrix C(A.ml,A.mh,B.nl,B.nh);
  int i,j,k;
  for (i=A.ml; i<= A.mh; i++) {
    for (j=B.nl; j<= B.nh; j++) {
      C.ptr[i][j] = 0;
      for (k=A.nl;k<=A.nh;k++) C.ptr[i][j] += A.ptr[i][k] * B.ptr[k][j];
    }
  }
  return C; 
}

DMatrix dyad(const DVector & A, const DVector & B){
  DMatrix  C(A.nl,A.nh,B.nl,B.nh);
  int i,j;
  for (i=A.nl; i<= A.nh; i++)
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A.ptr[i] * B.ptr[j];
  return C; 
}


/*  methods for inverse,  determinant, etc.     */

DMatrix DMatrix::transpose() const {
  DMatrix A(nl,nh,ml,mh);
  int i,j;
  for (i=ml; i<= mh; i++)
    for (j=nl; j<= nh; j++) A.ptr[j][i] = ptr[i][j];
  return A;
}

double DMatrix::trace() const {
  if (mh != nh || ml != nl)
         therror(" attempt to take the trace of an asymmetric DMatrix ");
  int i;   double tr = 0;
  for (i = nl; i <= nh; i++)  tr += ptr[i][i];
  return tr;
}

void DMatrix::rowint(int i, int j) {
  int k;
  double el;
  for (k= nl; k<= nh; k++) {
    el = ptr[i][k];
    ptr[i][k] = ptr[j][k];
    ptr[j][k] = el;
  }
}

void DMatrix::colint(int i, int j) {
  int k;
  double el;
  for (k= ml; k<= mh; k++) {
    el = ptr[k][i];
    ptr[k][i] = ptr[k][j];
    ptr[k][j] = el;
  }
}

double DMatrix::maxRow(int i) const { 
        /* returns the largest ABS element in row i  */
  int k;  double big, temp;
  big = 0.0; 
  for (k= nl; k<= nh; k++) {
    temp = ABS(ptr[i][k]);
    if (temp > big )  big = temp;
  }
  if (big == 0.0) therror(" no result found by maxRow ");
  return big;
}   

double DMatrix::maxCol(int i) const { 
        /* returns the largest ABS element in column i  */
  int k;  double  big, temp;
  big = 0.0; 
  for (k= ml; k<= mh; k++) {
    temp = ABS(ptr[k][i]);
    if (temp > big )  big = temp;
  }
  if (big == 0.0) therror(" no result found by maxCol ");
  return big;
}   

/*   the following matrix inversion and determinant routines are based on
         LU decomposition, as described in Numerical Recipes,
            by Press, Teukolsky, Vetterling, and Flannery       */


DMatrix DMatrix::LUdecomp(int *  index, double & d) const {
  if (mh != nh || ml != nl)
         therror(" attempt to invert an asymmetric DMatrix ");
  int i,j,k, imax;
  double sum, t;
  double big,temp;
  DMatrix C = *this;   /* copy the matrix into C, and work on C */
  d = 1;                   /* d = (-1)^permutation     */
  DVector V(nl,nh);
  for (i = nl; i<= nh; i++)  V[i] = 1.0/maxRow(i);
  for (j = nl;j<= nh; j++) { /* this is the main loop for Crout's algorithm */
    for (i=nl; i< j; i++) {
      sum = C[i][j];
      for (k=nl;k<i;k++) sum -= C[i][k] * C[k][j];
      C[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i<=nh; i++) {
      sum = C[i][j];
      for (k=nl;k<j;k++) sum -= C[i][k] * C[k][j];
      C[i][j] = sum;
      temp = V[i]* ABS(sum);
      if (temp > big) {
	big = temp;
	imax = i;
      }
    }	
    if (j != imax) {
      C.rowint(j,imax);
      d = -d;
      V[imax] = V[j];
    }
    index[j-nl] = imax;
    t = C[j][j];
    if (t == 0.0) t = TINY;   /* recommended cheat in Numerical Recipes */
    if (j != nh){
      for (i = j+1; i<= nh; i++) C[i][j] /= t;
    }
  }
  return C;
}  
 
void DMatrix::LUbacksub(DMatrix & B,int * index) {
  /*   note that B is destroyed and replaced in this operation;
           on the other hand, the inverted matrix is preserved    */
  if (nh != B.mh || nl != B.ml)
        therror(" attempt to multiply incompatible DMatrices in LUbacksub ");
  int i,j,k, ii, ip;
  double sum;
  for (i= ml; i< mh; i++) {
    ip = index[i-ml];
    B.rowint(i,ip);
  }
  for (k = B.nl; k<= B.nh; k++) {
    ii = 0;                    /* ii makes this algorithm more efficient */
    for (i=ml; i<= mh; i++) {  /*  if B contains many zeros              */
      sum = B[i][k];
      if (ii)
	for (j = ii; j< i; j++)  sum -= ptr[i][j] * B[j][k];
      else if (ABS(sum)) ii = i;
      B[i][k] = sum;
    }
    for (i= mh; i>= ml; i--) {
      sum = B[i][k];
      for (j = mh; j> i; j--)  sum -= ptr[i][j] * B[j][k];
      B[i][k] = sum/ptr[i][i];
    }
  }
}

DMatrix iproduct(const DMatrix & A, const DMatrix & B){
  int * index = new int[1+A.mh-A.ml];
  double d;
  DMatrix  C = A.LUdecomp(index,d);
  DMatrix  D = B;
  C.LUbacksub(D,index);
  delete [] index;
  return D;
}
  
DMatrix DMatrix::inverse() const {
  int * index = new int[1+mh-ml];
  double d;
  DMatrix C = LUdecomp(index,d);
  DMatrix D(ml,mh,1);
  C.LUbacksub(D,index);
  delete [] index;  
  return D;
}
  
double DMatrix::det() const {
  int * index = new int[1+mh-ml];
  double d;
  DMatrix C = LUdecomp(index,d);
  for (int i=ml; i<= mh; i++) d *= C[i][i];
  delete [] index; 
  return d;
}
  
double DMatrix::logdet() const {
  int * index = new int[1+mh-ml];
  double d; 
  double  sum=0;
  DMatrix C = LUdecomp(index,d);
  for (int i=ml; i<= mh; i++) sum += log(ABS(C[i][i]));
  delete [] index; 
  return sum;
}
  

/*  RSMatrix --  easy function definitions    */

void RSMatrix::check() const {
  if (ml != nl || mh != nh) therror(" non-square RSMatrix ");
  int i,j;
  for (i = ml;i<=mh;i++)
    for (j= i+1; j<= mh; j++)
      if (ptr[i][j] != ptr[j][i])
	   therror(" non-symmetric RSMatrix ");
} 


/*  RSMatrix -- hard function definitions   */

/*  this routine diagonalizes a real symmetric matrix, returning a
       matrix of eigenvectors:
          V[i][j] is the ith element of the jth eigenvector
     the DVector d, which should be predefined, returns the 
        eigenvalues.

     The algorithm is based on the subroutine jacobi, from 

       Numerical Recipes, by Press, Teukolsky, Vetterling, and Flannery

             C version,   2nd  ed.  1992

    The matrix is diagonalized so that the 
       final off-diagonal elements are smaller than the original sum
      of absolute values of matrix elements times eps.
                                                       */

DMatrix RSMatrix::diagonalize(DVector & d ,double acc) const {
  check();    /* check that this matrix is indeed RS */
  if ((mh-ml) == 0) {    /*    1 X 1 case    */
    d[ml] = ptr[ml][ml];
    DMatrix C(ml,ml,1.0);
    return C;
  }
  double thresh, t= 0.0;
  double c,s,ta;
  int i,j,iq,ip;
  if ((mh-ml) == 1) {    /*    2 X 2 case    */
    double aa = ptr[ml][ml];
    double bb = ptr[ml][mh];
    double cc = ptr[mh][mh];
    double gg = 100.0*ABS(bb);
    double hh = cc-aa;
    if ((ABS(hh) + gg) == ABS(hh)) {
      t = bb/hh; }
    else {
      double theta = 0.5*hh/bb;
      t = 1.0/(ABS(theta) + sqrt(1.0 + theta*theta));
      if (theta < 0.0)   t = -t ;
    } 
    c = 1.0/sqrt(1 + t*t);
    s = t*c;
    d[ml] = aa - t* bb;
    d[mh] = cc + t* bb;
    DMatrix C(ml,mh,c);  
    C[ml][mh] = s;
    C[mh][ml] = -s;
    return C;
  }     
  DMatrix A = *this;
      /*  compute the sum of the ABS of off-diagonal elements */
   double target = 0.0;
   for (ip= ml;ip<= mh; ip++) {
      for (iq=ip;iq <= mh; iq++)  target +=ABS(A[ip][iq]);
   }
   target *= acc;
   /*   initialize b and d to the diagonal of A */
   DVector b(ml,mh), z(ml,mh);   /*  z is initialized to zero  */
   for (ip = ml; ip <= mh; ip++) b[ip] = A[ip][ip];
   d = b;
   /* initialize v to the identity  */
   DMatrix V(ml,mh,1.0);
   int nrot = 0 ;
   for (i=1; i<= 50; i++) {       /* i = no. of sweeps */
      double sm = 0.0 ;      /* sum of off-diagonal elements  */
      for (ip = ml; ip <= mh-1; ip++) {
         for (iq = ip+1; iq <= mh; iq++)
            sm += ABS(A[ip][iq]);
      }
      if (sm < target)   return V;  
   /* normal return; a is diagonal to accuracy eps  */
      if (i < 4)
         thresh = 0.2 * sm/(SQR(mh-ml));
      else
         thresh = 0.0 ;
      for (ip= ml; ip<=mh-1; ip++) {
         for (iq = ip+1 ; iq <= mh; iq++) {
            double g = 100.0*ABS(A[ip][iq]);
            if (i> 4 &&
               (ABS(d[ip]) + g) ==  ABS(d[ip]) &&
               (ABS(d[iq]) + g) ==  ABS(d[iq])   )
                  A[ip][iq] = 0.0 ;
            else if (ABS(A[ip][iq]) > thresh) {
               double h = d[iq]-d[ip];
               if ( (ABS(h)+g) ==  ABS(h))
                  t = (A[ip][iq])/h ;
               else      {
                  double theta = 0.5*h/(A[ip][iq]);
                  t = 1.0/(ABS(theta) + sqrt(1.0 + theta*theta));
                  if (theta < 0.0)   t = -t ;
               }
               c = 1.0/sqrt(1 + t*t);
               s = t*c;
               ta =  s/(1.0+ c);
               h = t*A[ip][iq];
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               A[ip][iq] = 0.0;
               for (j= ml; j<=ip-1;j++)  {
                  A.ROTATE(s,ta,j,ip,j,iq);
               }
               for (j=ip+1; j<=iq-1;j++) {
                  A.ROTATE(s,ta,ip,j,j,iq);
               }
               for (j= iq+1; j<= mh; j++)  {
                  A.ROTATE(s,ta,ip,j,iq,j);
               }
               for (j=ml; j<=mh;j++)   {
                  V.ROTATE(s,ta,j,ip,j,iq);
               }
               ++(nrot);
            }
         }
      }
      b += z;
      d = b;
      z *= 0;
   }
   therror("Too many interations in routine jacobi");
   return A;  /* never get here */
}

 /*   basic Jacobi rotation:      */

inline  void DMatrix::ROTATE(double s,double ta,
                                int i,int j,int k,int l){
   double g = ptr[i][j];
   double h = ptr[k][l];
   ptr[i][j] = g - s*(h+g*ta);
   ptr[k][l] = h + s*(g-h*ta);
 }




