
#include "CClasses.h"


using namespace std;

const int OFFSET = 1;
  
  /*  the indexing and storage allocation for vectors and matrices follow the 
        prescriptions of Numerical Recipes in C, by 
            Press, Teukolsky, Vetterling, and Flannery  */

CVector::CVector() { 
    therror(" to construct a CVector, use the constructor with parameters");
}

CVector::CVector(int xnl, int xnh){
       init(xnl,xnh);
       zero();
}

CVector::CVector(const CVector & A){
       init(A.nl, A.nh);
       load(A);
}

CVector::CVector(const DVector & A){
       init(A.nl, A.nh);
       loadD(A);
}

CVector::CVector(const IVector & A){
       init(A.nl, A.nh);
       loadI(A);
}

CVector::~CVector(){ 
	destroy();
}

CVector & CVector::operator = (const CVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible CVector");
        }
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = A.ptr[i];
	return * this;
}


int CVector::l() const {return nl;}
int CVector::h() const {return nh;} 

void CVector::fill(int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        n *= 2;
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++) {
                double A = va_arg(a,double);
                double B = va_arg(a,double);
                ptr[i] = A + B*I;
        }
        va_end(a);
}

void CVector::zero() { 
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = 0;
}


CVector CVector::operator - () const {
        CVector A(nl,nh);
        int i;
        for (i = nl; i<= nh; i++) A.ptr[i] = -ptr[i]; 
        return A;
}

CVector & CVector::operator *= (complex  A) {
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] * A;
	return * this;
}

CVector & CVector::operator /= (complex  A) {
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] / A;
	return *this;
}

CVector & CVector::operator += (const CVector & A) {
        if (A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible CVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] +  A.ptr[i];
	return * this;
}
 
CVector & CVector::operator -= (const CVector & A) {
	if (A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible CVectors ");}
	int i;
	for (i = nl; i<=nh; i++) ptr[i] = ptr[i] - A.ptr[i];
	return * this;
}

CVector CVector::adjoint() const {
        CVector A(nl,nh);
	int i;
	for (i = nl; i<=nh; i++) A.ptr[i] = conj(ptr[i]);
	return A;
}

 
void CVector::init(int xnl,int xnh){
      nl = xnl; nh = xnh;
      if ((nh-nl)<0) 
         therror(" the number of components of a CVector should be positive ");
      ptr = new complex[nh-nl+1 + OFFSET];
      if (!ptr) therror("allocation failure in CVector()");
      ptr -= (nl - OFFSET); 
}

void CVector::load(const CVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void CVector::loadD(const DVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void CVector::loadI(const IVector & A) {
	int i;
	for (i=nl; i<= nh; i++)  ptr[i] = A.ptr[i];
}

void CVector::destroy() {
	ptr += (nl-OFFSET);
	delete [] ptr;
}	



std::ostream & operator << (std::ostream & os, const CVector & A){
    int i;
    os << "<";
    for (i=A.nl;i<= A.nh;i++) os << "\t" << A.ptr[i];
    os << ">" << endl;
    return os;
}


int operator ==  (const CVector &A,const CVector &B) {
  if ( A.nl != B.nl  || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i;
  for (i=A.nl; i<= A.nh; i++) {  if (A.ptr[i] != B.ptr[i]) return 0;};
  return 1;
}


int operator !=  (const CVector &A,const CVector &B) {
  return (A == B ? 0 : 1);
}

int eqacc(const CVector & A, const CVector &B, double acc) {
  if ( A.nl != B.nl  || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i;
  for (i=A.nl; i<= A.nh; i++) {
    if (ABS(A.ptr[i] - B.ptr[i]) > acc) return 0; }; 
  return 1;
}

CVector operator + (const CVector & A, const CVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to add incompatible CVectors ");}
  CVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] + B.ptr[i];
  return C; 
}
    
CVector operator - (const CVector & A, const CVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to subtract incompatible CVectors ");}
  CVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] - B.ptr[i];
  return C; 
}

CVector operator * (complex A, const CVector & B){
  CVector C(B.nl,B.nh);
  int i;
  for (i=B.nl; i<= B.nh; i++) C.ptr[i] = A * B.ptr[i];
  return C; 
}

CVector operator * (const CVector & A, complex  B){
  CVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] * B;
  return C; 
}

CVector  operator / (const CVector & A, complex  B){
  CVector C(A.nl,A.nh);
  int i;
  for (i=A.nl; i<= A.nh; i++) C.ptr[i] = A.ptr[i] / B;
  return C; 
}

complex  operator * (const CVector & A, const CVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to dot incompatible CVectors ");}
  complex s = 0;
  int i;
  for (i=A.nl; i<= A.nh;i++) s += A.ptr[i] * B.ptr[i];
  return s;
}

complex norm(const CVector & A, const CVector & B){
  if (A.nl != B.nl || A.nh != B.nh){
    therror(" attempt to dot incompatible CVectors ");}
  complex s = 0;
  int i;
  for (i=A.nl; i<= A.nh;i++) s += conj(A.ptr[i]) * B.ptr[i];
  return s;
}



CMatrix::CMatrix() { 
      therror("to construct a CMatrix, use the constructor with parameters");
}

CMatrix::CMatrix(int xml, int xmh, int xnl, int xnh){
       init(xml,xmh,xnl,xnh);
       zero();
}

CMatrix::CMatrix(const CMatrix & A){
       init(A.ml,A.mh,A.nl, A.nh);
       load(A);
}

CMatrix::CMatrix(const IMatrix & A){
       init(A.ml,A.mh,A.nl, A.nh);
       loadI(A);
}

CMatrix::CMatrix(int xnl, int xnh, complex a) {
       init(xnl,xnh,xnl,xnh);
       int i,j;
       for (i=ml; i<= mh; i++) {
	 for (j=nl;j<=nh; j++) ptr[i][j] = 0;
	 ptr[i][i] =  a;
       }
}

CMatrix::~CMatrix(){ 
	destroy();
}

CMatrix &  CMatrix::operator = (const CMatrix & A) {
	if (A.ml != ml || A.mh != mh ||A.nl != nl || A.nh != nh){
	  therror(" attempt to assign from an incompatible CMatrix");
        }
	int i,j;
	for (i = ml; i<=mh; i++) 
           for (j = nl; j<=nh; j++)  ptr[i][j] = A.ptr[i][j];
	return * this;
}

int CMatrix::rl() const {return ml;}
int CMatrix::rh() const {return mh;}
int CMatrix::cl() const {return nl;}
int CMatrix::ch() const {return nh;} 

void CMatrix::fillRow(int m, int n, ...) {
	if (n != (1+ nh-nl)) therror(" wrong size argument for fillVector ");
        n *= 2;
        va_list a;
        va_start(a,n);
	int i;
        for (i = nl; i<= nh; i++){
           double A = va_arg(a,double);
           double B = va_arg(a,double);
           ptr[m][i] = A + B*I;
        }
        va_end(a);
}

void CMatrix::readRow(int n,const CVector & A){
        if (A.nl != nl || A.nh != nh)
         therror(" attempt to fill a CMatrix with an incompatible CVector");
        int i;
        for (i= nl; i<= nh; i++)  ptr[n][i] = A.ptr[i];
}

void CMatrix::readColumn(int n, const CVector & A){
        if (A.nl != ml || A.nh != mh)
         therror(" attempt to fill a CMatrix with an incompatible CVector");
        int i;
        for (i= ml; i<= mh; i++)  ptr[i][n] = A.ptr[i];
}

void CMatrix::zero() {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = 0;
}

CVector CMatrix::row(int i) const{
	CVector A(nl,nh);
	int j;
        for (j=nl;j<=nh;j++) A[j] = ptr[i][j];
	return A;
}

CVector CMatrix::column(int j) const{
	CVector A(ml,mh);
	int i;
        for (i=ml;i<=mh;i++) A[i] = ptr[i][j];
	return A;
}

CMatrix CMatrix::operator - () const {
        CMatrix A(ml,mh,nl,nh);
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) A.ptr[i][j] = - ptr[i][j];
	return A;

}

CMatrix & CMatrix::operator *= (complex A) {
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) ptr[i][j] = ptr[i][j] * A;
	return * this; 
}

CMatrix & CMatrix::operator /= (complex A) {
	int i,j;
	for (i = ml; i<=mh; i++)
	  for (j=nl;j<=nh;j++) ptr[i][j] = ptr[i][j] / A;
	return *this;
}

CMatrix & CMatrix::operator += (const CMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to add incompatible CMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] +  A.ptr[i][j];
	return * this;
}

CMatrix & CMatrix::operator -= (const CMatrix & A) {
        if (A.ml != ml || A.mh != mh ||  A.nl != nl || A.nh != nh){
	  therror(" attempt to subtract incompatible CMatrices ");}
	int i,j;
	for (i = ml; i<=mh; i++) 
	  for (j= nl;j<= nh; j++) ptr[i][j] = ptr[i][j] -  A.ptr[i][j];
	return * this;
}

void CMatrix::init(int xml, int xmh, int xnl,int xnh){
	ml = xml; mh = xmh; nl = xnl; nh = xnh;
        if ((mh-ml)<0) 
           therror(" the number of rows of a CMatrix should be positive ");
        if ((nh-nl)<0) 
           therror(" the number of columns of a CMatrix should be positive ");
	int nrow = mh-ml+1;
        int ncol = nh-nl+1;
	ptr = new complex * [nrow + OFFSET];
	if (!ptr)  therror ("allocation failure in CMatrix()");
	ptr += OFFSET;
        ptr -= ml;    /*  now ptr has the value s.t. ptr + ml is the correct
                              starting address   */
	ptr[ml] = new complex[nrow*ncol + OFFSET];
	if (!ptr[ml]) therror("allocation failure, stage 2, in CMatrix()");
	ptr[ml] += OFFSET;
	ptr[ml] -= nl;
             /* now evaluate the pointers to rows  */
	for (int i = ml+1; i<= mh; i++) ptr[i] = ptr[i-1] + ncol;
} 

void CMatrix::load(const CMatrix & A) {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = A.ptr[i][j];
}

void CMatrix::loadI(const IMatrix & A) {
	int i,j;
	for (i=ml; i<= mh; i++)
	  for (j=nl; j<= nh; j++)  ptr[i][j] = A.ptr[i][j];
}

void CMatrix::destroy() {
	complex * w = ptr[ml] + nl - OFFSET;
	delete [] w;
	ptr += (ml - OFFSET);
	delete [] ptr;
}

ostream & operator << (ostream & os, const CMatrix & A){
    int i,j;
    for (i=A.ml;i<=A.mh;i++) {
       os << endl <<  "<";
       for (j=A.nl;j<= A.nh;j++) os << "\t" << A.ptr[i][j];
    os << "\t>" << endl;
    }
    os << endl;
    return os;
}

int operator ==  (const CMatrix &A,const CMatrix &B) {
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) if (A.ptr[i][j] != B.ptr[i][j]) return 0;
  return 1;
}

int operator !=  (const CMatrix &A,const CMatrix &B) {
  return (A == B ? 0 : 1);
}


int eqacc(const CMatrix &A,const CMatrix &B, double acc) {
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) return 0;
  if ( A.ptr == B.ptr) return 1;
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) 
       if (ABS(A.ptr[i][j] -B.ptr[i][j])> acc) return 0;
  return 1;
}

CMatrix operator + (const CMatrix & A, const CMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to add incompatible CMatrices ");
  CMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j]+ B.ptr[i][j];
  return C; 
}
    
CMatrix operator - (const CMatrix & A, const CMatrix & B){
  if (A.ml != B.ml || A.mh != B.mh ||  A.nl != B.nl || A.nh != B.nh) 
    therror(" attempt to subtract incompatible CMatrices ");
  CMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] - B.ptr[i][j];
  return C;
} 

CMatrix operator * (complex A, const CMatrix & B){
  CMatrix C(B.ml,B.mh,B.nl,B.nh);
  int i,j;
  for (i=B.ml; i<= B.mh; i++)  
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A * B.ptr[i][j];
  return C; 
}

CMatrix operator * (const CMatrix & A, complex B){
  CMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] * B;
  return C; 
}

CMatrix operator / (const CMatrix & A, complex B){
  CMatrix C(A.ml,A.mh,A.nl,A.nh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++)  
    for (j=A.nl;j<=A.nh; j++) C.ptr[i][j] = A.ptr[i][j] / B;
  return C; 
}

CVector operator * (const CVector & A, const CMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible CVector and CMatrix ");
  CVector C(B.nl,B.nh);
  int i,j;
  for (j=B.nl; j<= B.nh; j++) {
    C.ptr[j] = 0;
    for (i=B.ml;i<=B.mh; i++) C.ptr[j] += A.ptr[i] * B.ptr[i][j];
  }
  return C; 
}

CVector operator * (const CMatrix & A, const CVector & B){
  if (A.nl != B.nl || A.nh != B.nh)
    therror(" attempt to multiply incompatible CMatrix and CVector ");
  CVector C(A.ml,A.mh);
  int i,j;
  for (i=A.ml; i<= A.mh; i++) {
    C.ptr[i] = 0;
    for (j=B.nl;j<=B.nh; j++) C.ptr[i] += A.ptr[i][j] * B.ptr[j];
  }
  return C; 
}

CMatrix operator * (const CMatrix & A, const CMatrix & B){
  if (A.nl != B.ml || A.nh != B.mh)
    therror(" attempt to multiply incompatible DMatrices ");
  CMatrix C(A.ml,A.mh,B.nl,B.nh);
  int i,j,k;
  for (i=A.ml; i<= A.mh; i++) {
    for (j=B.nl; j<= B.nh; j++) {
      C.ptr[i][j] = 0;
      for (k=A.nl;k<=A.nh;k++) C.ptr[i][j] += A.ptr[i][k] * B.ptr[k][j];
    }
  }
  return C; 
}

CMatrix dyad(const CVector & A, const CVector & B){
  CMatrix  C(A.nl,A.nh,B.nl,B.nh);
  int i,j;
  for (i=A.nl; i<= A.nh; i++)
    for (j=B.nl;j<=B.nh; j++) C.ptr[i][j] = A.ptr[i] * B.ptr[j];
  return C; 
}


/*  methods for inverse,  determinant, etc.     */

CMatrix CMatrix::transpose() const {
  CMatrix A(nl,nh,ml,mh);
  int i,j;
  for (i=ml; i<= mh; i++)
    for (j=nl; j<= nh; j++) A.ptr[j][i] = ptr[i][j];
  return A;
}

CMatrix CMatrix::adjoint() const {
  CMatrix A(nl,nh,ml,mh);
  int i,j;
  for (i=ml; i<= mh; i++)
    for (j=nl; j<= nh; j++) A.ptr[j][i] = conj(ptr[i][j]);
  return A;
}

complex CMatrix::trace() const {
  if (mh != nh || ml != nl)
         therror(" attempt to take the trace of an asymmetric CMatrix ");
  int i;   complex tr = 0;
  for (i = nl; i <= nh; i++)  tr += ptr[i][i];
  return tr;
}

void CMatrix::rowint(int i, int j) {
  int k;
  complex el;
  for (k= nl; k<= nh; k++) {
    el = ptr[i][k];
    ptr[i][k] = ptr[j][k];
    ptr[j][k] = el;
  }
}

void CMatrix::colint(int i, int j) {
  int k;
  complex el;
  for (k= ml; k<= mh; k++) {
    el = ptr[k][i];
    ptr[k][i] = ptr[k][j];
    ptr[k][j] = el;
  }
}

double CMatrix::maxRow(int i) const { 
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

double CMatrix::maxCol(int i) const { 
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


CMatrix CMatrix::LUdecomp(int *  index, complex & d) const {
  if (mh != nh || ml != nl)
         therror(" attempt to invert an asymmetric CMatrix ");
  int i,j,k, imax;
  complex sum, t;
  double big,temp;
  CMatrix C = *this;   /* copy the matrix into C, and work on C */
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
 
void CMatrix::LUbacksub(CMatrix & B,int * index) {
  /*   note that B is destroyed and replaced in this operation;
           on the other hand, the inverted matrix is preserved    */
  if (nh != B.mh || nl != B.ml)
        therror(" attempt to multiply incompatible CMatrices in LUbacksub ");
  int i,j,k, ii, ip;
  complex sum;
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

CMatrix iproduct(const CMatrix & A, const CMatrix & B){
  int * index = new int[1+A.mh-A.ml];
  complex d;
  CMatrix  C = A.LUdecomp(index,d);
  CMatrix  D = B;
  C.LUbacksub(D,index);
  delete [] index;
  return D;
}
  
CMatrix CMatrix::inverse() const {
  int * index = new int[1+mh-ml];
  complex d;
  CMatrix C = LUdecomp(index,d);
  CMatrix D(ml,mh,1);
  C.LUbacksub(D,index);
  delete [] index;  
  return D;
}
  
complex CMatrix::det() const {
  int * index = new int[1+mh-ml];
  complex d;
  CMatrix C = LUdecomp(index,d);
  for (int i=ml; i<= mh; i++) d *= C[i][i];
  delete [] index; 
  return d;
}
  
double CMatrix::logdet() const {
  int * index = new int[1+mh-ml];
  complex d; 
  double  sum=0;
  CMatrix C = LUdecomp(index,d);
  for (int i=ml; i<= mh; i++) sum += log(ABS(C[i][i]));
  delete [] index; 
  return sum;
}
  


void HMatrix::check() const{
  if (ml != nl || mh != nh) therror(" non-square HMatrix ");
  int i,j;
  for (i = ml;i<=mh;i++)
    for (j= i; j<= mh; j++)
      if (ptr[i][j] != conj(ptr[j][i]))
	   therror(" non-Hermitian HMatrix ");
} 

/*  this routine diagonalizes a Hermitian matrix, returning a
       matrix of eigenvectors:
          V[i][j] is the ith element of the jth eigenvector
     the DVector d, which should be predefined, returns the 
        eigenvalues.

     The algorithm is based on the subroutine jacobi, from 

       Numerical Recipes, by Press, Teukolsky, Vetterling, and Flannery

             C version,   2nd  ed.  1992

    carried into complex analysis in a suitable way.
    The matrix is diagonalized so that the 
       final off-diagonal elements are smaller than the original sum
      of absolute values of matrix elements times eps.
                                                       */

CMatrix HMatrix::diagonalize(DVector & d ,double acc) const {
  check();    /* check that this matrix is indeed H */
  if ((mh-ml) == 0) {    /*    1 X 1 case    */
    d[ml] = real(ptr[ml][ml]);
    CMatrix C(ml,ml,1.0);
    return C;
  }
  complex t= 0.0;
  complex s,ta;
  double c,thresh;
  int i,j,iq,ip;
  if ((mh-ml) == 1) {    /*    2 X 2 case    */
    double aa = real(ptr[ml][ml]);
    complex bb = ptr[ml][mh];
    double cc = real(ptr[mh][mh]);
    double gg = 100.0*ABS(bb);
    double hh = cc-aa;
    if ((ABS(hh) + gg) == ABS(hh)) {
      t = bb/hh; }
    else {
      t = 2.0* bb/(ABS(hh) + sqrt(4.0 * norm(bb) + hh*hh));
      if (hh < 0.0)   t = -t ;
    } 
    c = 1.0/sqrt(1 + norm(t));
    s = t*c;
    d[ml] = aa - real(conj(t)* bb);
    d[mh] = cc + real(conj(t)* bb);
    CMatrix C(ml,mh,c);  
    C[ml][mh] = s;
    C[mh][ml] = -conj(s);
    return C;
  }     
  CMatrix A = *this;
      /*  compute the sum of the ABS of off-diagonal elements */
   double target = 0.0;
   for (ip= ml;ip<= mh; ip++) {
      for (iq=ip;iq <= mh; iq++)  target +=ABS(A[ip][iq]);
   }
   target *= acc;
   /*   initialize b and d to the diagonal of A */
   DVector b(ml,mh), z(ml,mh);   /*  z is initialized to zero  */
   for (ip = ml; ip <= mh; ip++) b[ip] = real(A[ip][ip]);
   d = b;
   /* initialize V to the identity  */
   CMatrix V(ml,mh,1.0);
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
	       complex B = A[ip][iq];
               if ( (ABS(h)+g) ==  ABS(h))
                  t = B/h ;
               else      {
		  t = 2.0* B/(ABS(h) + sqrt(4.0 * norm(B) + h*h));
		  if (h < 0.0)   t = -t ;
               }
               c = 1.0/sqrt(1 + norm(t));
               s = t*c;
               ta =  s/(1.0+ c);
               h = real(conj(t)*B);
               z[ip] -= h;
               z[iq] += h;
               d[ip] -= h;
               d[iq] += h;
               A[ip][iq] = 0.0;
               for (j= ml; j<=ip-1;j++)  {
                  A.ROTATE(s,ta,j,ip,j,iq);
               }
               for (j=ip+1; j<=iq-1;j++) {
                  A.ROTATE2(s,ta,ip,j,j,iq);
               }
               for (j= iq+1; j<= mh; j++)  {
                  A.ROTATE3(s,ta,ip,j,iq,j);
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

 /*   basic Jacobi rotations:      */

inline  void CMatrix::ROTATE(complex s,complex  ta,
                                int i,int j,int k,int l){
   complex g = ptr[i][j];
   complex h = ptr[k][l];
   ptr[i][j] = g - conj(s)*(h+g*ta);
   ptr[k][l] = h + s*(g-h*conj(ta));
 }

inline  void CMatrix::ROTATE2(complex s,complex  ta,
                                int i,int j,int k,int l){
   complex g = conj(ptr[i][j]);
   complex h = ptr[k][l];
   ptr[i][j] = conj(g - conj(s)*(h+g*ta));
   ptr[k][l] = h + s*(g-h*conj(ta));
 }

inline  void CMatrix::ROTATE3(complex s,complex  ta,
                                int i,int j,int k,int l){
   complex g = conj(ptr[i][j]);
   complex h = conj(ptr[k][l]);
   ptr[i][j] = conj(g - conj(s)*(h+g*ta));
   ptr[k][l] = conj(h + s*(g-h*conj(ta)));
 }

